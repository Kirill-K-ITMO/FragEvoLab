"""Content-addressed, transactional cache for a single-host docking workflow.

Only successful validated computations are reusable. Names and run IDs are
labels, never cache keys. SQLite contains both metadata and artifact bytes so
an interrupted commit cannot expose half of a result. Keep this DB on local disk.
"""
from __future__ import annotations

from contextlib import closing
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import tempfile
from typing import Callable
from filelock import FileLock

def json_text(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False)

def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def fingerprint(value) -> str:
    return digest(json_text(value).encode("utf-8"))

def atomic_write(path: Path, data: bytes) -> None:
    """Never make partial artifacts visible, including on a cache hit."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix=".writing-", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp, path)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)

def atomic_json(path: Path, value) -> None:
    atomic_write(path, json_text(value).encode("utf-8"))

def molecule_identity(smiles: str) -> tuple[str, str]:
    from rdkit import Chem, rdBase
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("Invalid or empty SMILES")
    if mol.GetStereoGroups():
        raise ValueError("Enhanced stereochemistry (CXSMILES stereo groups) is unsupported")
    # Atom-map labels are bookkeeping, not chemical identity. Do not merge
    # stereoisomers, protonation states, isotopes or tautomers automatically.
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    return canonical, fingerprint({"canonical_smiles": canonical,
                                  "rdkit": rdBase.rdkitVersion, "policy": 2})

def validate_pdbqt(data: bytes) -> None:
    """Structural sanity check, not a chemical/synthetic feasibility filter."""
    lines = data.decode("utf-8").splitlines()
    atoms = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    if not atoms:
        raise ValueError("PDBQT contains no atoms")
    for line in atoms:
        try:
            numbers = [float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]
            numbers.append(float(line[70:76]))
            if not all(math.isfinite(x) for x in numbers) or not line[77:].strip():
                raise ValueError()
        except (ValueError, IndexError):
            raise ValueError("Malformed PDBQT coordinates, charge or atom type") from None

@dataclass(frozen=True)
class CacheResult:
    key: str
    payload: dict
    artifact: bytes
    hit: bool

class CalculationCache:
    SCHEMA = 1

    def __init__(self, root: Path, lock_timeout: float = 86400):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.locks = self.root / "locks"
        self.locks.mkdir(exist_ok=True)
        self.database = self.root / "cache.sqlite3"
        self.lock_timeout = lock_timeout
        with FileLock(str(self.locks / "schema.lock"), timeout=60):
            with closing(self._connect()) as db, db:
                version = db.execute("PRAGMA user_version").fetchone()[0]
                if version not in (0, self.SCHEMA):
                    raise RuntimeError(f"Unsupported cache schema: {version}")
                db.execute("PRAGMA journal_mode=WAL")
                db.executescript("""
                    CREATE TABLE IF NOT EXISTS entries (
                        key TEXT PRIMARY KEY, kind TEXT NOT NULL, spec TEXT NOT NULL,
                        payload TEXT NOT NULL, payload_sha TEXT NOT NULL,
                        artifact BLOB NOT NULL, artifact_sha TEXT NOT NULL,
                        created TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP);
                    CREATE TABLE IF NOT EXISTS attempts (
                        id INTEGER PRIMARY KEY, key TEXT NOT NULL, kind TEXT NOT NULL,
                        status TEXT NOT NULL, error TEXT,
                        started TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP, finished TEXT);
                    CREATE INDEX IF NOT EXISTS attempts_key ON attempts(key);
                """)
                db.execute(f"PRAGMA user_version={self.SCHEMA}")

    def _connect(self):
        db = sqlite3.connect(self.database, timeout=60)
        db.execute("PRAGMA busy_timeout=60000")
        db.execute("PRAGMA synchronous=FULL")
        return db

    @staticmethod
    def key_for(kind: str, spec: dict) -> str:
        return fingerprint({"schema": 1, "kind": kind, "spec": spec})

    def get(self, kind: str, spec: dict) -> CacheResult | None:
        key = self.key_for(kind, spec)
        with closing(self._connect()) as db:
            row = db.execute("SELECT spec,payload,payload_sha,artifact,artifact_sha "
                             "FROM entries WHERE key=? AND kind=?", (key, kind)).fetchone()
        if row is None:
            return None
        stored_spec, payload, payload_sha, artifact, artifact_sha = row
        if (stored_spec != json_text(spec) or digest(payload.encode()) != payload_sha
                or digest(artifact) != artifact_sha):
            return None  # Corrupt entry is a miss; next successful compute repairs it.
        try:
            value = json.loads(payload)
        except (ValueError, TypeError):
            return None
        if not isinstance(value, dict):
            return None
        return CacheResult(key, value, artifact, True)

    def get_or_compute(self, kind: str, spec: dict,
                       compute: Callable[[], tuple[dict, bytes]], *,
                       force: bool = False) -> CacheResult:
        key = self.key_for(kind, spec)
        # Different keys compute in parallel. Re-check after acquiring this lock
        # so simultaneous requests for the same calculation do not duplicate work.
        with FileLock(str(self.locks / f"{key}.lock"), timeout=self.lock_timeout):
            if not force:
                previous = self.get(kind, spec)
                if previous is not None:
                    return previous
            with closing(self._connect()) as db, db:
                # The OS released a dead worker's key lock. Its incomplete attempt
                # is history, not a result and not a reason to suppress a retry.
                db.execute("UPDATE attempts SET status='interrupted', finished=CURRENT_TIMESTAMP "
                           "WHERE key=? AND status='running'", (key,))
                attempt = db.execute("INSERT INTO attempts(key,kind,status) VALUES (?,?,'running')",
                                     (key, kind)).lastrowid
            try:
                payload, artifact = compute()
                encoded = json_text(payload)
                if not isinstance(payload, dict) or not isinstance(artifact, bytes) or not artifact:
                    raise ValueError("A successful cache entry requires metadata and artifact bytes")
                with closing(self._connect()) as db, db:
                    db.execute("INSERT OR REPLACE INTO entries "
                               "(key,kind,spec,payload,payload_sha,artifact,artifact_sha) VALUES (?,?,?,?,?,?,?)",
                               (key, kind, json_text(spec), encoded, digest(encoded.encode()),
                                artifact, digest(artifact)))
                    db.execute("UPDATE attempts SET status='success',finished=CURRENT_TIMESTAMP WHERE id=?",
                               (attempt,))
            except BaseException as error:
                with closing(self._connect()) as db, db:
                    db.execute("UPDATE attempts SET status='failed',error=?,finished=CURRENT_TIMESTAMP WHERE id=?",
                               (str(error)[:4000], attempt))
                raise
            return CacheResult(key, payload, artifact, False)

    def clear(self, kind: str | None = None) -> int:
        """Clear reusable entries only. Call after workers stop; never delete runs."""
        with closing(self._connect()) as db, db:
            cursor = (db.execute("DELETE FROM entries WHERE kind=?", (kind,)) if kind
                      else db.execute("DELETE FROM entries"))
            return cursor.rowcount

    def statistics(self) -> dict:
        with closing(self._connect()) as db:
            return {kind: {"entries": count, "artifact_bytes": size}
                    for kind, count, size in db.execute(
                        "SELECT kind,count(*),sum(length(artifact)) FROM entries GROUP BY kind")}
