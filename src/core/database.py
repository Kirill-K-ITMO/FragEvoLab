"""database - FragEvoLab module."""

import sqlite3
from pathlib import Path
from typing import List, Dict, Any, Optional
from datetime import datetime
import json

DATABASE_FILE = "output/experiments.db"

class DockingDatabase:
    """
    Manage the docking database. Save,
    update, and retrieve experiments, ligands, receptors, and docking results.
    """
    def __init__(self, db_file: Optional[Path] = None):
        self.db_file = Path(db_file) if db_file else Path(DATABASE_FILE)
        self.connection = sqlite3.connect(self.db_file)
        self.connection.row_factory = sqlite3.Row
        self._init_schema()

    def _init_schema(self):
        """
        Initialize tables if they do not already exist.
        """
        with self.connection:
            self.connection.execute("""
            CREATE TABLE IF NOT EXISTS experiments (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT,
                started_at TEXT,
                finished_at TEXT,
                user_params TEXT, -- JSON string containing run parameters
                description TEXT
            )""")
            self.connection.execute("""
            CREATE TABLE IF NOT EXISTS ligands (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT,
                smiles TEXT,
                formula TEXT,
                mw REAL,
                logp REAL,
                tpsa REAL,
                hbd INTEGER,
                hba INTEGER,
                sa_score REAL
            )""")
            self.connection.execute("""
            CREATE TABLE IF NOT EXISTS proteins (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT,
                pdbqt_file TEXT
            )""")
            self.connection.execute("""
            CREATE TABLE IF NOT EXISTS dockings (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                experiment_id INTEGER,
                ligand_id INTEGER,
                protein_id INTEGER,
                delta_g REAL,
                ki REAL,
                pki REAL,
                rank INTEGER,
                method TEXT,
                exhaustiveness INTEGER,
                extra_json TEXT,
                FOREIGN KEY (experiment_id) REFERENCES experiments(id),
                FOREIGN KEY (ligand_id) REFERENCES ligands(id),
                FOREIGN KEY (protein_id) REFERENCES proteins(id)
            )""")

    def add_experiment(self, name:str, user_params:dict, description:str="") -> int:
        """
        Add an experiment and return its ID.
        """
        with self.connection:
            cursor = self.connection.execute("""
                INSERT INTO experiments (name, started_at, user_params, description)
                VALUES (?, ?, ?, ?)
            """, (
                name,
                datetime.now().isoformat(),
                json.dumps(user_params, ensure_ascii=False),
                description
            ))
            return cursor.lastrowid

    def finish_experiment(self, experiment_id:int):
        """
        Record the experiment end time.
        """
        with self.connection:
            self.connection.execute(
                "UPDATE experiments SET finished_at=? WHERE id=?",
                (datetime.now().isoformat(), experiment_id)
            )

    def add_or_get_ligand(self, ligand:dict) -> int:
        """
        Add a ligand if absent; otherwise return the existing ID.
        ligand: {'name':..., 'smiles':..., 'formula':..., ...}
        """
        # Look up by name and SMILES
        res = self.connection.execute(
            "SELECT id FROM ligands WHERE name=? AND smiles=?",
            (ligand['name'], ligand['smiles'])
        ).fetchone()
        if res:
            return res['id']
        # Insert
        with self.connection:
            cursor = self.connection.execute("""
                INSERT INTO ligands (name, smiles, formula, mw, logp, tpsa, hbd, hba, sa_score)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                ligand['name'], ligand['smiles'], ligand.get('formula'),
                ligand.get('mw'), ligand.get('logp'), ligand.get('tpsa'),
                ligand.get('hbd'), ligand.get('hba'), ligand.get('sa_score')
            ))
            return cursor.lastrowid

    def add_or_get_protein(self, protein:dict) -> int:
        """
        Add a receptor or return its ID (unique by name and pdbqt_file).
        protein: {'name':..., 'pdbqt_file':...}
        """
        res = self.connection.execute(
            "SELECT id FROM proteins WHERE name=? AND pdbqt_file=?",
            (protein['name'], protein['pdbqt_file'])
        ).fetchone()
        if res:
            return res['id']
        with self.connection:
            cursor = self.connection.execute(
                "INSERT INTO proteins (name, pdbqt_file) VALUES (?, ?)",
                (protein['name'], protein['pdbqt_file'])
            )
            return cursor.lastrowid

    def add_docking_result(self, experiment_id:int, ligand_id:int, protein_id:int, result:dict):
        """
        Save a docking result for a ligand/receptor/experiment combination.
        result = {'delta_g':..., 'ki':..., 'pki':..., 'rank':..., 'method':..., 'exhaustiveness':..., ...}
        """
        with self.connection:
            self.connection.execute("""
                INSERT INTO dockings (
                    experiment_id, ligand_id, protein_id,
                    delta_g, ki, pki, rank, method, exhaustiveness, extra_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                experiment_id,
                ligand_id,
                protein_id,
                result.get('delta_g'),
                result.get('ki'),
                result.get('pki'),
                result.get('rank', 1),
                result.get('method', 'vina'),
                result.get('exhaustiveness', 8),
                json.dumps(result, ensure_ascii=False)
            ))
    
    def get_experiment_list(self) -> List[Dict[str, Any]]:
        """
        List all experiments.
        """
        res = self.connection.execute(
            "SELECT * FROM experiments ORDER BY started_at DESC"
        ).fetchall()
        return [dict(r) for r in res]

    def get_experiment_results(self, experiment_id: int) -> List[Dict[str, Any]]:
        """
        Retrieve all results for experiment_id.
        """
        sql = """
        SELECT d.*, 
               l.name as ligand_name, l.smiles, l.formula, l.mw, l.logp, l.tpsa, l.hbd, l.hba, l.sa_score,
               p.name as protein_name, p.pdbqt_file
        FROM dockings d
            JOIN ligands l ON d.ligand_id = l.id
            JOIN proteins p ON d.protein_id = p.id
        WHERE d.experiment_id = ?
        ORDER BY d.rank, l.name, p.name
        """
        res = self.connection.execute(sql, (experiment_id,)).fetchall()
        return [dict(r) for r in res]

    def get_full_history(self) -> List[Dict[str, Any]]:
        """
        Retrieve the full history (all docking calculations in all experiments).
        """
        sql = """
        SELECT e.id as experiment_id, e.name as experiment_name, e.started_at, e.finished_at, 
               d.*, 
               l.name as ligand_name, l.smiles, l.formula, l.mw, l.logp, l.tpsa, l.hbd, l.hba, l.sa_score,
               p.name as protein_name, p.pdbqt_file
        FROM experiments e
            JOIN dockings d ON d.experiment_id = e.id
            JOIN ligands l ON d.ligand_id = l.id
            JOIN proteins p ON d.protein_id = p.id
        ORDER BY e.id DESC, d.rank, l.name, p.name
        """
        res = self.connection.execute(sql).fetchall()
        return [dict(r) for r in res]

    def export_for_report(self) -> Dict:
        """
        Export all data for generating db-report.html.
        Include SMILES for structure depiction.
        
        Returns:
            Dictionary containing results, ligands (including SMILES), and receptors
        """
        rows = self.get_full_history()
        results = [{**row, "ligand_smiles": row["smiles"],
                    "timestamp": row["started_at"]} for row in rows]
        return {"results": results,
                "total_experiments": self.connection.execute("SELECT COUNT(*) FROM experiments").fetchone()[0],
                "total_ligands": self.connection.execute("SELECT COUNT(*) FROM ligands").fetchone()[0],
                "total_proteins": self.connection.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]}

    def close(self):
        self.connection.close()

