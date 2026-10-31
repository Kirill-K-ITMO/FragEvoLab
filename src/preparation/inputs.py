"""One CSV schema for preparation, menus, reports and evolution seeds."""
import csv
from collections import Counter
from pathlib import Path

def read_ligand_table(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as stream:
        sample=stream.read(8192);stream.seek(0)
        try:dialect=csv.Sniffer().sniff(sample,delimiters=",;\t")
        except csv.Error:dialect=csv.excel
        reader=csv.DictReader(stream,dialect=dialect)
        if not reader.fieldnames:raise ValueError("Empty CSV")
        fields=[f.strip() for f in reader.fieldnames]
        if len(fields)!=len(set(fields)):raise ValueError("Duplicate CSV columns")
        reader.fieldnames=fields
        if "smiles" not in fields:raise ValueError("CSV requires a smiles column")
        if 'название' in fields and 'name' in fields:raise ValueError("Use either название or name, not both")
        name_field='название' if 'название' in fields else 'name' if 'name' in fields else None
        rows=[]
        for i,row in enumerate(reader,1):
            if None in row:raise ValueError(f"Extra columns on CSV row {i+1}")
            name=(row.get(name_field) or '').strip() if name_field else f'ligand_{i:06d}'
            row["name"] = name
            row["smiles"] = (row.get("smiles") or "").strip()
            if name_field == "название": row.pop("название", None)
            rows.append(row)
        fields = ["name" if f == "название" else f for f in fields]
        if "name" not in fields: fields.insert(0, "name")
        return fields, rows

def read_ligand_rows(path):
    _, rows = read_ligand_table(path)
    return [(row["name"], row["smiles"]) for row in rows]

def unambiguous_rows(rows):
    counts=Counter(name.casefold() for name,smiles in rows)
    return [(name,smiles) for name,smiles in rows if name and smiles and counts[name.casefold()]==1]
