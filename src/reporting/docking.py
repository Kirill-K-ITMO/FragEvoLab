from src.version import APP_NAME, __version__
from src.reporting.branding import brand_html

AFFINITY_NOTE = 'Estimated Ki/pKi from docking score, not experimental affinity. Standard FragEvoLab workflow: 298.15 K (25 C). Ki is in mol/L; pKi is dimensionless. Imported results or custom API temperatures must be verified separately.'


import csv
import base64
from io import BytesIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime
import logging

try:
    from rdkit import Chem
    from rdkit.Chem import Draw, Descriptors, Crippen
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit is not installed. Descriptors will not be calculated.")

try:
    from rdkit.Chem import RDConfig
    import os
    import sys
    sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
    import sascorer
    SA_AVAILABLE = True
except:
    SA_AVAILABLE = False
    print("Warning: SA Score is unavailable.")

from jinja2 import Template

class ReportGenerator:
    """
    Generate molecular docking reports.
    """
    
    def __init__(self, project_root: Optional[Path] = None):
        """
        Initialize the report generator.
        
        Args:
            project_root: Project root directory
        """
        self.project_root = project_root or Path.cwd()
        
        self.data_dir = self.project_root / "data"
        self.reports_dir = self.project_root / "output" / "reports"
        self.cache_dir = self.project_root / "output" / "cache"
        
        self.reports_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger = logging.getLogger(__name__)
        
        self.mol_img_size = (200, 200)
    
    def calculate_molecular_descriptors(self, smiles: str) -> Dict:
        """Calculate molecular descriptors from SMILES."""
        descriptors = {
            'formula': 'N/A',
            'mw': None,
            'logp': None,
            'tpsa': None,
            'hbd': None,
            'hba': None,
            'sa_score': None
        }
        
        if not RDKIT_AVAILABLE:
            return descriptors
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return descriptors
            
            descriptors['formula'] = Chem.rdMolDescriptors.CalcMolFormula(mol)
            descriptors['mw'] = Descriptors.MolWt(mol)
            descriptors['logp'] = Crippen.MolLogP(mol)
            descriptors['tpsa'] = Descriptors.TPSA(mol)
            descriptors['hbd'] = Descriptors.NumHDonors(mol)
            descriptors['hba'] = Descriptors.NumHAcceptors(mol)
            
            if SA_AVAILABLE:
                descriptors['sa_score'] = sascorer.calculateScore(mol)
        
        except Exception as e:
            self.logger.error(f"Failed to calculate descriptors for {smiles}: {e}")
        
        return descriptors
    
    def load_smiles_data(self, csv_path: Optional[Path] = None) -> Dict[str, str]:
        """Load SMILES from a CSV file."""
        if csv_path is None:
            csv_path = self.data_dir / "ligands.csv"
        
        smiles_dict = {}
        
        if not csv_path.exists():
            self.logger.warning(f"CSV file not found: {csv_path}")
            return smiles_dict
        
        try:
            from src.preparation.inputs import read_ligand_rows, unambiguous_rows
            smiles_dict=dict(unambiguous_rows(read_ligand_rows(csv_path)))

            self.logger.info(f"Loaded {len(smiles_dict)} SMILES from {csv_path}")
        
        except Exception as e:
            self.logger.error(f"Failed to read CSV: {e}")
        
        return smiles_dict
    
    def smiles_to_image_base64(self, smiles: str, size: Tuple[int, int] = None) -> Optional[str]:
        """Convert SMILES to a base64-encoded image for HTML."""
        if not RDKIT_AVAILABLE:
            return None
        
        if size is None:
            size = self.mol_img_size
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            img = Draw.MolToImage(mol, size=size)
            
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_base64 = base64.b64encode(buffered.getvalue()).decode('utf-8')
            
            return f"data:image/png;base64,{img_base64}"
        
        except Exception as e:
            self.logger.error(f"Failed to generate an image for {smiles}: {e}")
            return None
    
    def prepare_report_data(self, 
                           enriched_results: Dict,
                           smiles_dict: Dict[str, str]) -> Tuple[List[Dict], List[str]]:
        """Prepare report data with molecular descriptors."""
        report_data = []
        
        all_proteins = set()
        for ligand_data in enriched_results.values():
            all_proteins.update(ligand_data.keys())
        
        proteins_sorted = sorted(all_proteins)
        
        for ligand_name, proteins_data in enriched_results.items():
            smiles = next((row['smiles'] for row in proteins_data.values() if row.get('smiles')), smiles_dict.get(ligand_name, 'N/A'))
            
            descriptors = self.calculate_molecular_descriptors(smiles)
            
            row = {
                'name': ligand_name,
                'smiles': smiles,
                'structure_img': None,
                'formula': descriptors['formula'],
                'mw': descriptors['mw'],
                'logp': descriptors['logp'],
                'tpsa': descriptors['tpsa'],
                'hbd': descriptors['hbd'],
                'hba': descriptors['hba'],
                'sa_score': descriptors['sa_score'],
                'proteins': {}
            }
            
            if smiles != 'N/A':
                row['structure_img'] = self.smiles_to_image_base64(smiles)
            
            for protein in proteins_sorted:
                if protein in proteins_data:
                    data = proteins_data[protein]
                    row['proteins'][protein] = {
                        'delta_g': data.get('affinity', data.get('delta_g')),
                        'ki': data.get('ki'),
                        'ki_formatted': data.get('ki_formatted', 'N/A'),
                        'pki': data.get('pki')
                    }
                else:
                    row['proteins'][protein] = {
                        'delta_g': None,
                        'ki': None,
                        'ki_formatted': 'N/A',
                        'pki': None
                    }
            
            import math
            for field in ('mw', 'logp', 'tpsa', 'hbd', 'hba', 'sa_score'):
                if row[field] is not None and not math.isfinite(row[field]):
                    row[field] = None
            report_data.append(row)
        
        return report_data, proteins_sorted
    
    def generate_html_report(self,
                            enriched_results: Dict,
                            smiles_dict: Optional[Dict[str, str]] = None,
                            output_file: Optional[Path] = None,
                            title: str = "Docking results",
                            exhaustiveness: int = None) -> Path:
        """
        Generate an HTML results report.
        Includes theme switching and filters.
        """
        self.logger.info("Generating HTML report...")
        
        if smiles_dict is None:
            smiles_dict = self.load_smiles_data()
        
        report_data, proteins = self.prepare_report_data(enriched_results, smiles_dict)
        
        if output_file is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = self.reports_dir / f"docking_report_{timestamp}.html"
        
        # Complete HTML template
        html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{ title }}</title>
    <style>
        :root {
            --bg-primary: #ffffff; --bg-secondary: #f8f9fa; --bg-card: #ffffff;
            --text-primary: #1a1a1a; --text-secondary: #6c757d; --border-color: #e9ecef;
            --accent-color: #2d5f3f; --accent-hover: #1f4029; --accent-light: rgba(45, 95, 63, 0.1);
            --shadow: rgba(0,0,0,0.1);
        }
        [data-theme="dark"] {
            --bg-primary: #1a1a1a; --bg-secondary: #2d2d2d; --bg-card: #252525;
            --text-primary: #e9ecef; --text-secondary: #adb5bd; --border-color: #444;
            --accent-color: #4ade80; --accent-hover: #22c55e; --accent-light: rgba(74, 222, 128, 0.15);
            --shadow: rgba(0,0,0,0.5);
        }
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background: var(--bg-primary); color: var(--text-primary); padding: 20px; min-height: 100vh; transition: background 0.3s, color 0.3s; }
        .container { max-width: 98%; margin: 0 auto; background: var(--bg-card); border-radius: 15px; box-shadow: 0 20px 60px var(--shadow); overflow: hidden; }
        .header { background: linear-gradient(135deg, var(--accent-color) 0%, var(--accent-hover) 100%); color: white; padding: 30px; text-align: center; position: relative; }
        .header h1 { font-size: 2.5em; margin-bottom: 10px; text-shadow: 2px 2px 4px rgba(0,0,0,0.2); }
        .header p { font-size: 1.1em; opacity: 0.9; }
        .theme-toggle { position: absolute; top: 20px; right: 20px; background: rgba(255,255,255,0.2); border: none; padding: 10px 20px; border-radius: 20px; cursor: pointer; color: white; font-weight: bold; transition: all 0.3s; }
        .theme-toggle:hover { background: rgba(255,255,255,0.3); transform: scale(1.05); }
        .stats { display: flex; justify-content: space-around; padding: 20px; background: var(--bg-secondary); border-bottom: 2px solid var(--border-color); }
        .stat-box { text-align: center; padding: 15px; }
        .stat-box .number { font-size: 2em; font-weight: bold; color: var(--accent-color); }
        .stat-box .label { color: var(--text-secondary); margin-top: 5px; }
        .controls { padding: 20px; background: var(--bg-secondary); border-bottom: 2px solid var(--border-color); display: flex; flex-wrap: wrap; gap: 15px; align-items: center; }
        .control-group { display: flex; align-items: center; gap: 10px; }
        .controls label { font-weight: bold; color: var(--text-primary); }
        .controls select, .controls input { padding: 8px 12px; border-radius: 5px; border: 1px solid var(--border-color); background: var(--bg-card); color: var(--text-primary); font-size: 14px; }
        .filter-btn { padding: 8px 16px; border-radius: 5px; border: 2px solid var(--accent-color); background: var(--bg-card); color: var(--accent-color); cursor: pointer; font-weight: bold; transition: all 0.3s; }
        .filter-btn:hover { background: var(--accent-color); color: white; }
        .filter-btn.active { background: var(--accent-color); color: white; }
        .table-container { padding: 20px; overflow-x: auto; }
        table { width: 100%; border-collapse: collapse; box-shadow: 0 2px 10px var(--shadow); }
        thead { background: linear-gradient(135deg, var(--accent-color) 0%, var(--accent-hover) 100%); color: white; }
        th { padding: 15px 10px; text-align: left; font-weight: 600; position: sticky; top: 0; z-index: 10; font-size: 0.9em; }
        td { padding: 12px 10px; border-bottom: 1px solid var(--border-color); font-size: 0.85em; color: var(--text-primary); }
        tbody tr { background: var(--bg-card); transition: background-color 0.3s; }
        tbody tr:hover { background-color: var(--accent-light); }
        tbody tr.filtered-out { display: none; }
        .mol-img { max-width: 150px; max-height: 150px; border: 2px solid var(--border-color); border-radius: 8px; padding: 5px; background: white; }
        .smiles-text { font-family: 'Courier New', monospace; font-size: 0.8em; color: var(--text-secondary); max-width: 200px; word-wrap: break-word; }
        .affinity-cell { text-align: center; }
        .affinity-good { color: #22c55e; font-weight: bold; }
        .affinity-moderate { color: #eab308; font-weight: bold; }
        .affinity-weak { color: #ef4444; font-weight: bold; }
        .ki-value { font-size: 0.8em; color: var(--text-secondary); display: block; margin-top: 3px; }
        .descriptor-value { text-align: center; font-size: 0.9em; }
        .sa-good { color: #22c55e; font-weight: bold; }
        .sa-moderate { color: #eab308; font-weight: bold; }
        .sa-hard { color: #ef4444; font-weight: bold; }
        .footer { text-align: center; padding: 20px; background: var(--bg-secondary); color: var(--text-secondary); border-top: 2px solid var(--border-color); }
        .no-data { color: var(--text-secondary); font-style: italic; }
        .info-badge { display: inline-block; padding: 4px 8px; background: var(--accent-light); color: var(--accent-color); border-radius: 4px; font-size: 0.85em; font-weight: bold; margin: 0 5px; }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <button class="theme-toggle" onclick="toggleTheme()">🌓 Theme</button>
            <h1>{{ title }}</h1>
            <p>{{ subtitle }}</p>
            <p class="affinity-method-note">{{ affinity_note }}</p>
            {% if exhaustiveness %}<p><span class="info-badge">Exhaustiveness: {{ exhaustiveness }}</span></p>{% endif %}
        </div>
        <div class="stats">
            <div class="stat-box"><div class="number">{{ total_ligands }}</div><div class="label">Ligands</div></div>
            <div class="stat-box"><div class="number">{{ total_proteins }}</div><div class="label">Receptors</div></div>
            <div class="stat-box"><div class="number">{{ total_dockings }}</div><div class="label">Docking runs</div></div>
        </div>
        <div class="controls">
            <div class="control-group">
                <label for="sortProtein">Receptor:</label>
                <select id="sortProtein" onchange="sortTable()">
                    <option value="">-- Select --</option>
                    {% for protein in proteins %}<option value="p{{ loop.index0 }}">{{ protein }}</option>{% endfor %}
                </select>
            </div>
            <div class="control-group">
                <label for="sortCriteria">Sort by:</label>
                <select id="sortCriteria" onchange="sortTable()">
                    <option value="ki">Estimated Ki (ascending)</option>
                    <option value="pki">Estimated pKi (descending)</option>
                    <option value="sa">SA Score (↓ easier to synthesize)</option>
                </select>
            </div>
            <div class="control-group">
                <label>Filters:</label>
                <button class="filter-btn" onclick="toggleFilter('lipinski')" id="filterLipinski">Lipinski Ro5</button>
                <button class="filter-btn" onclick="toggleFilter('ki')" id="filterKi">Estimated Ki &lt; 1 µM</button>
                <button class="filter-btn" onclick="toggleFilter('sa')" id="filterSa">SA &lt; 4</button>
                <button class="filter-btn" onclick="clearFilters()">Reset</button>
            </div>
        </div>
        <div class="table-container">
            <table id="resultsTable">
                <thead>
                    <tr>
                        <th>№</th><th>Name</th><th>Formula</th><th>MW</th><th>LogP</th><th>TPSA</th><th>HBD</th><th>HBA</th><th>SA Score</th><th>SMILES</th><th>Structure</th>
                        {% for protein in proteins %}<th>{{ protein }}<br><span style="font-size:0.7em; font-weight:normal;">Docking score (kcal/mol) / estimated Ki / estimated pKi</span></th>{% endfor %}
                    </tr>
                </thead>
                <tbody id="tableBody">
                    {% for row in data %}
                    <tr data-ligand="{{ row.name }}" data-mw="{{ row.mw if row.mw is not none else 9999 }}" data-logp="{{ row.logp if row.logp is not none else 9999 }}" data-hbd="{{ row.hbd if row.hbd is not none else 99 }}" data-hba="{{ row.hba if row.hba is not none else 99 }}" {% for protein in proteins %}data-p{{ loop.index0 }}-ki="{{ row.proteins[protein].ki if row.proteins[protein].ki is not none else 999999 }}" data-p{{ loop.index0 }}-pki="{{ row.proteins[protein].pki if row.proteins[protein].pki is not none else 0 }}" {% endfor %}data-sa="{{ row.sa_score if row.sa_score is not none else 999 }}">
                        <td class="row-number">{{ loop.index }}</td>
                        <td><strong>{{ row.name }}</strong></td>
                        <td class="descriptor-value">{{ row.formula }}</td>
                        <td class="descriptor-value">{{ "%.2f"|format(row.mw) if row.mw is not none else "N/A" }}</td>
                        <td class="descriptor-value">{{ "%.2f"|format(row.logp) if row.logp is not none else "N/A" }}</td>
                        <td class="descriptor-value">{{ "%.2f"|format(row.tpsa) if row.tpsa is not none else "N/A" }}</td>
                        <td class="descriptor-value">{{ row.hbd if row.hbd is not none else "N/A" }}</td>
                        <td class="descriptor-value">{{ row.hba if row.hba is not none else "N/A" }}</td>
                        <td class="descriptor-value">{% if row.sa_score is not none %}{% if row.sa_score < 3.0 %}<span class="sa-good">{{ "%.2f"|format(row.sa_score) }}</span>{% elif row.sa_score < 6.0 %}<span class="sa-moderate">{{ "%.2f"|format(row.sa_score) }}</span>{% else %}<span class="sa-hard">{{ "%.2f"|format(row.sa_score) }}</span>{% endif %}{% else %}<span class="no-data">N/A</span>{% endif %}</td>
                        <td><div class="smiles-text">{{ row.smiles }}</div></td>
                        <td>{% if row.structure_img %}<img src="{{ row.structure_img }}" alt="{{ row.name }}" class="mol-img">{% else %}<span class="no-data">No image</span>{% endif %}</td>
                        {% for protein in proteins %}
                        <td class="affinity-cell">{% set pdata = row.proteins[protein] %}{% if pdata.delta_g is not none %}{% if pdata.delta_g < -8.0 %}<span class="affinity-good">{{ "%.2f"|format(pdata.delta_g) }}</span>{% elif pdata.delta_g < -6.0 %}<span class="affinity-moderate">{{ "%.2f"|format(pdata.delta_g) }}</span>{% else %}<span class="affinity-weak">{{ "%.2f"|format(pdata.delta_g) }}</span>{% endif %}<span class="ki-value">{{ pdata.ki_formatted }}</span>{% if pdata.pki is not none %}<span class="ki-value">pKi: {{ "%.2f"|format(pdata.pki) }}</span>{% endif %}{% else %}<span class="no-data">—</span>{% endif %}</td>
                        {% endfor %}
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
        </div>
        <div class="footer">
            <p>Report generated: {{ timestamp }}</p>
            <p>Project: {{ app_name }} {{ app_version }} | AutoDock Vina Docking</p>
            <p>SA Score: 1-3 (easy), 3-6 (moderate), &gt;6 (difficult) | Lipinski Ro5: MW≤500, LogP≤5, HBD≤5, HBA≤10</p>
        </div>
    </div>
    <script>
        function toggleTheme() { const html = document.documentElement; const currentTheme = html.getAttribute('data-theme'); const newTheme = currentTheme === 'dark' ? 'light' : 'dark'; html.setAttribute('data-theme', newTheme); localStorage.setItem('theme', newTheme); }
        (function() { const savedTheme = localStorage.getItem('theme') || 'light'; document.documentElement.setAttribute('data-theme', savedTheme); })();
        function sortTable() { const protein = document.getElementById('sortProtein').value; const criteria = document.getElementById('sortCriteria').value; if (!protein && criteria !== 'sa') { return; } const tbody = document.getElementById('tableBody'); const rows = Array.from(tbody.querySelectorAll('tr')).filter(r => !r.classList.contains('filtered-out')); rows.sort((a, b) => { let valA, valB; if (criteria === 'sa') { valA = parseFloat(a.getAttribute('data-sa')); valB = parseFloat(b.getAttribute('data-sa')); return valA - valB; } else if (criteria === 'ki') { valA = parseFloat(a.getAttribute(`data-${protein}-ki`)); valB = parseFloat(b.getAttribute(`data-${protein}-ki`)); return valA - valB; } else if (criteria === 'pki') { valA = parseFloat(a.getAttribute(`data-${protein}-pki`)); valB = parseFloat(b.getAttribute(`data-${protein}-pki`)); return valB - valA; } return 0; }); rows.forEach(row => tbody.appendChild(row)); updateRowNumbers(); }
        let activeFilters = new Set();
        function toggleFilter(filterType) { const btn = document.getElementById(`filter${filterType.charAt(0).toUpperCase() + filterType.slice(1)}`); if (activeFilters.has(filterType)) { activeFilters.delete(filterType); btn.classList.remove('active'); } else { activeFilters.add(filterType); btn.classList.add('active'); } applyFilters(); }
        function clearFilters() { activeFilters.clear(); document.querySelectorAll('.filter-btn.active').forEach(btn => btn.classList.remove('active')); applyFilters(); }
        function applyFilters() { const rows = document.querySelectorAll('#tableBody tr'); rows.forEach(row => { let show = true; if (activeFilters.has('lipinski')) { const mw = parseFloat(row.getAttribute('data-mw')); const logp = parseFloat(row.getAttribute('data-logp')); const hbd = parseInt(row.getAttribute('data-hbd')); const hba = parseInt(row.getAttribute('data-hba')); if (mw > 500 || logp > 5 || hbd > 5 || hba > 10) { show = false; } } if (activeFilters.has('ki') && show) { const kiAttrs = Array.from(row.attributes).filter(attr => attr.name.endsWith('-ki')); const hasGoodKi = kiAttrs.some(attr => parseFloat(attr.value) < 1e-6); if (!hasGoodKi) { show = false; } } if (activeFilters.has('sa') && show) { const sa = parseFloat(row.getAttribute('data-sa')); if (sa >= 4) { show = false; } } row.classList.toggle('filtered-out', !show); }); updateRowNumbers(); }
        function updateRowNumbers() { const visibleRows = Array.from(document.querySelectorAll('#tableBody tr')).filter(r => !r.classList.contains('filtered-out')); visibleRows.forEach((row, index) => { row.querySelector('.row-number').textContent = index + 1; }); }
    </script>
</body>
</html>"""
        
        # Rendering
        template = Template(html_template, autoescape=True)
        html_content = template.render(app_name=APP_NAME, app_version=__version__, affinity_note=AFFINITY_NOTE, title=title, subtitle=f"Generated {datetime.now().strftime('%Y-%m-%d at %H:%M')}", total_ligands=len(report_data), total_proteins=len(proteins), total_dockings=sum(len(values) for values in enriched_results.values()), proteins=proteins, data=report_data, timestamp=datetime.now().strftime('%d.%m.%Y %H:%M:%S'), exhaustiveness=exhaustiveness)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(brand_html(html_content))
        
        self.logger.info(f"✓ HTML report saved: {output_file}")
        return output_file
    
    def generate_csv_report(self,
                           enriched_results: Dict,
                           smiles_dict: Optional[Dict[str, str]] = None,
                           output_file: Optional[Path] = None) -> Path:
        """Generate a CSV report with results and descriptors."""
        self.logger.info("Generating CSV report...")
        
        if smiles_dict is None:
            smiles_dict = self.load_smiles_data()
        
        report_data, proteins = self.prepare_report_data(enriched_results, smiles_dict)
        
        if output_file is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = self.reports_dir / f"docking_report_{timestamp}.csv"
        
        # Build column headers
        headers = ['name', 'smiles', 'formula', 'mw', 'logp', 'tpsa', 'hbd', 'hba', 'sa_score']
        
        for protein in proteins:
            headers.extend([
                f'{protein}_delta_g',
                f'{protein}_ki',
                f'{protein}_ki_formatted',
                f'{protein}_pki'
            ])
        
        # Write the CSV file
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            
            for row in report_data:
                csv_row = [
                    row['name'],
                    row['smiles'],
                    row['formula'],
                    row['mw'] if row['mw'] is not None else '',
                    row['logp'] if row['logp'] is not None else '',
                    row['tpsa'] if row['tpsa'] is not None else '',
                    row['hbd'] if row['hbd'] is not None else '',
                    row['hba'] if row['hba'] is not None else '',
                    row['sa_score'] if row['sa_score'] is not None else ''
                ]
                
                for protein in proteins:
                    pdata = row['proteins'][protein]
                    csv_row.extend([
                        pdata['delta_g'] if pdata['delta_g'] is not None else '',
                        pdata['ki'] if pdata['ki'] is not None else '',
                        pdata['ki_formatted'],
                        pdata['pki'] if pdata['pki'] is not None else ''
                    ])
                
                writer.writerow(csv_row)
        
        Path(str(output_file) + ".notes.txt").write_text(AFFINITY_NOTE + "\nDocking score columns (delta_g) are in kcal/mol.\n", encoding="utf-8")
        self.logger.info(f"✓ CSV report saved: {output_file}")
        return output_file
    
    def generate_reports(self,
                        enriched_results: Dict,
                        smiles_dict: Optional[Dict[str, str]] = None,
                        title: str = "Docking results",
                        exhaustiveness: int = None) -> Tuple[Path, Path]:
        """Generate both reports (HTML and CSV)."""
        self.logger.info("REPORT GENERATION")
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        html_file = self.reports_dir / f"docking_report_{timestamp}.html"
        csv_file = self.reports_dir / f"docking_report_{timestamp}.csv"
        
        html_path = self.generate_html_report(
            enriched_results, smiles_dict, html_file, title, exhaustiveness
        )
        
        csv_path = self.generate_csv_report(
            enriched_results, smiles_dict, csv_file
        )
        
        self.logger.info("✓ Reports generated successfully")
        self.logger.info(f"  HTML: {html_path}")
        self.logger.info(f"  CSV: {csv_path}")
        
        return html_path, csv_path
    
    def generate_db_report(self, db_data: Dict) -> Path:
        """
        Generate an HTML report from the database (all experiment history).
        Use SMILES from the database to generate structures.
        
        Args:
            db_data: Data from database.export_for_report()
            
        Returns:
            Path to the generated db-report.html
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = self.reports_dir / f"db_report_{timestamp}.html"
        
        # Convert database data to the report format
        enriched_results = {}
        smiles_dict = {}
        
        for ordinal, result in enumerate(db_data['results'], 1):
            ligand_name = (f"Experiment {result.get('experiment_id', '?')} / "
                           f"record {result.get('id', ordinal)} / {result['ligand_name']}")
            protein_name = result['protein_name']
            ligand_smiles = result.get('ligand_smiles', '')
            
            # Store SMILES
            if ligand_smiles and ligand_name not in smiles_dict:
                smiles_dict[ligand_name] = ligand_smiles
            
            if ligand_name not in enriched_results:
                enriched_results[ligand_name] = {}
            
            enriched_results[ligand_name][protein_name] = {
                'affinity': result['delta_g'],
                'ki': result['ki'],
                'ki_formatted': self._format_ki(result['ki']) if result['ki'] is not None else 'N/A',
                'pki': result['pki']
            }
        
        # Generate HTML using SMILES from the database
        return self.generate_html_report(
            enriched_results,
            smiles_dict,
            output_file,
            title="Docking history (database)",
            exhaustiveness=None
        )

    def _format_ki(self, ki: float) -> str:
        """Format Ki as a human-readable string."""
        if ki < 1e-9:
            return f"{ki*1e12:.2f} pM"
        elif ki < 1e-6:
            return f"{ki*1e9:.2f} nM"
        elif ki < 1e-3:
            return f"{ki*1e6:.2f} µM"
        else:
            return f"{ki*1e3:.2f} mM"

