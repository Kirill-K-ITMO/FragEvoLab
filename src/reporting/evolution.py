from src.reporting.branding import brand_html
"""Offline HTML/SVG evidence report. No CDN, browser network or executable input."""
from html import escape
from collections import Counter
import math
from statistics import mean
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, Draw
from src.version import APP_NAME, __version__

def esc(value):
    return escape(str(value), quote=True)

def fmt(value):
    return f'{value:.3f}' if isinstance(value, (float,int)) and not isinstance(value,bool) else esc(value if value is not None else 'no data')

def table(headers, rows):
    return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+esc(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+fmt(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'

def line_chart(title, series, x_label='Generation'):
    """Shared axes with explicit units; None breaks the polyline."""
    points=[(x,y) for _, values in series for x,y in values if y is not None and math.isfinite(y)]
    if not points: return '<section><h2>'+esc(title)+'</h2><p>No data</p></section>'
    xs,ys=zip(*points); xmin,xmax=min(xs),max(xs); ymin,ymax=min(ys),max(ys)
    if xmin==xmax: xmax=xmin+1
    if ymin==ymax: ymin-=.5; ymax+=.5
    def px(x): return 75+600*(x-xmin)/(xmax-xmin)
    def py(y): return 235-185*(y-ymin)/(ymax-ymin)
    svg='<svg viewBox="0 0 730 290" role="img" aria-label="'+esc(title)+'">'
    for i in range(5):
        y=ymin+(ymax-ymin)*i/4
        svg+=f'<line x1="75" x2="675" y1="{py(y):.2f}" y2="{py(y):.2f}" stroke="#dfe6ed"/><text x="65" y="{py(y)+4:.2f}" text-anchor="end">{y:.2f}</text>'
    for x in sorted(set(xs))[::max(1,len(set(xs))//10)]:
        svg+=f'<text x="{px(x):.2f}" y="255" text-anchor="middle">{x}</text>'
    colors=['#087e8b','#d75f26','#7055ad','#257a3e','#b52c65']
    legend=''
    for i,(name,values) in enumerate(series):
        color=colors[i%len(colors)]; segment=[]
        for x,y in [*values,(None,None)]:
            if y is None or not math.isfinite(y):
                if segment: svg+=f'<polyline points="{" ".join(segment)}" fill="none" stroke="{color}" stroke-width="2"/>'
                segment=[]; continue
            segment.append(f'{px(x):.2f},{py(y):.2f}')
            svg+=f'<circle cx="{px(x):.2f}" cy="{py(y):.2f}" r="3" fill="{color}"><title>{esc(name)}: G{x}, {y:.4f}</title></circle>'
        legend+=f'<span style="color:{color}">{esc(name)}</span> '
    return '<section><h2>'+esc(title)+'</h2><p>'+legend+'</p>'+svg+f'<text x="375" y="282" text-anchor="middle">{esc(x_label)}</text></svg></section>'

def bars(title, values):
    if not values: return ''
    maximum=max(abs(v) for _,v in values) or 1
    return '<section><h2>'+esc(title)+'</h2>'+''.join('<div class="barrow"><span>'+esc(k)+'</span><meter min="0" max="'+str(maximum)+'" value="'+str(max(0,v))+'"></meter><b>'+fmt(v)+'</b></div>' for k,v in values)+'</section>'

def structure(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None: return '<p>Structure unavailable</p>'
    drawer=Draw.rdMolDraw2D.MolDraw2DSVG(330,200)
    drawer.DrawMolecule(mol); drawer.FinishDrawing()
    svg=drawer.GetDrawingText()
    return svg[svg.index('<svg'):]

def diversity(smiles):
    """Deterministic bounded subset, documented as a sample above 200 structures."""
    smiles=sorted(set(smiles))[:200]
    fp=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    fps=[fp.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]
    distances=[1-DataStructs.TanimotoSimilarity(a,b) for i,a in enumerate(fps) for b in fps[:i]]
    return mean(distances) if distances else None

def write_report(report, path):
    generations=report.get('generations',[])
    proteins=[p['name'] if isinstance(p,dict) else p for p in report['config']['targets']+report['config'].get('off_targets',[])]
    best=[]; avg=[]; div=[]; complete=[]; feasible=[]; caches=[]; dockings=[]
    for g in generations:
        n=g['generation']; rows=list(g['results'].values())
        scores=[r['score'].get('fitness') for r in rows if r['score'].get('complete') and r['score'].get('fitness') is not None]
        best.append((n,max(scores) if scores else None)); avg.append((n,mean(scores) if scores else None))
        div.append((n,diversity(g['molecules'])))
        complete.append((n,sum(bool(r['score'].get('complete')) for r in rows)))
        feasible.append((n,sum(bool(r['score'].get('feasible')) for r in rows)))
        caches.append((n,g['cache_hits']));dockings.append((n,g['new_dockings']))
    sections=[line_chart('Fitness: higher is better',[('Best fitness',best),('Mean fitness',avg)]),
        line_chart('Diversity: mean 1 - Tanimoto similarity',[('Morgan r=2, 2048 bits',div)]),
        '<p>Diversity: up to 200 unique structures per generation, ordered by SMILES; for larger populations this is a sample estimate. Undefined for a single structure.</p>',
        line_chart('Complete and feasible candidates',[('Complete docking results',complete),('Target thresholds satisfied',feasible)]),
        line_chart('Calculations and cache',[('New docking calculations',dockings),('Cache hits',caches)])]
    energy=[]
    for p in proteins:
        series=[]
        for g in generations:
            vals=[r['energies'][p] for r in g['results'].values() if p in r['energies']]
            series.append((g['generation'],mean(vals) if vals else None))
        energy.append((p,series))
    sections.append(line_chart('Mean docking score, kcal/mol',energy))
    sections.append(bars('Generation attempt outcomes',sorted(Counter(a['status'] for a in report.get('attempts',[])).items())))
    reasons=Counter(e.get('reason',e.get('status')) for a in report.get('attempts',[]) for e in a.get('events',[]) if e.get('reason') or e.get('status')!='product')
    sections.append(bars('Operator rejection reasons',sorted((k,v) for k,v in reasons.items() if k)))
    records=report.get('fragments',{}).get('records',{})
    ordered=sorted(records.items(),key=lambda x:(-x[1]['evaluated'],-x[1]['weight'],x[0]))
    sections.append('<section><h2>Fragment weight history</h2><p>The top 10 fragments by evaluated replacement count are shown. A weight is not a probability of activity.</p>'+line_chart('Adaptive weights',[(f'F{i}',[(g['generation'],g.get('fragment_weights',{}).get(key)) for g in generations]) for i,(key,r) in enumerate(ordered[:10],1)])+'</section>')
    sections.append('<section><h2>Fragment statistics</h2><p>The improvement rate includes evaluated single replacements only. Association between a fragment and a favorable docking score does not establish the contribution of that fragment.</p>'+table(['ID','SMILES','Attempts','Valid','Duplicates','Evaluated','Improved','Improvement rate','Mean reward','Weight'],[(f'F{i}',key,r['attempted'],r['valid'],r['duplicate'],r['evaluated'],r['improved'],r['improved']/r['evaluated'] if r['evaluated'] else None,r['reward_sum']/r['evaluated'] if r['evaluated'] else None,r['weight']) for i,(key,r) in enumerate(ordered,1)])+'</section>')
    sections.append('<section><h2>Fragment structures: first 20</h2><div class="cards">'+''.join('<article><h3>F'+str(i)+'</h3>'+structure(key)+'<code>'+esc(key)+'</code></article>' for i,(key,r) in enumerate(ordered[:20],1))+'</div></section>')
    sections.append('<section><h2>Operators</h2>'+table(['Operator','Attempts','Valid','Duplicates','Evaluated','Improved','Evaluation failures'],[(k,*[v.get(field,0) for field in ['attempted','valid','duplicate','evaluated','improved','evaluation_failed']]) for k,v in sorted(report.get('operators',{}).items())])+'</section>')
    sections.append('<section><h2>Top molecules</h2><div class="cards">'+''.join('<article><h3>G'+str(c['generation'])+' / '+esc(c['origin'])+'</h3>'+structure(c['smiles'])+'<code>'+esc(c['smiles'])+'</code><p>Fitness: '+fmt(c['score'].get('fitness'))+'</p>'+table(['Receptor','kcal/mol'],[(p,c['energies'].get(p)) for p in proteins])+'<details><summary>Parent and mutations</summary><p>'+esc(c.get('parent'))+'</p><pre>'+esc(c.get('events',[]))+'</pre></details></article>' for c in report.get('top_molecules',[]))+'</div></section>')
    sections.append('<section><h2>All candidates</h2><input id="filter" placeholder="Search SMILES" aria-label="Search SMILES"><div id="candidates">'+table(['SMILES','G','Origin','Fitness','Feasible',*proteins],[(c['smiles'],c['generation'],c['origin'],c['score'].get('fitness'),c['score'].get('feasible'),*[c['energies'].get(p) for p in proteins]) for c in report.get('candidates',report.get('top_molecules',[]))])+'</div></section>')
    style='body{font:15px system-ui;background:#f0f4f8;color:#193349;margin:0 auto;max-width:1250px;padding:24px}h1{font-size:32px}section,article{background:white;border-radius:12px;padding:20px;margin:18px 0}svg{max-width:100%;height:auto}svg text{font:12px system-ui;fill:#193349}.scroll{overflow:auto}table{border-collapse:collapse;width:100%}td,th{border-bottom:1px solid #ddd;padding:10px;text-align:left;max-width:400px;overflow-wrap:anywhere}th{background:#eaf1f7}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:16px}article{border:1px solid #ccd9e1;margin:0}code,pre{overflow-wrap:anywhere;white-space:pre-wrap}input{padding:12px;width:80%}.barrow{display:flex;gap:12px;margin:12px}.barrow span{width:240px}.barrow meter{flex:1}p{line-height:1.6}'
    document='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>FragEvoLab - Evolution report</title><style>'+style+'</style></head><body><h1>'+APP_NAME+' '+__version__+' - Evolution</h1><p>Run: '+esc(report['run_id'])+' | Stop reason: '+esc(report['stop_reason'])+'</p><p>Unique molecules: '+str(report['unique_molecules'])+' | New docking calculations: '+str(report['total_dockings'])+' | Cache: '+str(report['cache_hits'])+'</p><p>Docking scores are computational estimates, not measured activity. The underlying data are available in report.json, candidates.csv, and checkpoint.json.</p>'+''.join(sections)+'<script>document.getElementById("filter").addEventListener("input",function(){const q=this.value.toLowerCase();document.querySelectorAll("#candidates tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));});</script></body></html>'
    path.write_text(brand_html(document),encoding='utf-8')
