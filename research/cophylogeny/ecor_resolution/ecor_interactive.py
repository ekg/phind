#!/usr/bin/env python3
"""Interactive SVG tanglegram for the ECOR-resolution prophage vs host comparison."""
import re, os, csv, html, json

RES = "research"
OUT = os.path.join(RES, "cophylogeny", "ecor_resolution")
os.makedirs(OUT, exist_ok=True)

def parse_newick(s):
    s = re.sub(r"\s+", "", s.strip())
    def rec(i):
        if s[i] == "(":
            i += 1; children = []
            while True:
                ch, i = rec(i); children.append(ch)
                c = s[i]
                if c == ",": i += 1
                else:
                    assert c == ")"; i += 1
                    name, i = read_name(i)
                    return {"children": children, "name": name}, i
        else:
            name, i = read_name(i)
            return {"name": name, "children": []}, i
    def read_name(i):
        j = i
        while j < len(s) and s[j] not in "(),:":
            j += 1
        name = s[i:j]
        if j < len(s) and s[j] == ":":
            j += 1
            while j < len(s) and s[j] not in "(),":
                j += 1
        return name, j
    return rec(0)[0]

def count_leaves(n):
    n["_nl"] = sum(count_leaves(c) for c in n["children"]) if n["children"] else 1
    return n["_nl"]

def ladderize(n):
    for c in n["children"]: ladderize(c)
    if n["children"]: n["children"].sort(key=lambda c: c["_nl"])
    return n

def nodes(n):
    yield n
    for c in n["children"]: yield from nodes(c)

def assign_y(n, ys):
    if not n["children"]:
        i = len(ys); ys.append(n["name"]); n["_y"] = float(i); return
    for c in n["children"]: assign_y(c, ys)
    n["_y"] = (n["children"][0]["_y"] + n["children"][-1]["_y"]) / 2.0

def assign_depth(n, d=0.0):
    n["_d"] = d
    for c in n["children"]: assign_depth(c, d + 1.0)

def prune(n, keep):
    if not n["children"]:
        return n if n["name"] in keep else None
    kids = [prune(c, keep) for c in n["children"]]
    kids = [k for k in kids if k is not None]
    if not kids: return None
    if len(kids) == 1: return kids[0]
    n["children"] = kids
    return n

# ---------------- read + prune ----------------
host_tree = parse_newick(open(os.path.join(RES, "host_mash_tree", "host_tree.nwk")).read())
phage_tree = parse_newick(open(os.path.join(RES, "mash_tree", "full_prophages_tree.nwk")).read())

ecor = list(csv.DictReader(open(os.path.join(RES, "ecor", "ecor_phylogroups_known.tsv")), delimiter="\t"))
ecor_hosts = [r["gcf_accession"] for r in ecor]
pg_of = {r["gcf_accession"]: r["known_phylogroup"] for r in ecor}
mani = list(csv.DictReader(open(os.path.join(RES, "ecor", "ecor_manifest.csv"))))
pros = [r["prophage_id"] for r in mani]

hpr = prune(host_tree, set(ecor_hosts))
ppr = prune(phage_tree, set(pros))
count_leaves(hpr); count_leaves(ppr)
ladderize(hpr); ladderize(ppr)
hys = []; assign_y(hpr, hys)
pys = []; assign_y(ppr, pys)
assign_depth(hpr); assign_depth(ppr)
max_depth_h = max(n["_d"] for n in nodes(hpr))
max_depth_p = max(n["_d"] for n in nodes(ppr))
Nhost = len(hys); Npro = len(pys)

# ---------------- layout ----------------
tipgap = 15.0
W = 2000
H = 60 + max(Nhost, Npro) * tipgap
box = {"x0h": 80, "x1h": 600, "x0p": 1380, "x1p": 1920}
ybase = 36

def host_y(name): return ybase + hys.index(name) * tipgap
def pro_y(name): return ybase + pys.index(name) * tipgap
def hx(n): return box["x0h"] + (n["_d"] / max_depth_h) * (box["x1h"] - box["x0h"])
def px(n): return box["x1p"] - (n["_d"] / max_depth_p) * (box["x1p"] - box["x0p"])

pg_color = {"A": "#e41a1c", "B1": "#377eb8", "B2": "#4daf4a", "C": "#984ea3",
            "D": "#ff7f00", "E": "#ffffb0", "F": "#a65628", "G": "#f781bf"}
def hostcolor(h): return pg_color.get(pg_of.get(h, ""), "#888")
def procolor(p): return hostcolor(re.sub("_prophage_.*", "", p))

# ---------------- SVG tree emitter ----------------
svg_parts = []
def emit_tree(n, orient):
    if not n["children"]:
        if orient == "h":
            y = host_y(n["name"]); x1 = box["x1h"]
            n["_ymin"] = n["_ymax"] = y
            col = hostcolor(n["name"])
            svg_parts.append(f'<line x1="{hx(n)}" y1="{y}" x2="{x1}" y2="{y}" stroke="#444" stroke-width="1.6"/>')
            svg_parts.append(f'<circle class="hhl" data-h="{html.escape(n["name"])}" cx="{x1-4}" cy="{y}" r="6" fill="{col}" stroke="#fff" pointer-events="all"/>')
            svg_parts.append(f'<text x="{x1+7}" y="{y+4}" font-size="11" fill="{col}">{html.escape(n["name"])}</text>')
        else:
            y = pro_y(n["name"]); x0 = box["x0p"]
            n["_ymin"] = n["_ymax"] = y
            col = procolor(n["name"])
            svg_parts.append(f'<line x1="{px(n)}" y1="{y}" x2="{x0}" y2="{y}" stroke="#444" stroke-width="1.6"/>')
            svg_parts.append(f'<circle cx="{x0+4}" cy="{y}" r="2.6" fill="{col}"/>')
        return y
    ymin = emit_tree(n["children"][0], orient)
    ymax = emit_tree(n["children"][-1], orient)
    yself = (ymin + ymax) / 2
    n["_ymin"] = min(c["_ymin"] for c in n["children"])
    n["_ymax"] = max(c["_ymax"] for c in n["children"])
    if orient == "h":
        x = hx(n)
        svg_parts.append(f'<line x1="{x}" y1="{n["_ymin"]}" x2="{x}" y2="{n["_ymax"]}" stroke="#444" stroke-width="1.6"/>')
        for c in n["children"]:
            cy = (c["_ymin"] + c["_ymax"]) / 2
            svg_parts.append(f'<line x1="{x}" y1="{cy}" x2="{hx(c)}" y2="{cy}" stroke="#444" stroke-width="1.6"/>')
    else:
        x = px(n)
        svg_parts.append(f'<line x1="{x}" y1="{n["_ymin"]}" x2="{x}" y2="{n["_ymax"]}" stroke="#444" stroke-width="1.6"/>')
        for c in n["children"]:
            cy = (c["_ymin"] + c["_ymax"]) / 2
            svg_parts.append(f'<line x1="{x}" y1="{cy}" x2="{px(c)}" y2="{cy}" stroke="#444" stroke-width="1.6"/>')
    return yself

emit_tree(hpr, "h")
emit_tree(ppr, "p")

# ---------------- links ----------------
link_from = box["x1h"] + 30
link_to = box["x0p"] - 50
link_svg = []
for p in pros:
    h = re.sub("_prophage_.*", "", p)
    yh = host_y(h); yp = pro_y(p)
    col = procolor(p)
    d = f'M{link_from},{yh} C {(link_from+link_to)/2},{yh} {(link_from+link_to)/2},{yp} {link_to},{yp}'
    link_svg.append(f'<path class="link" data-h="{html.escape(h)}" d="{d}" fill="none" stroke="{col}" stroke-opacity="0.3" stroke-width="1.1"/>')

SVG = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Arial,Helvetica,sans-serif">
<style>
 .link {{ transition: stroke-opacity .12s, stroke-width .12s; }}
 .hhl {{ cursor: pointer; }}
</style>
{chr(10).join(svg_parts)}
{chr(10).join(link_svg)}
<text x="{(box["x0h"]+box["x1h"])/2}" y="16" text-anchor="middle" font-size="14" font-weight="bold" fill="#222">Host phylogeny (72 ECOR strains, colored by phylogroup)</text>
<text x="{(box["x0p"]+box["x1p"])/2}" y="16" text-anchor="middle" font-size="14" font-weight="bold" fill="#222">Prophage phylogeny (300 ECOR prophages)</text>
<text id="tip" x="24" y="{H-12}" font-size="13" fill="#111"></text>
</svg>
'''

json_phylo = json.dumps(pg_of)
html_doc = f'''<!DOCTYPE html><html><head><meta charset="utf-8">
<title>ECOR tanglegram: prophage vs host cophylogeny</title>
<style>
 body {{ background:#fafafa; font-family:Arial; margin:0; padding:14px; }}
 h1 {{ font-size:19px; margin:4px 0 6px; }}
 .legend span {{ display:inline-block; width:13px; height:13px; margin-right:4px; vertical-align:middle; }}
 .legend {{ margin-bottom:8px; }}
</style></head><body>
<h1>Prophage &ndash; host cophylogeny &mdash; ECOR resolution (72 strains / 300 prophages)</h1>
<div class="legend">
{''.join(f'<span style="background:{pg_color[g]}"></span>{g}&nbsp;&nbsp;' for g in pg_color)}
</div>
{SVG}
<script>
var P = {json_phylo};
var links = document.querySelectorAll('.link');
var tip = document.getElementById('tip');
function clear(){{ links.forEach(function(l){{ l.setAttribute('stroke-opacity','0.3'); l.setAttribute('stroke-width','1.1'); }}); tip.textContent=''; }}
function highlight(h){{
  clear(); var any=false;
  links.forEach(function(l){{ if(l.getAttribute('data-h')===h){{ l.setAttribute('stroke-opacity','0.95'); l.setAttribute('stroke-width','3'); any=true; }} }});
  if(any) tip.textContent = 'Host strain '+h+' (phylogroup '+(P[h]||'?')+') associates with its resident prophage lineage(s)';
}}
document.querySelectorAll('.hhl').forEach(function(g){{
  g.addEventListener('mouseenter', function(){{ highlight(g.getAttribute('data-h')); }});
  g.addEventListener('mouseleave', clear);
}});
</script>
</body></html>
'''

with open(os.path.join(OUT, "ecor_tanglegram_interactive.html"), "w") as f:
    f.write(html_doc)
print("Wrote interactive HTML, hosts:", Nhost, "prophages:", Npro)
