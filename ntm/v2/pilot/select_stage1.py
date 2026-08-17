#!/usr/bin/env python3
"""Preregistered Stage-1 pilot selection (execute-bounded-ntm).

Policy fixed BEFORE any Stage-1 results are seen (2026-08-17, post Stage-0
submission, pre Stage-1 submission). Deterministic; reruns byte-identical.

Selection (24 submissions, GenBank_RefSeq, threshold 0.7):
  A. one representative interior_module bait per A_strong_candidate clade
     (first by manifest order); for clade 0_0021 (no interior_module bait)
     its member_interior bait instead  -> 15 baits
  B. top-5 junction baits restricted to purely reconstructed joins
     (observed_adjacency_count == 0), ranked by co_occurrence_count DESC,
     then clade_id ASC, then bait_id ASC                        ->  5 baits
  C. two positive controls: mid-genome 1200 bp slices of complete public
     reference genomes D29 (AF022214.2) and L5 (NC_001335.1)    ->  2 baits
  D. two negative controls: first two control_shuffled_negative baits
     by bait_id ASC                                             ->  2 baits
"""
import csv, hashlib, json, os, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
PANEL = os.path.join(REPO, 'ntm', 'v2', 'bait', 'panel', 'bait_manifest.tsv')
BAITS_FA = os.path.join(os.path.dirname(PANEL), 'baits.fa')
EXT_DIR = '/mnt/nvme3n1/erikg/phind-genome-work/ntm/v2/external_validation/panel'
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT_DIR, exist_ok=True)

def read_fasta(path):
    recs, name, chunks = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if name: recs[name] = ''.join(chunks)
            name, chunks = line[1:].split()[0], []
        elif line:
            chunks.append(line)
    if name: recs[name] = ''.join(chunks)
    return recs

rows = list(csv.DictReader(open(PANEL), delimiter='\t'))
baits = read_fasta(BAITS_FA)

# A: one interior_module per A-tier clade (first in manifest order)
sel_a, seen = [], set()
for r in rows:
    if r['tier'] != 'A_strong_candidate': continue
    if r['clade_id'] in seen: continue
    if r['bait_class'] == 'interior_module':
        sel_a.append(r); seen.add(r['clade_id'])
# clade 0_0021 fallback: member_interior
for r in rows:
    if r['tier'] == 'A_strong_candidate' and r['clade_id'] not in seen and r['bait_class'] == 'member_interior':
        sel_a.append(r); seen.add(r['clade_id'])
        break

# B: top-5 purely reconstructed junctions by (co_occurrence DESC, clade ASC, bait_id ASC)
juncs = [r for r in rows if r['bait_class'] == 'junction' and r['observed_adjacency_count'] == '0']
juncs.sort(key=lambda r: (-int(r['co_occurrence_count']), r['clade_id'], r['bait_id']))
sel_b = juncs[:5]

# C: positive controls — D29, L5 mid-genome slices
ext_recs = read_fasta(os.path.join(EXT_DIR, 'panel.mycobacteriophage.fasta'))
man = {r['panel_id']: r for r in csv.DictReader(open(os.path.join(EXT_DIR, 'manifest_panel.tsv')), delimiter='\t')}
pos = []
for panel_id in ('PHPUB-000587', 'PHPUB-001403'):
    seq = ext_recs[panel_id]
    m = man[panel_id]
    L = len(seq)
    start, end = L // 2 - 600, L // 2 + 600
    pos.append({'bait_id': f'POSCON_{m["phage_name"].upper()}_{m["primary_accession"].replace(".", "_")}_MID1200',
                'seq': seq[start:end], 'source_panel_id': panel_id,
                'source_accession': m['primary_accession'], 'phage': m['phage_name'],
                'source_len_bp': L, 'slice_start0': start, 'slice_end0': end,
                'seq_sha256_source_record': m['seq_sha256']})

# D: negative controls — first two shuffled by bait_id ASC
negs = sorted([r for r in rows if r['bait_class'] == 'control_shuffled_negative'], key=lambda r: r['bait_id'])[:2]

# ---- emit stage1 panel FASTA + manifest ----
out_rows = []
with open(os.path.join(OUT_DIR, 'stage1_panel.fa'), 'w') as fa:
    def emit(bait_id, seq):
        fa.write(f'>{bait_id}\n')
        for i in range(0, len(seq), 70):
            fa.write(seq[i:i+70] + '\n')
    for r in sel_a + sel_b:
        bid = r['bait_id']
        seq = baits[bid]
        emit(bid, seq)
        out_rows.append({'bait_id': bid, 'role': 'ntm_bait', 'bait_class': r['bait_class'],
                         'clade_id': r['clade_id'], 'tier': r['tier'], 'module': r.get('module',''),
                         'length_bp': len(seq), 'seq_sha256': hashlib.sha256(seq.encode()).hexdigest(),
                         'source': 'ntm/v2/bait/panel/baits.fa'})
    for p in pos:
        emit(p['bait_id'], p['seq'])
        out_rows.append({'bait_id': p['bait_id'], 'role': 'positive_control', 'bait_class': 'control_public_reference',
                         'clade_id': '', 'tier': '', 'module': '',
                         'length_bp': len(p['seq']), 'seq_sha256': hashlib.sha256(p['seq'].encode()).hexdigest(),
                         'source': f"ext_panel:{p['source_panel_id']}:{p['source_accession']}[{p['slice_start0']}:{p['slice_end0']}] len={p['source_len_bp']}"})
    for r in negs:
        bid = r['bait_id']
        seq = baits[bid]
        emit(bid, seq)
        out_rows.append({'bait_id': bid, 'role': 'negative_control', 'bait_class': r['bait_class'],
                         'clade_id': '', 'tier': '', 'module': '',
                         'length_bp': len(seq), 'seq_sha256': hashlib.sha256(seq.encode()).hexdigest(),
                         'source': 'ntm/v2/bait/panel/baits.fa'})

fields = ['bait_id','role','bait_class','clade_id','tier','module','length_bp','seq_sha256','source']
with open(os.path.join(OUT_DIR, 'stage1_selection.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
    w.writeheader(); w.writerows(out_rows)

meta = {'policy': __doc__.strip(), 'n_submissions': len(out_rows),
        'group': 'GenBank_RefSeq', 'threshold': 0.7,
        'junction_b_support': [(r['bait_id'], r['clade_id'], r['co_occurrence_count']) for r in sel_b],
        'generated_at_utc': '2026-08-17T21:00:00Z'}
json.dump(meta, open(os.path.join(OUT_DIR, 'stage1_selection_meta.json'), 'w'), indent=2)
print(f"selected {len(out_rows)} submissions: A={len(sel_a)} B={len(sel_b)} C={len(pos)} D={len(negs)}")
for r in out_rows: print(' ', r['bait_id'], r['role'], r['length_bp'])
