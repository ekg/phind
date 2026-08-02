#!/usr/bin/env python3
"""
pggb phage-parameter retry pilot — cluster_6 (single-cluster validation).

Runs pggb 0.6.0 (certified tool, `pggb-env`) on
`pggb_analysis/cluster_6/cluster_6.fa` with the corrected phage-oriented
parameter set from the user (chat 3):

    p=75  s=250  l=500  k=11  ani-diff=80

and compares the resulting graph / alignment / consensus outputs against the
historical run (`pggb_analysis/cluster_6/`, pipeline log: wfmash -p 85 -n 5
-k 19 -l 2000 -X, seqwish -k 19, smoothxg -X 100 -r 53, gfaffix).

pggb 0.6.0 CLI mapping (see artifacts/consumer_compatibility/tool_versions.json
certified help + pggb wrapper source):
    -p, --map-pct-id        -> wfmash identity threshold (%)            [p=75]
    -s, --segment-length    -> wfmash seed segment length (bp)          [s=250]
    -l, --block-length      -> wfmash min block length (bp)             [l=500]
    -k, --min-match-len     -> seqwish exact-match filter (bp)          [k=11]
    -g, --hg-filter-ani-diff-> wfmash hypergeometric ANI-diff filter    [ani-diff=80]
    -c, --n-mappings        -> wfmash mappings per segment (parity with
                               historical -n 5)                         [c=5]

Outputs (written under pggb_analysis/phage_params_retry/):
    *.final.gfa             final normalized graph GFA
    *.smooth.gfa            smoothxg-smoothed graph GFA (kept via -A)
    *.seqwish.gfa           seqwish induced graph GFA
    *.alignments.wfmash.paf wfmash base-level alignments
    *.log                   pggb run log
    paths.fa                all graph paths as FASTA (odgi paths -f)
    aln.fa                  MAFFT alignment of prophage paths
    consensus.fa            smoothxg Consensus_* path sequences
    ancestral.fa            majority-rule consensus + position confidence
    pilot_stats.json        machine-readable stats (new run + old run)
    comparison_table.md     old-vs-new comparison table

Usage:
    python3 scripts/pggb_phage_params_pilot.py [--skip-pggb] [--threads N]
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

BASE_DIR = Path('/home/erikg/phind/.wg-worktrees/agent-9')
WORK_DIR = BASE_DIR / 'pggb_analysis'
CLUSTER_DIR = WORK_DIR / 'cluster_6'
OUT_DIR = WORK_DIR / 'phage_params_retry'

# Tools — pggb 0.6.0 certified env; mafft from the historical pipeline env.
PGGB_ENV = '/home/erikg/micromamba/envs/pggb-env/bin'
PGGB_ENV_OLD = '/home/erikg/micromamba/envs/pggb_env/bin'
TOOLS = {
    'pggb': f'{PGGB_ENV}/pggb',
    'odgi': f'{PGGB_ENV}/odgi',
    'mafft': f'{PGGB_ENV_OLD}/mafft',
    'pigz': f'{PGGB_ENV}/pigz',
}

# Corrected (user, chat 3) parameter set.
PARAMS = {
    'p': 75,            # wfmash map-pct-id (% identity)
    's': 250,           # wfmash segment length (bp)
    'l': 500,           # wfmash block length (bp)
    'k': 11,            # seqwish min-match length (bp)
    'ani-diff': 80,     # wfmash hg-filter-ani-diff
    'c': 5,             # n-mappings (parity with historical run)
    'n': 53,            # n haplotypes
}

MIN_CORE_ALIGNMENT = 50


def run_cmd(cmd, desc=None, timeout=None, cwd=None, stdout_file=None, env=None):
    desc = desc or cmd[0]
    print(f"\n  [{desc}]")
    print(f"    {' '.join(str(c) for c in cmd)}", flush=True)
    t0 = time.time()
    try:
        if stdout_file:
            with open(stdout_file, 'w') as fout:
                res = subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE,
                                     text=True, timeout=timeout, cwd=cwd, env=env)
        else:
            res = subprocess.run(cmd, capture_output=True, text=True,
                                 timeout=timeout, cwd=cwd, env=env)
        dt = time.time() - t0
        print(f"    exit={res.returncode} time={dt:.1f}s", flush=True)
        if res.returncode != 0:
            err = (res.stderr or res.stdout or '')[-500:]
            print(f"    STDERR/OUT tail: {err}", flush=True)
            raise RuntimeError(f"Command failed (exit {res.returncode}): {cmd[0]}")
        return res
    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT after {timeout}s", flush=True)
        raise RuntimeError(f"Command timed out: {cmd[0]}")


def tool_versions():
    out = {}
    for name, path in TOOLS.items():
        if not os.path.isfile(path):
            out[name] = 'MISSING'
            continue
        try:
            r = subprocess.run([path, '--version'], capture_output=True, text=True, timeout=20)
            v = (r.stdout or r.stderr).strip().split('\n')[0]
            out[name] = v
        except Exception as e:
            out[name] = f'err:{e}'
    return out


def run_pggb(input_fasta, out_dir, threads, skip=False):
    """Run pggb 0.6.0 with the corrected params. Idempotent on final.gfa."""
    out_dir.mkdir(parents=True, exist_ok=True)
    final_gfas = list(out_dir.glob('*.final.gfa'))
    if skip or final_gfas:
        print(f"  pggb outputs present ({len(final_gfas)} final.gfa), skipping pggb.")
        return final_gfas[0] if final_gfas else None

    # pggb requires a .fai for -n auto-detection only; we pass -n explicitly.
    cmd = [
        TOOLS['pggb'],
        '-i', str(input_fasta),
        '-o', str(out_dir),
        '-p', str(PARAMS['p']),
        '-s', str(PARAMS['s']),
        '-l', str(PARAMS['l']),
        '-k', str(PARAMS['k']),
        '-g', str(PARAMS['ani-diff']),
        '-c', str(PARAMS['c']),
        '-n', str(PARAMS['n']),
        '-t', str(threads),
        '-A',          # keep intermediate graphs (seqwish gfa, PAFs)
        '--skip-viz',  # no PNG visualizations (not needed for validation)
    ]
    # The pggb wrapper invokes sibling tools (vg, multiqc, seqwish, ...) as
    # bare commands, so its env bin dir must be on PATH.
    env = dict(os.environ)
    env['PATH'] = PGGB_ENV + os.pathsep + env.get('PATH', '')
    run_cmd(cmd, desc='pggb 0.6.0', timeout=7200, env=env)
    final_gfas = sorted(out_dir.glob('*.final.gfa'))
    if not final_gfas:
        raise RuntimeError('pggb produced no *.final.gfa')
    return final_gfas[0]


def parse_gfa_stats(gfa_path):
    """Count S/L/P lines, total segment bp, N50 of segment lengths."""
    n_seg = n_edges = n_paths = 0
    total_bp = 0
    seg_lens = []
    with open(gfa_path) as f:
        for line in f:
            if line.startswith('S\t'):
                n_seg += 1
                parts = line.split('\t')
                ln = len(parts[1])
                seg_lens.append(ln)
                total_bp += ln
            elif line.startswith('L\t'):
                n_edges += 1
            elif line.startswith('P\t') or line.startswith('W\t'):
                n_paths += 1
    seg_lens.sort(reverse=True)
    n50 = n50_of(seg_lens)
    return {
        'segments': n_seg,
        'edges': n_edges,
        'paths': n_paths,
        'total_bp': total_bp,
        'segment_n50': n50,
        'segment_l50': l50_of(seg_lens, total_bp),
    }


def n50_of(lengths):
    total = sum(lengths)
    if not total:
        return 0
    half = total / 2
    acc = 0
    for ln in lengths:
        acc += ln
        if acc >= half:
            return ln
    return 0


def l50_of(lengths, total):
    if not total:
        return 0
    half = total / 2
    acc = 0
    n = 0
    for ln in lengths:
        acc += ln
        n += 1
        if acc >= half:
            return n
    return n


def extract_paths_fasta(gfa_path, out_fasta):
    """odgi paths -f -> FASTA of all graph paths."""
    run_cmd([TOOLS['odgi'], 'paths', '-i', str(gfa_path), '-f'],
            desc='odgi paths -f', timeout=1800, stdout_file=str(out_fasta))
    return out_fasta


def read_fasta(path):
    seqs = {}
    cur = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if cur is not None:
                    seqs[cur] = ''.join(buf)
                cur = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = ''.join(buf)
    return seqs


def path_stats(paths_fa):
    """Path statistics. N50/total/mean are computed over the prophage paths
    (non-Consensus_) so the metric is like-for-like with the historical
    cluster_6.paths.fa, which contains the 53 prophage paths only."""
    seqs = read_fasta(paths_fa)
    prophage = {k: v for k, v in seqs.items() if not k.startswith('Consensus_')}
    consensus = {k: v for k, v in seqs.items() if k.startswith('Consensus_')}
    plens = sorted((len(v) for v in prophage.values()), reverse=True)
    ptotal = sum(plens)
    clens = sorted((len(v) for v in consensus.values()), reverse=True)
    ctotal = sum(clens)
    return {
        'n_paths': len(seqs),
        'n_prophage_paths': len(prophage),
        'n_consensus_paths': len(consensus),
        'prophage_total_bp': ptotal,
        'prophage_n50': n50_of(plens),
        'prophage_l50': l50_of(plens, ptotal),
        'prophage_mean_len': round(ptotal / len(prophage), 1) if prophage else 0,
        'prophage_min_len': min(plens) if plens else 0,
        'prophage_max_len': max(plens) if plens else 0,
        'consensus_total_bp': ctotal,
        'consensus_n50': n50_of(clens),
    }


def alignment_rate(paf_path, n_queries_expected=None):
    """Per-query union coverage fraction from a wfmash PAF."""
    intervals = defaultdict(list)   # qname -> [(start,end)]
    qlens = {}
    total_span = 0
    n_records = 0
    with open(paf_path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 11:
                continue
            qname, qlen, qstart, qend = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
            qlens[qname] = qlen
            intervals[qname].append((qstart, qend))
            total_span += (qend - qstart)
            n_records += 1

    covs = []
    for qname, ivs in intervals.items():
        ivs.sort()
        merged = []
        for s, e in ivs:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        covered = sum(e - s for s, e in merged)
        qlen = qlens[qname]
        covs.append(covered / qlen if qlen else 0.0)

    return {
        'n_paf_records': n_records,
        'n_queries_with_hits': len(covs),
        'mean_query_coverage': round(sum(covs) / len(covs), 4) if covs else 0.0,
        'median_query_coverage': round(sorted(covs)[len(covs)//2], 4) if covs else 0.0,
        'min_query_coverage': round(min(covs), 4) if covs else 0.0,
        'max_query_coverage': round(max(covs), 4) if covs else 0.0,
        'total_aligned_span_bp': total_span,
    }


def mafft_align(paths_fa, out_aln, threads):
    run_cmd([TOOLS['mafft'], '--auto', '--thread', str(threads), str(paths_fa)],
            desc='mafft', timeout=3600, stdout_file=str(out_aln))
    return out_aln


def majority_consensus(seqs, out_fa, out_conf):
    """Majority-rule consensus over path sequences (mirrors pipeline
    reconstruct_ancestral_genome) with per-position confidence."""
    prophage = {k: v for k, v in seqs.items() if not k.startswith('Consensus_')}
    if not prophage:
        return None
    max_len = max(len(v) for v in prophage.values())
    consensus = []
    pos_counts = []
    for i in range(max_len):
        counts = Counter()
        for seq in prophage.values():
            if i < len(seq):
                b = seq[i].upper()
                if b in 'ACGT':
                    counts[b] += 1
        if counts:
            total = sum(counts.values())
            base, n = counts.most_common(1)[0]
            consensus.append(base)
            pos_counts.append(n / total)
        else:
            consensus.append('N')
            pos_counts.append(0.0)
    cseq = ''.join(consensus)
    conf = sum(pos_counts) / len(pos_counts) if pos_counts else 0.0
    with open(out_fa, 'w') as f:
        f.write(f'>cluster_6_ancestral_consensus n_paths={len(prophage)}\n{cseq}\n')
    with open(out_conf, 'w') as f:
        f.write('position\tbase\tfrequency\n')
        for i, (b, fr) in enumerate(zip(consensus, pos_counts)):
            f.write(f'{i}\t{b}\t{fr:.4f}\n')
    return {'length': len(cseq), 'confidence': round(conf, 4), 'n_seqs': len(prophage)}


def alignment_consensus(aln_fa, out_fa, out_conf):
    """Majority-rule consensus over an aligned FASTA (MAFFT). Column-wise
    majority ignoring gaps; per-column support = fraction of non-gap bases
    agreeing with the majority base."""
    seqs = read_fasta(aln_fa)
    seqs = {k: v for k, v in seqs.items() if not k.startswith('Consensus_')}
    if not seqs:
        return None
    names = list(seqs.keys())
    n = len(names)
    length = len(seqs[names[0]])
    consensus = []
    pos_counts = []
    for i in range(length):
        counts = Counter()
        for name in names:
            b = seqs[name][i].upper()
            if b in 'ACGT':
                counts[b] += 1
        if counts:
            base, cnt = counts.most_common(1)[0]
            consensus.append(base)
            total = sum(counts.values())
            pos_counts.append(cnt / total)
        else:
            consensus.append('N')
            pos_counts.append(0.0)
    cseq = ''.join(consensus)
    conf = sum(pos_counts) / len(pos_counts) if pos_counts else 0.0
    with open(out_fa, 'w') as f:
        f.write(f'>cluster_6_aln_consensus n_seqs={n}\n{cseq}\n')
    with open(out_conf, 'w') as f:
        f.write('position\tbase\tfrequency\n')
        for i, (b, fr) in enumerate(zip(consensus, pos_counts)):
            f.write(f'{i}\t{b}\t{fr:.4f}\n')
    return {'length': len(cseq), 'confidence': round(conf, 4), 'n_seqs': n,
            'columns': length}


def find_files(out_dir, patterns):
    found = {}
    for key, pat in patterns.items():
        hits = sorted(out_dir.glob(pat), key=lambda p: p.stat().st_mtime)
        found[key] = str(hits[-1]) if hits else None
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--skip-pggb', action='store_true',
                    help='skip the pggb run (use existing outputs)')
    ap.add_argument('--threads', type=int, default=64)
    args = ap.parse_args()

    input_fasta = CLUSTER_DIR / 'cluster_6.fa'
    if not input_fasta.exists():
        sys.exit(f'Missing input: {input_fasta}')

    print('=' * 70)
    print('pggb phage-parameter retry pilot (cluster_6)')
    print('=' * 70)
    print('Corrected params (user chat 3):', json.dumps(PARAMS))
    print('Tool versions:')
    for name, ver in tool_versions().items():
        print(f'  {name}: {ver}')

    # 1. pggb run
    final_gfa = run_pggb(input_fasta, OUT_DIR, args.threads, skip=args.skip_pggb)
    print(f'\nFinal GFA: {final_gfa}')

    files = find_files(OUT_DIR, {
        'final_gfa': '*.final.gfa',
        'smooth_gfa': '*.smooth.gfa',
        'seqwish_gfa': '*.seqwish.gfa',
        'alignments_paf': '*.alignments.wfmash.paf',
        'mappings_paf': '*.mappings.wfmash.paf',
        'pggb_log': '*.log',
    })
    print('Found outputs:', {k: (Path(v).name if v else None) for k, v in files.items()})
    if not files['alignments_paf']:
        print('WARNING: no alignments PAF kept; alignment rate will be missing.')
    if not files['smooth_gfa']:
        print('WARNING: no smooth GFA kept; consensus extraction will be limited.')

    # 2. Paths FASTA (new run)
    paths_fa = OUT_DIR / 'paths.fa'
    extract_paths_fasta(files['final_gfa'], paths_fa)
    new_path_stats = path_stats(paths_fa)

    # 3. MAFFT alignment of prophage paths (mirror historical aln.fa)
    seqs_all = read_fasta(paths_fa)
    prophage_seqs = {k: v for k, v in seqs_all.items() if not k.startswith('Consensus_')}
    lens_set = {len(v) for v in prophage_seqs.values()}
    aln_fa = OUT_DIR / 'aln.fa'
    new_aln = {}
    if aln_fa.exists() and aln_fa.stat().st_size > 0:
        aln_seqs = read_fasta(aln_fa)
        new_aln = {'n_seqs': len(aln_seqs),
                   'columns': len(next(iter(aln_seqs.values()))) if aln_seqs else 0,
                   'method': 'mafft'}
        print(f'  aln.fa exists ({aln_fa.stat().st_size} B), reusing MAFFT alignment.')
    elif len(lens_set) <= 1:
        aln_len = next(iter(lens_set)) if lens_set else 0
        with open(aln_fa, 'w') as f:
            for name, seq in prophage_seqs.items():
                f.write(f'>{name}\n{seq}\n')
        new_aln = {'n_seqs': len(prophage_seqs), 'columns': aln_len, 'method': 'graph-paths'}
    else:
        mafft_align(paths_fa, aln_fa, max(1, args.threads // 4))
        aln_seqs = read_fasta(aln_fa)
        new_aln = {'n_seqs': len(aln_seqs),
                   'columns': len(next(iter(aln_seqs.values()))) if aln_seqs else 0,
                   'method': 'mafft'}

    # 4. Consensus / ancestral sequence
    #    pggb 0.6.0 default run does not emit Consensus_* paths (consensus-spec
    #    is off), so two pipeline-style consensuses are produced:
    #    consensus.fa        majority-rule over the MAFFT-aligned paths
    #    ancestral.fa        majority-rule over the raw graph paths (mirrors
    #                        scripts/pggb_per_cluster_pipeline.py
    #                        reconstruct_ancestral_genome)
    consensus = {}
    anc = None
    if aln_fa.exists() and aln_fa.stat().st_size > 0:
        anc = alignment_consensus(aln_fa, OUT_DIR / 'consensus.fa',
                                  OUT_DIR / 'consensus_confidence.tsv')
    if anc is not None:
        consensus = {'method': 'majority-rule (MAFFT-aligned paths)'}
    anc_raw = majority_consensus(seqs_all, OUT_DIR / 'ancestral.fa',
                                 OUT_DIR / 'ancestral_confidence.tsv')

    # 5. Stats for NEW run
    new_gfa_stats = parse_gfa_stats(files['final_gfa'])
    new_aln_rate = (alignment_rate(files['alignments_paf'])
                    if files['alignments_paf'] else None)
    new_seqwish_stats = (parse_gfa_stats(files['seqwish_gfa'])
                         if files['seqwish_gfa'] else None)

    # 6. Stats for OLD run (historical cluster_6)
    old_final_gfa = CLUSTER_DIR / 'cluster_6.fixed.gfa'
    old_align_paf = CLUSTER_DIR / 'cluster_6.alignments.paf'
    old_paths_fa = CLUSTER_DIR / 'cluster_6.paths.fa'
    old_gfa_stats = parse_gfa_stats(old_final_gfa)
    old_aln_rate = alignment_rate(old_align_paf)
    old_path_stats = path_stats(old_paths_fa)
    old_seqwish_stats = parse_gfa_stats(CLUSTER_DIR / 'cluster_6.seqwish.gfa')
    old_consensus_fa = CLUSTER_DIR / 'ancestral' / '6_ancestral.fa'
    old_anc = None
    if old_consensus_fa.exists():
        s = read_fasta(old_consensus_fa)
        if s:
            v = next(iter(s.values()))
            old_anc = {'length': len(v), 'confidence': None, 'n_seqs': None}
    old_aln_fa = CLUSTER_DIR / 'cluster_6.aln.fa'
    old_aln = None
    if old_aln_fa.exists():
        s = read_fasta(old_aln_fa)
        if s:
            old_aln = {'n_seqs': len(s),
                       'columns': len(next(iter(s.values()))),
                       'method': 'mafft-or-paths'}

    # 7. Assemble stats JSON
    stats = {
        'params_corrected': PARAMS,
        'params_historical': {
            'note': 'scripts/pggb_per_cluster_pipeline.py log (pggb_analysis/pipeline.log)',
            'wfmash': '-p 85 -n 5 -k 19 -l 2000 -X',
            'seqwish': '-k 19',
            'smoothxg': '-X 100 -r 53',
            'gfaffix': True,
        },
        'tool_versions': tool_versions(),
        'new_run': {
            'dir': str(OUT_DIR),
            'final_gfa': files['final_gfa'],
            'graph': new_gfa_stats,
            'seqwish': new_seqwish_stats,
            'paths': new_path_stats,
            'alignment_rate': new_aln_rate,
            'alignment': new_aln,
            'consensus': consensus,
            'ancestral': anc,          # majority-rule over MAFFT alignment
            'ancestral_raw_paths': anc_raw,  # majority-rule over raw graph paths
        },
        'old_run': {
            'dir': str(CLUSTER_DIR),
            'graph': old_gfa_stats,
            'seqwish': old_seqwish_stats,
            'paths': old_path_stats,
            'alignment_rate': old_aln_rate,
            'alignment': old_aln,
            'ancestral': old_anc,
        },
    }
    with open(OUT_DIR / 'pilot_stats.json', 'w') as f:
        json.dump(stats, f, indent=2)

    # 8. Comparison table (markdown)
    def g(d, k):
        return d.get(k) if isinstance(d, dict) else None

    rows = [
        ('Graph segments (final)', g(stats['new_run']['graph'], 'segments'),
         g(stats['old_run']['graph'], 'segments'), 'lower = less fragmentation'),
        ('Graph edges (final)', g(stats['new_run']['graph'], 'edges'),
         g(stats['old_run']['graph'], 'edges'), None),
        ('Graph total bp (final)', g(stats['new_run']['graph'], 'total_bp'),
         g(stats['old_run']['graph'], 'total_bp'), None),
        ('Segment N50 (final, bp)', g(stats['new_run']['graph'], 'segment_n50'),
         g(stats['old_run']['graph'], 'segment_n50'), 'higher = longer clean nodes'),
        ('Seqwish segments (pre-smooth)', g(stats['new_run'].get('seqwish') or {}, 'segments'),
         g(stats['old_run']['seqwish'], 'segments'), None),
        ('Seqwish paths', g(stats['new_run'].get('seqwish') or {}, 'paths'),
         g(stats['old_run']['seqwish'], 'paths'), None),
        ('Final graph paths (P/W lines)', g(stats['new_run']['graph'], 'paths'),
         g(stats['old_run']['graph'], 'paths'), 'old had chimeric 69kb traversal path'),
        ('Prophage paths', g(stats['new_run']['paths'], 'n_prophage_paths'),
         g(stats['old_run']['paths'], 'n_prophage_paths'), 'should = 53'),
        ('Prophage path total bp', g(stats['new_run']['paths'], 'prophage_total_bp'),
         g(stats['old_run']['paths'], 'prophage_total_bp'), None),
        ('Prophage path N50 (bp)', g(stats['new_run']['paths'], 'prophage_n50'),
         g(stats['old_run']['paths'], 'prophage_n50'), None),
        ('Prophage path mean len (bp)', g(stats['new_run']['paths'], 'prophage_mean_len'),
         g(stats['old_run']['paths'], 'prophage_mean_len'), None),
        ('Consensus paths in graph', g(stats['new_run']['paths'], 'n_consensus_paths'),
         g(stats['old_run']['paths'], 'n_consensus_paths'), 'smoothxg POA consensus paths'),
        ('Mean query coverage (PAF)', g(stats['new_run'].get('alignment_rate') or {}, 'mean_query_coverage'),
         g(stats['old_run']['alignment_rate'], 'mean_query_coverage'), 'higher = more prophage captured'),
        ('PAF records', g(stats['new_run'].get('alignment_rate') or {}, 'n_paf_records'),
         g(stats['old_run']['alignment_rate'], 'n_paf_records'), None),
        ('Alignment columns', g(stats['new_run'].get('alignment') or {}, 'columns'),
         g(stats['old_run'].get('alignment') or {}, 'columns'), None),
        ('Ancestral consensus len (raw paths, bp)', g(stats['new_run'].get('ancestral_raw_paths') or {}, 'length'),
         g(stats['old_run'].get('ancestral') or {}, 'length'), 'like-for-like with old pipeline ancestral'),
        ('Ancestral confidence (raw paths)', g(stats['new_run'].get('ancestral_raw_paths') or {}, 'confidence'),
         g(stats['old_run'].get('ancestral') or {}, 'confidence'), None),
        ('Consensus len (MAFFT-majority, bp)', g(stats['new_run'].get('ancestral') or {}, 'length'),
         g(stats['old_run'].get('ancestral') or {}, 'length'), 'new-only (alignment-based)'),
        ('Consensus method', g(stats['new_run'].get('consensus') or {}, 'method'),
         'pipeline majority-rule (raw paths)', None),
    ]
    lines = [
        '# cluster_6 pggb parameter comparison (old vs corrected phage params)',
        '',
        f"Corrected params (this run): `p=75 s=250 l=500 k=11 ani-diff=80` via pggb 0.6.0 (`-p 75 -s 250 -l 500 -k 11 -g 80 -c 5 -n 53`).",
        '',
        f"Historical params: `wfmash -p 85 -n 5 -k 19 -l 2000 -X`, `seqwish -k 19`, `smoothxg -X 100 -r 53`, `gfaffix`.",
        '',
        '| metric | new (p=75,s=250,l=500,k=11,ani80) | old (p=85,l=2000,k=19) | notes |',
        '|---|---|---|---|',
    ]
    for label, nv, ov, note in rows:
        nvs = f'{nv:,}' if isinstance(nv, int) else str(nv)
        ovs = f'{ov:,}' if isinstance(ov, int) else str(ov)
        lines.append(f'| {label} | {nvs} | {ovs} | {note or ""} |')
    table = '\n'.join(lines)
    print('\n' + table)
    with open(OUT_DIR / 'comparison_table.md', 'w') as f:
        f.write(table + '\n')

    # 9. Copy report skeleton path for the writeup
    print(f'\nStats JSON: {OUT_DIR / "pilot_stats.json"}')
    print(f'Comparison table: {OUT_DIR / "comparison_table.md"}')
    print('DONE')


if __name__ == '__main__':
    main()
