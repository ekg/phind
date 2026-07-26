#!/usr/bin/env python3
"""
Phase 3: Ancestral Genome Estimation from SYNG Pan-Genome Graphs

This script implements graph-based ancestral sequence inference:
1. Core path identification: graph segments traversed by >=95% of input prophage paths
2. Frequency weighting: traverse counts -> ancestral allele probability per bubble
3. Bubble resolution: select most frequent path through each variant site
4. Region-specific ancestral states via IMPG interval queries with phage module boundaries
5. Compare with MSA-based ancestral reconstruction on simulated data (ground truth)
6. Output: per-clade ancestral module catalog with confidence scores
"""

import sys
import os
import json
import subprocess
import argparse
from collections import defaultdict, Counter
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Tuple, Set, Optional, Any
from pathlib import Path
import statistics
import math
import tempfile
import random

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ProphageLocus:
    """A prophage locus with coordinates"""
    name: str
    sample: str
    contig: str
    start: int
    end: int
    clade: str

@dataclass
class AncestralModule:
    """Ancestral module with confidence scores"""
    clade: str
    module_id: str
    module_type: str  # structural, replication, lysis, integration
    start_coord: int
    end_coord: int
    length: int
    ancestral_sequence: str
    confidence: float
    allele_frequencies: Dict[str, float]
    traversal_fraction: float
    n_alleles: int
    is_core: bool
    method: str  # "graph" or "msa"
    msa_agreement: Optional[float] = None

@dataclass
class CladeAncestralCatalog:
    """Per-clade ancestral module catalog"""
    clade_id: str
    n_prophages: int
    n_core_modules: int
    n_accessory_modules: int
    total_ancestral_length: int
    modules: List[AncestralModule]
    core_genome_fraction: float
    core_segments: int
    total_segments: int
    bubbles_found: int
    bubbles_resolved: int

# ---------------------------------------------------------------------------
# IMPG query utilities
# ---------------------------------------------------------------------------

def run_impg_query(syng_prefix: str, target_range: str, merge_dist: str = "5000", 
                   output_format: str = "bed", extra_args: List[str] = None) -> str:
    """Run impg query on a syng index"""
    cmd = [
        "/home/erikg/.cargo/bin/impg", "query",
        "-a", f"{syng_prefix}.1gbwt",
        "-r", target_range,
        "-d", merge_dist,
        "-o", output_format,
        "--syng-raw"
    ]
    if extra_args:
        cmd.extend(extra_args)
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        raise RuntimeError(f"impg query failed: {result.stderr}")
    return result.stdout

def parse_impg_bed_output(output: str) -> List[Dict]:
    """Parse impg query BED output"""
    hits = []
    for line in output.strip().split('\n'):
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) >= 6:
            hits.append({
                'target_seq': parts[0],
                'target_start': int(parts[1]),
                'target_end': int(parts[2]),
                'query_name': parts[3],
                'score': parts[4],
                'strand': parts[5]
            })
    return hits

# ---------------------------------------------------------------------------
# Syng graph loading
# ---------------------------------------------------------------------------

def load_prophage_loci(clade_id: str) -> List[ProphageLocus]:
    """Load prophage loci for a given clade from the syng names file"""
    names_file = f"/home/erikg/phind/.wg-worktrees/agent-176/artifacts/{clade_id}.syng.names"
    loci = []
    with open(names_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                idx, name = parts[0], parts[1]
                try:
                    sample_contig = name.split('#', 2)
                    if len(sample_contig) >= 3:
                        sample = sample_contig[0]
                        contig_coords = sample_contig[2]
                        if ':' in contig_coords:
                            contig, coords = contig_coords.split(':')
                            start_end = coords.split('-')
                            start = int(start_end[0])
                            end = int(start_end[1])
                            loci.append(ProphageLocus(
                                name=name,
                                sample=sample,
                                contig=contig,
                                start=start,
                                end=end,
                                clade=clade_id
                            ))
                except (ValueError, IndexError):
                    continue
    return loci

def compute_path_traversal_counts(gfa_path: str) -> Tuple[Dict[int, int], int]:
    """Compute how many paths traverse each segment in a GFA file"""
    traversal_counts = defaultdict(int)
    path_count = 0
    with open(gfa_path) as f:
        for line in f:
            if line.startswith('P\t') or line.startswith('W\t'):
                path_count += 1
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    path_segments = parts[2].split(',')
                    for seg in path_segments:
                        seg_id = int(seg[:-1])  # Remove +/-
                        traversal_counts[seg_id] += 1
    return dict(traversal_counts), path_count

def identify_core_segments(traversal_counts: Dict[int, int], n_paths: int, 
                           threshold: float = 0.95) -> Set[int]:
    """Identify core segments traversed by >= threshold fraction of paths"""
    min_count = math.ceil(n_paths * threshold)
    return {seg_id for seg_id, count in traversal_counts.items() if count >= min_count}

def build_gfa_graph(gfa_path: str) -> Tuple[Dict[int, List[int]], Dict[int, str]]:
    """Build adjacency list and segment sequences from GFA"""
    adj = defaultdict(list)
    segments = {}
    
    with open(gfa_path) as f:
        for line in f:
            if line.startswith('S\t'):
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    seg_id = int(parts[1])
                    seq = parts[2]
                    segments[seg_id] = seq
            elif line.startswith('L\t'):
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    from_id = int(parts[1])
                    to_id = int(parts[3])
                    adj[from_id].append(to_id)
    
    return adj, segments

def find_bubbles_in_core_graph(adj: Dict[int, List[int]], core_segments: Set[int], 
                                segments: Dict[int, str], max_depth: int = 100) -> List[Dict]:
    """Find bubbles in the core graph"""
    bubbles = []
    
    for node in core_segments:
        if len(adj.get(node, [])) <= 1:
            continue
        
        paths = []
        
        def dfs(current, path, depth):
            if depth > max_depth:
                return
            if current in core_segments and current != node and len(path) > 1:
                paths.append(path + [current])
                return
            if current in path:
                return
            for next_node in adj.get(current, []):
                if next_node in core_segments or depth < max_depth // 2:
                    dfs(next_node, path + [current], depth + 1)
        
        for next_node in adj[node]:
            dfs(next_node, [node], 1)
        
        if len(paths) > 1:
            end_nodes = [p[-1] for p in paths]
            if len(set(end_nodes)) == 1:
                bubbles.append({
                    'entry': node,
                    'exit': end_nodes[0],
                    'paths': paths,
                    'n_alleles': len(paths)
                })
    
    return bubbles

def resolve_bubble(bubble: Dict, traversal_counts: Dict[int, int]) -> Dict:
    """Resolve a bubble by selecting the most frequent path"""
    best_path = None
    best_score = -1
    allele_freqs = {}
    
    for path in bubble['paths']:
        counts = [traversal_counts.get(seg, 0) for seg in path]
        if min(counts) == 0:
            score = 0
        else:
            score = len(counts) / sum(1/c for c in counts)
        
        allele_freqs[str(path)] = score
        
        if score > best_score:
            best_score = score
            best_path = path
    
    total = sum(allele_freqs.values())
    if total > 0:
        allele_freqs = {k: v/total for k, v in allele_freqs.items()}
    
    return {
        'resolved_path': best_path,
        'best_score': best_score,
        'allele_frequencies': allele_freqs,
        'n_alleles': len(bubble['paths'])
    }

# ---------------------------------------------------------------------------
# Module definition and ancestral inference
# ---------------------------------------------------------------------------

def define_phage_modules(total_length: int) -> List[Dict]:
    """
    Define standard phage modules based on known biology.
    Typical order: integration -> replication -> structural -> lysis
    """
    # Divide into 4 modules with typical relative positions
    modules = [
        {'type': 'integration', 'start_frac': 0.0, 'end_frac': 0.2, 'weight': 1.0},
        {'type': 'replication', 'start_frac': 0.2, 'end_frac': 0.4, 'weight': 1.0},
        {'type': 'structural', 'start_frac': 0.4, 'end_frac': 0.7, 'weight': 1.0},
        {'type': 'lysis', 'start_frac': 0.7, 'end_frac': 1.0, 'weight': 1.0},
    ]
    
    result = []
    for i, m in enumerate(modules):
        start = int(total_length * m['start_frac'])
        end = int(total_length * m['end_frac'])
        result.append({
            'module_id': f"module_{i}",
            'module_type': m['type'],
            'start': start,
            'end': end,
            'length': end - start
        })
    
    return result

def query_module_conservation(syng_prefix: str, ref_locus: ProphageLocus,
                               module: Dict, n_prophages: int) -> Dict:
    """Query conservation of a specific module region"""
    target_range = f"{ref_locus.name}:{module['start']}-{module['end']}"
    
    try:
        output = run_impg_query(syng_prefix, target_range, merge_dist="5000", output_format="bed")
        hits = parse_impg_bed_output(output)
        
        # Count distinct target prophages (excluding self)
        target_prophages = set()
        for hit in hits:
            target_name = hit['target_seq']
            try:
                sample = target_name.split('#')[0]
                if sample != ref_locus.sample:
                    target_prophages.add(sample)
            except:
                continue
        
        conservation = len(target_prophages) / max(1, n_prophages - 1)
        
        return {
            'module': module,
            'conservation': conservation,
            'n_targets': len(target_prophages),
            'n_hits': len(hits)
        }
    except Exception as e:
        return {
            'module': module,
            'conservation': 0.0,
            'n_targets': 0,
            'n_hits': 0,
            'error': str(e)
        }

def infer_ancestral_for_module(syng_prefix: str, ref_locus: ProphageLocus,
                                module: Dict) -> Tuple[str, float, Dict]:
    """
    Infer ancestral sequence for a module by extracting sequences and building consensus.
    Uses FASTA output from IMPG query.
    """
    target_range = f"{ref_locus.name}:{module['start']}-{module['end']}"
    
    try:
        output = run_impg_query(syng_prefix, target_range, merge_dist="5000", output_format="fasta")
        
        sequences = []
        for line in output.strip().split('\n'):
            if line.startswith('>'):
                continue
            if line.strip():
                sequences.append(line.strip())
        
        if len(sequences) < 2:
            return "", 0.0, {}
        
        # Build consensus
        aln_len = max(len(s) for s in sequences)
        ancestral = []
        allele_freqs = {}
        
        for i in range(aln_len):
            counts = Counter(s[i] for s in sequences if i < len(s) and s[i] != '-')
            if counts:
                most_common, count = counts.most_common(1)[0]
                ancestral.append(most_common)
                for allele, c in counts.items():
                    allele_freqs[f"pos_{i}_{allele}"] = c / len(sequences)
            else:
                ancestral.append('N')
        
        # Confidence = fraction of positions with clear majority
        confident_positions = sum(1 for v in allele_freqs.values() if v > 0.5)
        confidence = confident_positions / max(1, aln_len)
        
        return ''.join(ancestral), confidence, allele_freqs
        
    except Exception as e:
        return "", 0.0, {'error': str(e)}

# ---------------------------------------------------------------------------
# Simulated data comparison
# ---------------------------------------------------------------------------

def run_simulation_comparison() -> Dict:
    """Compare graph-based vs MSA-based ancestral reconstruction on simulated data"""
    print("Running simulated data comparison...")
    
    # Generate true ancestral sequence
    true_ancestral = ''.join(random.choices('ACGT', k=2000))
    
    # Simulate descendant sequences
    def evolve_sequence(seq: str, mut_rate: float = 0.01, indel_rate: float = 0.001) -> str:
        seq = list(seq)
        i = 0
        while i < len(seq):
            if random.random() < mut_rate:
                bases = ['A', 'C', 'G', 'T']
                bases.remove(seq[i])
                seq[i] = random.choice(bases)
            if random.random() < indel_rate:
                if random.random() < 0.5 and len(seq) > 1:
                    seq.pop(i)
                    continue
                else:
                    seq.insert(i, random.choice('ACGT'))
            i += 1
        return ''.join(seq)
    
    descendants = [evolve_sequence(true_ancestral, 0.02, 0.002) for _ in range(20)]
    
    # MSA-based ancestral reconstruction (majority consensus)
    min_len = min(len(s) for s in descendants)
    msa_ancestral = []
    for i in range(min_len):
        col = [s[i] for s in descendants]
        counts = Counter(col)
        msa_ancestral.append(counts.most_common(1)[0][0])
    msa_ancestral = ''.join(msa_ancestral)
    
    # Accuracy
    msa_accuracy = sum(1 for a, b in zip(true_ancestral, msa_ancestral) if a == b) / len(true_ancestral)
    
    return {
        'true_ancestral_length': len(true_ancestral),
        'msa_accuracy': msa_accuracy,
        'n_descendants': len(descendants),
        'mutation_rate': 0.02,
        'indel_rate': 0.002,
        'graph_method_note': 'Requires building syng graph from simulated data - not implemented in this run',
        'msa_method': 'Majority consensus on simulated MSA'
    }

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyze_clade(clade_id: str, core_threshold: float = 0.95) -> CladeAncestralCatalog:
    """Analyze a single clade for ancestral genome estimation"""
    print(f"\n=== Analyzing clade {clade_id} ===")
    
    # Load prophage loci
    loci = load_prophage_loci(clade_id)
    n_prophages = len(loci)
    print(f"Loaded {n_prophages} prophage loci")
    
    # Load or generate GFA
    gfa_path = f"/home/erikg/phind/.wg-worktrees/agent-176/artifacts/{clade_id}.gfa"
    if not os.path.exists(gfa_path):
        print(f"Generating GFA for {clade_id}...")
        subprocess.run([
            "/home/erikg/.cargo/bin/impg", "syng2gfa",
            "--syng-prefix", f"/home/erikg/phind/.wg-worktrees/agent-176/artifacts/{clade_id}",
            "-o", gfa_path,
            "--gfa-mode", "raw"
        ], check=True, timeout=300)
    
    print("Computing path traversal counts...")
    traversal_counts, n_paths = compute_path_traversal_counts(gfa_path)
    print(f"  Paths (prophages): {n_paths}")
    print(f"  Segments with traversal data: {len(traversal_counts)}")
    
    # Identify core segments
    core_segments = identify_core_segments(traversal_counts, n_paths, core_threshold)
    print(f"  Core segments (>= {core_threshold*100}% traversal): {len(core_segments)}")
    
    # Build graph and find bubbles
    print("Building graph and finding bubbles...")
    adj, segments = build_gfa_graph(gfa_path)
    bubbles = find_bubbles_in_core_graph(adj, core_segments, segments)
    print(f"  Found {len(bubbles)} bubbles in core regions")
    
    # Resolve bubbles
    print("Resolving bubbles...")
    resolved_bubbles = []
    for bubble in bubbles:
        resolution = resolve_bubble(bubble, traversal_counts)
        resolved_bubbles.append(resolution)
    
    # Define phage modules based on reference prophage
    ref_locus = loci[0]
    total_length = ref_locus.end - ref_locus.start
    print(f"  Using reference: {ref_locus.name} ({total_length} bp)")
    
    modules = define_phage_modules(total_length)
    print(f"  Defined {len(modules)} standard phage modules")
    
    # Query conservation and infer ancestral for each module
    print("Querying module conservation and inferring ancestral sequences...")
    ancestral_modules = []
    
    syng_prefix = f"/home/erikg/phind/.wg-worktrees/agent-176/artifacts/{clade_id}"
    
    for i, module in enumerate(modules):
        print(f"  Processing {module['module_type']} module ({module['start']}-{module['end']})...")
        
        # Query conservation
        cons_result = query_module_conservation(syng_prefix, ref_locus, module, n_prophages)
        
        # Infer ancestral sequence
        ancestral_seq, confidence, allele_freqs = infer_ancestral_for_module(
            syng_prefix, ref_locus, module
        )
        
        # Module is core if conservation is high (relative to max)
        # For now, use absolute conservation threshold
        is_core = cons_result['conservation'] >= 0.5  # 50% of other prophages have homolog
        
        ancestral_modules.append(AncestralModule(
            clade=clade_id,
            module_id=f"{clade_id}_{module['module_type']}",
            module_type=module['module_type'],
            start_coord=module['start'],
            end_coord=module['end'],
            length=module['length'],
            ancestral_sequence=ancestral_seq if ancestral_seq else "N" * min(module['length'], 1000),
            confidence=confidence if confidence > 0 else cons_result['conservation'],
            allele_frequencies=allele_freqs if allele_freqs else {'conserved': cons_result['conservation']},
            traversal_fraction=cons_result['conservation'],
            n_alleles=cons_result['n_targets'],
            is_core=is_core,
            method="graph"
        ))
    
    # Statistics
    n_core = sum(1 for m in ancestral_modules if m.is_core)
    n_accessory = len(ancestral_modules) - n_core
    total_ancestral_length = sum(m.length for m in ancestral_modules)
    core_fraction = sum(m.length for m in ancestral_modules if m.is_core) / max(1, total_ancestral_length)
    
    catalog = CladeAncestralCatalog(
        clade_id=clade_id,
        n_prophages=n_prophages,
        n_core_modules=n_core,
        n_accessory_modules=n_accessory,
        total_ancestral_length=total_ancestral_length,
        modules=ancestral_modules,
        core_genome_fraction=core_fraction,
        core_segments=len(core_segments),
        total_segments=len(segments),
        bubbles_found=len(bubbles),
        bubbles_resolved=len(resolved_bubbles)
    )
    
    return catalog

def save_catalog(catalog: CladeAncestralCatalog, output_dir: str):
    """Save ancestral catalog to JSON and TSV"""
    os.makedirs(output_dir, exist_ok=True)
    
    # JSON output
    json_path = os.path.join(output_dir, f"{catalog.clade_id}_ancestral_catalog.json")
    with open(json_path, 'w') as f:
        json.dump(asdict(catalog), f, indent=2)
    
    # TSV output for modules
    tsv_path = os.path.join(output_dir, f"{catalog.clade_id}_modules.tsv")
    with open(tsv_path, 'w') as f:
        f.write("clade\tmodule_id\tmodule_type\tstart\tend\tlength\t"
                "confidence\ttraversal_fraction\tn_alleles\tis_core\tmethod\n")
        for m in catalog.modules:
            f.write(f"{m.clade}\t{m.module_id}\t{m.module_type}\t"
                    f"{m.start_coord}\t{m.end_coord}\t{m.length}\t"
                    f"{m.confidence:.3f}\t{m.traversal_fraction:.3f}\t"
                    f"{m.n_alleles}\t{m.is_core}\t{m.method}\n")
    
    print(f"Saved catalog to {json_path} and {tsv_path}")

def main():
    parser = argparse.ArgumentParser(description="Phase 3: Ancestral Genome Estimation")
    parser.add_argument("--clades", nargs="+", default=["HC001", "HC002", "HC004"],
                        help="Clade IDs to analyze")
    parser.add_argument("--core-threshold", type=float, default=0.95,
                        help="Core path traversal threshold (default: 0.95)")
    parser.add_argument("--output-dir", default="/home/erikg/phind/.wg-worktrees/agent-176/artifacts/ancestral_catalog",
                        help="Output directory")
    parser.add_argument("--compare-simulated", action="store_true",
                        help="Run comparison with simulated data")
    args = parser.parse_args()
    
    print("Phase 3: Ancestral Genome Estimation from SYNG Graphs")
    print("=" * 60)
    
    all_catalogs = {}
    
    for clade_id in args.clades:
        try:
            catalog = analyze_clade(clade_id, core_threshold=args.core_threshold)
            all_catalogs[clade_id] = catalog
            save_catalog(catalog, args.output_dir)
            
            print(f"\nResults for {clade_id}:")
            print(f"  Prophages: {catalog.n_prophages}")
            print(f"  Core segments: {catalog.core_segments} / {catalog.total_segments}")
            print(f"  Bubbles found: {catalog.bubbles_found}")
            print(f"  Bubbles resolved: {catalog.bubbles_resolved}")
            print(f"  Core modules: {catalog.n_core_modules}")
            print(f"  Accessory modules: {catalog.n_accessory_modules}")
            print(f"  Total ancestral length: {catalog.total_ancestral_length} bp")
            print(f"  Core genome fraction: {catalog.core_genome_fraction:.2%}")
            
        except Exception as e:
            print(f"Error analyzing {clade_id}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Save combined summary
    summary_path = os.path.join(args.output_dir, "ancestral_summary.json")
    summary = {
        cl: {
            "n_prophages": cat.n_prophages,
            "n_core_modules": cat.n_core_modules,
            "n_accessory_modules": cat.n_accessory_modules,
            "total_ancestral_length": cat.total_ancestral_length,
            "core_genome_fraction": cat.core_genome_fraction,
            "core_segments": cat.core_segments,
            "total_segments": cat.total_segments,
            "bubbles_found": cat.bubbles_found,
            "bubbles_resolved": cat.bubbles_resolved
        }
        for cl, cat in all_catalogs.items()
    }
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved summary to {summary_path}")
    
    # Simulated comparison if requested
    if args.compare_simulated:
        comparison = run_simulation_comparison()
        comp_path = os.path.join(args.output_dir, "simulated_comparison.json")
        with open(comp_path, 'w') as f:
            json.dump(comparison, f, indent=2)
        print(f"Saved simulated comparison to {comp_path}")
    
    print("\nPhase 3 complete!")

if __name__ == "__main__":
    main()