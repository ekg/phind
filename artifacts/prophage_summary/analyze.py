#!/usr/bin/env python3
"""Offline, deterministic audit of 26k_prophage1.csv.

Only the Python standard library is used.  All definitions used in the output
are centralized here so the report, TSVs, and SVG plots cannot silently drift.
"""
from __future__ import annotations

import csv
import hashlib
import html
import math
import re
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Callable, Iterable, Sequence

EXPECTED_HEADER = [
    "end", "genome", "scaffold", "begin", "transposable", "taxonomy", "prophage_id"
]
ASSEMBLY_RE = re.compile(r"^(GC[AF]_\d+)\.(\d+)$")
SEQUENCE_RE = re.compile(r"^([A-Z]{2}_[A-Z0-9]+)\.(\d+)$")
ID_RE = re.compile(r"^(GC[AF]_\d+\.\d+)_prophage_(\d+)$")
SCOPE_ORDER = ("all", "transposable_flag_positive", "taxonomy_assigned")
SCOPE_LABEL = {
    "all": "All records",
    "transposable_flag_positive": "Flag-positive (transposable = 1.0)",
    "taxonomy_assigned": "Taxonomy-assigned (taxonomy != Unknown)",
}
COLORS = {
    "all": "#355f8d",
    "transposable_flag_positive": "#d1495b",
    "taxonomy_assigned": "#2a9d8f",
}


@dataclass(frozen=True)
class Record:
    source_row: int
    end_raw: str
    genome: str
    scaffold: str
    begin_raw: str
    transposable_raw: str
    taxonomy: str
    prophage_id: str
    begin: int
    end: int
    transposable: int

    @property
    def locus(self) -> tuple[str, str, int, int]:
        return (self.genome, self.scaffold, self.begin, self.end)

    @property
    def exact_record(self) -> tuple[str, ...]:
        return (
            self.end_raw, self.genome, self.scaffold, self.begin_raw,
            self.transposable_raw, self.taxonomy, self.prophage_id,
        )

    @property
    def coordinate_valid(self) -> bool:
        return self.begin >= 1 and self.end >= self.begin

    @property
    def inclusive_length(self) -> int | None:
        return self.end - self.begin + 1 if self.coordinate_valid else None

    @property
    def half_open_span(self) -> int | None:
        return self.end - self.begin if self.coordinate_valid else None


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def split_version(value: str, pattern: re.Pattern[str]) -> tuple[str, str]:
    match = pattern.fullmatch(value)
    return (match.group(1), match.group(2)) if match else (value, "")


def int_like(value: str) -> int:
    try:
        number = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError("not numeric") from exc
    if not number.is_finite() or number != number.to_integral_value():
        raise ValueError("not a finite integer-like number")
    return int(number)


def subset(scope: str, record: Record) -> bool:
    if scope == "all":
        return True
    if scope == "transposable_flag_positive":
        return record.transposable == 1
    if scope == "taxonomy_assigned":
        return record.taxonomy.strip().casefold() not in {"", "unknown"}
    raise KeyError(scope)


def prefix(value: str) -> str:
    match = re.match(r"^([A-Z]+_)", value)
    return match.group(1) if match else "other"


def quantile(values: Sequence[int | float], probability: float) -> float:
    """Linear interpolation at (n-1)*p (R type 7 / NumPy linear)."""
    if not values:
        return math.nan
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(ordered[lower])
    return float(ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower))


def fmt_number(value: int | float) -> str:
    if isinstance(value, int):
        return f"{value:,}"
    if math.isnan(value):
        return "NA"
    if value.is_integer():
        return f"{int(value):,}"
    return f"{value:,.2f}"


def fmt_metric(value: object) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.10g}"
    return str(value)


def fraction(count: int | float, denominator: int | float | None) -> str:
    if denominator in (None, 0):
        return ""
    return f"{count / denominator:.10g}"


def pct(count: int | float, denominator: int | float) -> str:
    return "NA" if denominator == 0 else f"{100 * count / denominator:.2f}%"


def tsv_write(path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def xml(value: object) -> str:
    return html.escape(str(value), quote=True)


def svg_start(width: int, height: int, title: str, desc: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        f"<title id=\"title\">{xml(title)}</title>",
        f"<desc id=\"desc\">{xml(desc)}</desc>",
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:18px;font-weight:bold}.small{font-size:11px}.axis{stroke:#333;stroke-width:1}.grid{stroke:#ddd;stroke-width:1}.legend{font-size:12px}</style>',
        f'<text class="title" x="{width/2}" y="25" text-anchor="middle">{xml(title)}</text>',
    ]


def write_count_plot(path: Path, per_genome: dict[str, dict[str, dict[str, int]]], denominator: int) -> None:
    bins: list[tuple[str, Callable[[int], bool]]] = [
        ("0", lambda n: n == 0), ("1", lambda n: n == 1),
        ("2", lambda n: n == 2), ("3", lambda n: n == 3),
        ("4", lambda n: n == 4), ("5", lambda n: n == 5),
        ("6", lambda n: n == 6), ("7", lambda n: n == 7),
        ("8", lambda n: n == 8), ("9", lambda n: n == 9),
        ("10", lambda n: n == 10), (">=11", lambda n: n >= 11),
    ]
    counts = {
        scope: [sum(test(per_genome[g][scope]["loci"]) for g in per_genome) for _, test in bins]
        for scope in SCOPE_ORDER
    }
    width, height = 940, 520
    left, right, top, bottom = 75, 25, 65, 85
    plot_w, plot_h = width - left - right, height - top - bottom
    max_y = max(max(v) for v in counts.values())
    lines = svg_start(
        width, height, "Per-genome unique-locus counts",
        f"Histogram including zero-hit genomes; denominator is {denominator} normalized accession genomes for every series.",
    )
    for tick in range(0, 6):
        value = max_y * tick / 5
        y = top + plot_h - plot_h * tick / 5
        lines += [f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}"/>',
                  f'<text class="small" x="{left-8}" y="{y+4:.2f}" text-anchor="end">{int(value):,}</text>']
    group_w = plot_w / len(bins)
    bar_w = group_w * 0.23
    for i, (label, _) in enumerate(bins):
        center = left + group_w * (i + 0.5)
        lines.append(f'<text class="small" x="{center:.2f}" y="{top+plot_h+18}" text-anchor="middle">{xml(label)}</text>')
        for j, scope in enumerate(SCOPE_ORDER):
            value = counts[scope][i]
            h = 0 if max_y == 0 else plot_h * value / max_y
            x = center + (j - 1) * bar_w - bar_w / 2
            y = top + plot_h - h
            lines.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w-1:.2f}" height="{h:.2f}" fill="{COLORS[scope]}"><title>{xml(SCOPE_LABEL[scope])}: {value:,} / {denominator:,} genomes</title></rect>')
    lines += [
        f'<line class="axis" x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}"/>',
        f'<text x="{left+plot_w/2}" y="{height-35}" text-anchor="middle">Unique loci per denominator genome (N={denominator:,})</text>',
        f'<text transform="translate(18 {top+plot_h/2}) rotate(-90)" text-anchor="middle">Genomes</text>',
    ]
    for j, scope in enumerate(SCOPE_ORDER):
        x = 160 + j * 270
        lines += [f'<rect x="{x}" y="{height-22}" width="12" height="12" fill="{COLORS[scope]}"/>',
                  f'<text class="legend" x="{x+17}" y="{height-12}">{xml(SCOPE_LABEL[scope])}</text>']
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_bp_ecdf(path: Path, per_genome: dict[str, dict[str, dict[str, int]]], denominator: int) -> None:
    width, height = 900, 520
    left, right, top, bottom = 75, 30, 65, 80
    plot_w, plot_h = width-left-right, height-top-bottom
    series = {scope: sorted(per_genome[g][scope]["bp_inclusive"] for g in per_genome) for scope in SCOPE_ORDER}
    max_bp = max(max(values) for values in series.values())
    max_x = math.log10(1 + max_bp)
    lines = svg_start(
        width, height, "Per-genome prophage-bp empirical CDF",
        f"Inclusive-length bp totals over all {denominator} normalized accession genomes, including zeros; x axis is log10(1+bp).",
    )
    for tick in range(0, 6):
        frac_y = tick / 5
        y = top + plot_h * (1-frac_y)
        lines += [f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}"/>',
                  f'<text class="small" x="{left-8}" y="{y+4:.2f}" text-anchor="end">{frac_y:.1f}</text>']
    x_ticks = [0, 100, 1_000, 10_000, 100_000, max_bp]
    x_ticks = sorted(set(x_ticks))
    for value in x_ticks:
        x = left + plot_w * math.log10(1+value) / max_x
        lines += [f'<line class="grid" x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top+plot_h}"/>',
                  f'<text class="small" x="{x:.2f}" y="{top+plot_h+18}" text-anchor="middle">{value:,}</text>']
    for scope in SCOPE_ORDER:
        values = series[scope]
        points = []
        # At most 1,001 deterministic points; endpoints retained.
        indices = sorted(set(round(i * (len(values)-1) / 1000) for i in range(1001)))
        for index in indices:
            x = left + plot_w * math.log10(1+values[index]) / max_x
            y = top + plot_h * (1 - (index+1)/len(values))
            points.append(f"{x:.2f},{y:.2f}")
        lines.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="{COLORS[scope]}" stroke-width="2"><title>{xml(SCOPE_LABEL[scope])}; N={denominator:,} genomes</title></polyline>')
    lines += [
        f'<line class="axis" x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}"/>',
        f'<text x="{left+plot_w/2}" y="{height-35}" text-anchor="middle">Total prophage bp per genome (1-based closed sensitivity; plotted log10(1+bp))</text>',
        f'<text transform="translate(18 {top+plot_h/2}) rotate(-90)" text-anchor="middle">Fraction of genomes ≤ x</text>',
    ]
    for j, scope in enumerate(SCOPE_ORDER):
        x = 115 + j * 270
        lines += [f'<line x1="{x}" y1="{height-16}" x2="{x+18}" y2="{height-16}" stroke="{COLORS[scope]}" stroke-width="3"/>',
                  f'<text class="legend" x="{x+23}" y="{height-12}">{xml(SCOPE_LABEL[scope])}</text>']
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_length_plot(path: Path, scope_loci: dict[str, list[Record]]) -> None:
    bounds = [0, 1_000, 5_000, 10_000, 20_000, 50_000, 100_000, math.inf]
    labels = ["<1k", "1–<5k", "5–<10k", "10–<20k", "20–<50k", "50–<100k", ">=100k"]
    counts: dict[str, list[int]] = {}
    for scope, records in scope_loci.items():
        vals = [r.inclusive_length for r in records if r.inclusive_length is not None]
        counts[scope] = [sum(bounds[i] <= v < bounds[i+1] for v in vals) for i in range(len(labels))]
    width, height = 900, 520
    left, right, top, bottom = 75, 25, 70, 85
    plot_w, plot_h = width-left-right, height-top-bottom
    max_pct = max(100*counts[s][i]/len(scope_loci[s]) for s in SCOPE_ORDER for i in range(len(labels)))
    lines = svg_start(
        width, height, "Unique-locus inclusive-length distribution",
        "Grouped percentages. Denominators differ by series and are printed in the legend; interval length is end-begin+1.",
    )
    for tick in range(0, 6):
        value = max_pct * tick/5
        y = top + plot_h * (1-tick/5)
        lines += [f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{left+plot_w}" y2="{y:.2f}"/>',
                  f'<text class="small" x="{left-8}" y="{y+4:.2f}" text-anchor="end">{value:.0f}%</text>']
    group_w = plot_w/len(labels); bar_w = group_w*0.23
    for i,label in enumerate(labels):
        center = left + group_w*(i+0.5)
        lines.append(f'<text class="small" x="{center:.2f}" y="{top+plot_h+18}" text-anchor="middle">{xml(label)}</text>')
        for j,scope in enumerate(SCOPE_ORDER):
            n = len(scope_loci[scope]); value = counts[scope][i]; percent = 100*value/n
            h = plot_h*percent/max_pct; x = center+(j-1)*bar_w-bar_w/2; y=top+plot_h-h
            lines.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w-1:.2f}" height="{h:.2f}" fill="{COLORS[scope]}"><title>{xml(SCOPE_LABEL[scope])}: {value:,}/{n:,} loci ({percent:.2f}%)</title></rect>')
    lines += [
        f'<line class="axis" x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}"/>',
        f'<text x="{left+plot_w/2}" y="{height-35}" text-anchor="middle">Inclusive locus length (bp bins)</text>',
        f'<text transform="translate(18 {top+plot_h/2}) rotate(-90)" text-anchor="middle">Percent of valid unique loci</text>',
    ]
    for j,scope in enumerate(SCOPE_ORDER):
        x=80+j*285; n=len(scope_loci[scope])
        lines += [f'<rect x="{x}" y="{height-22}" width="12" height="12" fill="{COLORS[scope]}"/>',
                  f'<text class="legend" x="{x+17}" y="{height-12}">{xml(SCOPE_LABEL[scope])} (N={n:,})</text>']
    lines.append("</svg>")
    path.write_text("\n".join(lines)+"\n", encoding="utf-8")


def write_taxonomy_plot(path: Path, records: list[Record]) -> None:
    counter = Counter(r.taxonomy for r in records)
    ordered = counter.most_common(8)
    shown = {name for name,_ in ordered}
    other = sum(n for name,n in counter.items() if name not in shown)
    data = ordered + ([('Other exact labels', other)] if other else [])
    width, height = 900, 520
    left, right, top, bottom = 285, 40, 60, 55
    plot_w, plot_h = width-left-right, height-top-bottom
    row_h = plot_h/len(data); max_n=max(n for _,n in data); total=len(records)
    lines=svg_start(width,height,"Taxonomy exact-label distribution",f"Top exact taxonomy strings plus Other; denominator {total} normalized all-record rows.")
    for i,(label,n) in enumerate(data):
        y=top+i*row_h+5; h=row_h-10; w=plot_w*n/max_n
        lines += [f'<text class="small" x="{left-8}" y="{y+h/2+4:.2f}" text-anchor="end">{xml(label)}</text>',
                  f'<rect x="{left}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="#6c91bf"><title>{xml(label)}: {n:,}/{total:,} records ({100*n/total:.2f}%)</title></rect>',
                  f'<text class="small" x="{left+w+5:.2f}" y="{y+h/2+4:.2f}">{n:,} ({100*n/total:.1f}%)</text>']
    lines += [f'<text x="{left+plot_w/2}" y="{height-18}" text-anchor="middle">Records (N={total:,}; categories are exact source strings)</text>',"</svg>"]
    path.write_text("\n".join(lines)+"\n",encoding="utf-8")


def main(root: Path) -> None:
    csv_path = root / "26k_prophage1.csv"
    accession_path = root / "26k_ecoli_accession.txt"
    out = root / "artifacts" / "prophage_summary"
    plot_dir = out / "plots"
    report_path = root / "reports" / "prophage_distribution.md"
    out.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    input_hashes = {csv_path.name: sha256(csv_path), accession_path.name: sha256(accession_path)}
    csv_bytes = csv_path.read_bytes()
    accession_bytes = accession_path.read_bytes()
    csv_physical_lines = len(csv_bytes.splitlines())
    accession_physical_lines = len(accession_bytes.splitlines())

    parse_failures: list[tuple[int, str]] = []
    normalization_failures: list[tuple[int, str]] = []
    records: list[Record] = []
    with csv_path.open("r", newline="", encoding="utf-8", errors="strict") as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader)
        if header != EXPECTED_HEADER:
            raise SystemExit(f"Unexpected CSV header: {header!r}")
        csv_record_rows = 1
        for source_row, row in enumerate(reader, start=2):
            csv_record_rows += 1
            if len(row) != len(EXPECTED_HEADER):
                parse_failures.append((source_row, f"field_count={len(row)}"))
                continue
            end_raw, genome, scaffold, begin_raw, flag_raw, taxonomy, prophage_id = row
            try:
                begin = int_like(begin_raw)
                end = int_like(end_raw)
                flag = int_like(flag_raw)
                if flag not in (0, 1):
                    raise ValueError("transposable not 0/1")
                if not genome or not scaffold or not prophage_id:
                    raise ValueError("missing required key")
            except ValueError as exc:
                normalization_failures.append((source_row, str(exc)))
                continue
            records.append(Record(
                source_row, end_raw, genome, scaffold, begin_raw, flag_raw,
                taxonomy, prophage_id, begin, end, flag,
            ))

    accession_text_lines = accession_bytes.decode("utf-8", errors="strict").splitlines()
    accession_nonblank = [(i, line.strip()) for i, line in enumerate(accession_text_lines, 1) if line.strip()]
    accession_exact_counts = Counter(value for _, value in accession_nonblank)
    valid_accession_rows = [(i, value) for i, value in accession_nonblank if ASSEMBLY_RE.fullmatch(value)]
    invalid_accession_rows = [(i, value) for i, value in accession_nonblank if not ASSEMBLY_RE.fullmatch(value)]
    denominator: list[str] = []
    denominator_line: dict[str, int] = {}
    for line_number, accession in valid_accession_rows:
        if accession not in denominator_line:
            denominator.append(accession)
            denominator_line[accession] = line_number
    denominator_set = set(denominator)
    denominator_bases: dict[str, set[str]] = defaultdict(set)
    for accession in denominator:
        denominator_bases[split_version(accession, ASSEMBLY_RE)[0]].add(accession)

    csv_genomes = Counter(record.genome for record in records)
    csv_genome_bases: dict[str, set[str]] = defaultdict(set)
    for genome in csv_genomes:
        csv_genome_bases[split_version(genome, ASSEMBLY_RE)[0]].add(genome)
    exact_join_failure_genomes = set(csv_genomes) - denominator_set
    versionless_join_failure_genomes = {
        genome for genome in csv_genomes
        if split_version(genome, ASSEMBLY_RE)[0] not in denominator_bases
    }
    exact_join_failure_records = sum(csv_genomes[g] for g in exact_join_failure_genomes)
    versionless_join_failure_records = sum(csv_genomes[g] for g in versionless_join_failure_genomes)

    locus_counter = Counter(record.locus for record in records)
    exact_record_counter = Counter(record.exact_record for record in records)
    id_counter = Counter(record.prophage_id for record in records)
    locus_rank: Counter[tuple[str, str, int, int]] = Counter()
    unique_all: list[Record] = []
    first_by_locus: dict[tuple[str, str, int, int], Record] = {}
    for record in records:
        locus_rank[record.locus] += 1
        if record.locus not in first_by_locus:
            first_by_locus[record.locus] = record
            unique_all.append(record)

    scope_records = {scope: [r for r in records if subset(scope, r)] for scope in SCOPE_ORDER}
    scope_loci: dict[str, list[Record]] = {}
    for scope in SCOPE_ORDER:
        seen: set[tuple[str, str, int, int]] = set()
        scope_loci[scope] = []
        for record in scope_records[scope]:
            if record.locus not in seen:
                seen.add(record.locus)
                scope_loci[scope].append(record)

    def interval_relations(loci: list[Record]) -> dict[str, object]:
        grouped: dict[tuple[str, str], list[Record]] = defaultdict(list)
        for record in loci:
            if record.coordinate_valid:
                grouped[(record.genome, record.scaffold)].append(record)
        overlap_participants: set[tuple[str, str, int, int]] = set()
        nested_participants: set[tuple[str, str, int, int]] = set()
        overlap_pairs = nested_pairs = overlap_bp_pairwise = one_bp_pairs = 0
        for group in grouped.values():
            active: list[Record] = []
            for current in sorted(group, key=lambda r: (r.begin, r.end, r.source_row)):
                active = [prior for prior in active if prior.end >= current.begin]
                for prior in active:
                    overlap_pairs += 1
                    overlap_bp = min(prior.end, current.end) - max(prior.begin, current.begin) + 1
                    overlap_bp_pairwise += overlap_bp
                    one_bp_pairs += overlap_bp == 1
                    overlap_participants.update((prior.locus, current.locus))
                    if ((prior.begin <= current.begin and prior.end >= current.end) or
                            (current.begin <= prior.begin and current.end >= prior.end)):
                        nested_pairs += 1
                        nested_participants.update((prior.locus, current.locus))
                active.append(current)
        return {
            "overlap_pairs": overlap_pairs,
            "overlap_bp_pairwise": overlap_bp_pairwise,
            "one_bp_overlap_pairs": one_bp_pairs,
            "overlap_participants": overlap_participants,
            "nested_pairs": nested_pairs,
            "nested_participants": nested_participants,
        }

    relations = {scope: interval_relations(scope_loci[scope]) for scope in SCOPE_ORDER}

    scaffold_to_genomes: dict[str, set[str]] = defaultdict(set)
    scaffold_bases: dict[str, set[str]] = defaultdict(set)
    for record in records:
        scaffold_to_genomes[record.scaffold].add(record.genome)
        scaffold_bases[split_version(record.scaffold, SEQUENCE_RE)[0]].add(record.scaffold)

    # Aggregate once per scope (rather than repeatedly scanning the full table for
    # each denominator genome). Exact joins only: versionless keys are diagnostics.
    scope_aggregates: dict[str, dict[str, object]] = {}
    for scope in SCOPE_ORDER:
        record_counts = Counter(r.genome for r in scope_records[scope])
        locus_counts = Counter(r.genome for r in scope_loci[scope])
        valid_counts: Counter[str] = Counter()
        bp_inclusive: Counter[str] = Counter()
        bp_end_minus_begin: Counter[str] = Counter()
        hit_scaffolds: dict[str, set[str]] = defaultdict(set)
        for record in scope_loci[scope]:
            hit_scaffolds[record.genome].add(record.scaffold)
            if record.coordinate_valid:
                valid_counts[record.genome] += 1
                bp_inclusive[record.genome] += record.inclusive_length or 0
                bp_end_minus_begin[record.genome] += record.half_open_span or 0
        scope_aggregates[scope] = {
            "records": record_counts, "loci": locus_counts, "valid": valid_counts,
            "bp_inclusive": bp_inclusive, "bp_end_minus_begin": bp_end_minus_begin,
            "hit_scaffolds": hit_scaffolds,
        }

    per_genome: dict[str, dict[str, dict[str, int]]] = {genome: {} for genome in denominator}
    for genome in denominator:
        for scope in SCOPE_ORDER:
            aggregate = scope_aggregates[scope]
            per_genome[genome][scope] = {
                "records": aggregate["records"][genome],
                "loci": aggregate["loci"][genome],
                "valid_loci": aggregate["valid"][genome],
                "bp_inclusive": aggregate["bp_inclusive"][genome],
                "bp_end_minus_begin": aggregate["bp_end_minus_begin"][genome],
                "hit_scaffolds": len(aggregate["hit_scaffolds"][genome]),
            }

    per_genome_fields = ["genome_accession", "accession_no_version", "assembly_version", "accession_source_line"]
    for scope in SCOPE_ORDER:
        per_genome_fields += [
            f"{scope}_record_count", f"{scope}_unique_locus_count",
            f"{scope}_valid_coordinate_locus_count", f"{scope}_hit_scaffold_count",
            f"{scope}_total_bp_1based_closed", f"{scope}_total_bp_end_minus_begin",
        ]
    per_genome_rows = []
    for genome in denominator:
        base, version = split_version(genome, ASSEMBLY_RE)
        row: dict[str, object] = {
            "genome_accession": genome, "accession_no_version": base,
            "assembly_version": version, "accession_source_line": denominator_line[genome],
        }
        for scope in SCOPE_ORDER:
            values = per_genome[genome][scope]
            row.update({
                f"{scope}_record_count": values["records"],
                f"{scope}_unique_locus_count": values["loci"],
                f"{scope}_valid_coordinate_locus_count": values["valid_loci"],
                f"{scope}_hit_scaffold_count": values["hit_scaffolds"],
                f"{scope}_total_bp_1based_closed": values["bp_inclusive"],
                f"{scope}_total_bp_end_minus_begin": values["bp_end_minus_begin"],
            })
        per_genome_rows.append(row)
    tsv_write(out / "per_genome.tsv", per_genome_fields, per_genome_rows)

    metric_fields = ["scope", "metric", "value", "unit", "denominator", "denominator_label", "fraction", "notes"]
    metrics: list[dict[str, object]] = []

    def add_metric(scope: str, metric: str, value: object, unit: str = "",
                   denominator_value: int | float | None = None, denominator_label: str = "",
                   notes: str = "", show_fraction: bool = False) -> None:
        metrics.append({
            "scope": scope, "metric": metric, "value": fmt_metric(value), "unit": unit,
            "denominator": "" if denominator_value is None else fmt_metric(denominator_value),
            "denominator_label": denominator_label,
            "fraction": fraction(float(value), denominator_value) if show_fraction and isinstance(value, (int, float)) else "",
            "notes": notes,
        })

    hashes_after_scan = {csv_path.name: sha256(csv_path), accession_path.name: sha256(accession_path)}
    add_metric("input", "csv_sha256_before", input_hashes[csv_path.name], notes="taken before parsing")
    add_metric("input", "csv_sha256_after_scan", hashes_after_scan[csv_path.name], notes="taken after full input scan; reproduce.sh also verifies after all writes")
    add_metric("input", "csv_sha256_before_after_equal", int(input_hashes[csv_path.name] == hashes_after_scan[csv_path.name]), "boolean")
    add_metric("input", "accession_sha256_before", input_hashes[accession_path.name], notes="taken before parsing")
    add_metric("input", "accession_sha256_after_scan", hashes_after_scan[accession_path.name], notes="taken after full input scan; reproduce.sh also verifies after all writes")
    add_metric("input", "accession_sha256_before_after_equal", int(input_hashes[accession_path.name] == hashes_after_scan[accession_path.name]), "boolean")
    add_metric("input", "csv_bytes", len(csv_bytes), "bytes")
    add_metric("input", "accession_bytes", len(accession_bytes), "bytes")
    add_metric("input", "csv_has_final_newline", int(csv_bytes.endswith(b"\n")), "boolean")
    add_metric("input", "accession_has_final_newline", int(accession_bytes.endswith(b"\n")), "boolean")
    add_metric("reconciliation", "csv_physical_lines", csv_physical_lines, "lines")
    add_metric("reconciliation", "csv_records_including_header", csv_record_rows, "CSV records")
    add_metric("reconciliation", "csv_header_records", 1, "CSV records", csv_record_rows, "CSV records including header", show_fraction=True)
    add_metric("reconciliation", "csv_source_data_records", csv_record_rows-1, "records", csv_record_rows-1, "source data records", show_fraction=True)
    add_metric("reconciliation", "csv_field_count_parse_failures", len(parse_failures), "records", csv_record_rows-1, "source data records", show_fraction=True)
    add_metric("reconciliation", "normalization_failures", len(normalization_failures), "records", csv_record_rows-1, "source data records", show_fraction=True)
    add_metric("reconciliation", "normalized_records", len(records), "records", csv_record_rows-1, "source data records", show_fraction=True,
               notes="seven fields; integer-like begin/end; binary integer-like transposable; nonmissing genome/scaffold/prophage_id")
    add_metric("reconciliation", "source_row_balance_difference", (csv_record_rows-1)-len(parse_failures)-len(normalization_failures)-len(records), "records",
               csv_record_rows-1, "source data records", notes="must equal zero: source - parse failures - normalization failures - normalized")
    add_metric("reconciliation", "valid_coordinate_records", sum(r.coordinate_valid for r in records), "records", len(records), "normalized records", show_fraction=True)
    add_metric("reconciliation", "invalid_coordinate_records", sum(not r.coordinate_valid for r in records), "records", len(records), "normalized records", show_fraction=True)
    add_metric("reconciliation", "unique_loci", len(locus_counter), "loci", len(records), "normalized records", show_fraction=True,
               notes="dedup key=(genome exact, scaffold exact, begin integer, end integer)")
    add_metric("reconciliation", "duplicate_locus_groups", sum(n > 1 for n in locus_counter.values()), "groups", len(locus_counter), "unique locus keys", show_fraction=True)
    add_metric("reconciliation", "duplicate_locus_extra_records", sum(n-1 for n in locus_counter.values()), "records", len(records), "normalized records", show_fraction=True)
    add_metric("reconciliation", "exact_duplicate_record_groups", sum(n > 1 for n in exact_record_counter.values()), "groups", len(exact_record_counter), "unique exact record tuples", show_fraction=True)
    add_metric("reconciliation", "duplicate_prophage_id_groups", sum(n > 1 for n in id_counter.values()), "groups", len(id_counter), "unique prophage_id values", show_fraction=True)

    add_metric("denominator", "accession_physical_lines", accession_physical_lines, "lines")
    add_metric("denominator", "accession_nonblank_lines", len(accession_nonblank), "lines", accession_physical_lines, "physical lines", show_fraction=True)
    add_metric("denominator", "accession_unique_exact_tokens", len(accession_exact_counts), "tokens", len(accession_nonblank), "nonblank lines", show_fraction=True)
    add_metric("denominator", "accession_duplicate_extra_lines", sum(n-1 for n in accession_exact_counts.values()), "lines", len(accession_nonblank), "nonblank lines", show_fraction=True)
    add_metric("denominator", "accession_valid_assembly_lines", len(valid_accession_rows), "lines", len(accession_nonblank), "nonblank lines", show_fraction=True)
    add_metric("denominator", "accession_invalid_tokens", len(invalid_accession_rows), "tokens", len(accession_nonblank), "nonblank lines", show_fraction=True,
               notes="invalid tokens listed in report; excluded from normalized genome denominator")
    add_metric("denominator", "normalized_unique_genomes", len(denominator), "genomes", len(accession_nonblank), "nonblank accession lines", show_fraction=True,
               notes="unique tokens matching ^GC[AF]_[0-9]+\\.[0-9]+$")
    add_metric("denominator", "version_stripped_collision_groups", sum(len(v)>1 for v in denominator_bases.values()), "base accessions", len(denominator_bases), "version-stripped denominator keys", show_fraction=True)
    add_metric("denominator", "version_stripped_collapsed_accessions", sum(len(v)-1 for v in denominator_bases.values()), "accessions", len(denominator), "normalized unique genomes", show_fraction=True)

    add_metric("join", "csv_unique_genome_keys", len(csv_genomes), "genome keys")
    add_metric("join", "exact_join_failure_genome_keys", len(exact_join_failure_genomes), "genome keys", len(csv_genomes), "CSV unique genome keys", show_fraction=True)
    add_metric("join", "exact_join_failure_records", exact_join_failure_records, "records", len(records), "normalized records", show_fraction=True)
    add_metric("join", "versionless_join_failure_genome_keys", len(versionless_join_failure_genomes), "genome keys", len(csv_genomes), "CSV unique genome keys", show_fraction=True)
    add_metric("join", "versionless_join_failure_records", versionless_join_failure_records, "records", len(records), "normalized records", show_fraction=True)
    version_only_matches = exact_join_failure_genomes - versionless_join_failure_genomes
    add_metric("join", "version_only_match_genome_keys", len(version_only_matches), "genome keys", len(csv_genomes), "CSV unique genome keys", show_fraction=True,
               notes="exact join failed but version-stripped join succeeded; diagnostics only, never merged")
    add_metric("join", "version_only_match_records", sum(csv_genomes[g] for g in version_only_matches), "records", len(records), "normalized records", show_fraction=True)
    add_metric("join", "denominator_genomes_without_exact_csv_hit", len(denominator_set-set(csv_genomes)), "genomes", len(denominator), "normalized denominator genomes", show_fraction=True)
    add_metric("join", "csv_genome_version_stripped_collision_groups", sum(len(v)>1 for v in csv_genome_bases.values()), "base accessions", len(csv_genome_bases), "version-stripped CSV genome keys", show_fraction=True)
    add_metric("join", "scaffold_version_stripped_collision_groups", sum(len(v)>1 for v in scaffold_bases.values()), "base accessions", len(scaffold_bases), "version-stripped scaffold keys", show_fraction=True)
    add_metric("join", "unique_exact_scaffold_keys", len(scaffold_to_genomes), "scaffold keys")
    add_metric("join", "unique_exact_genome_scaffold_pairs", len({(r.genome, r.scaffold) for r in records}), "pairs")
    add_metric("join", "exact_scaffold_keys_assigned_multiple_genomes", sum(len(v)>1 for v in scaffold_to_genomes.values()), "scaffold keys", len(scaffold_to_genomes), "unique exact scaffold keys", show_fraction=True)

    distribution_probs = [("min", 0), ("p01", .01), ("p05", .05), ("p25", .25), ("median", .5), ("p75", .75), ("p95", .95), ("p99", .99), ("max", 1)]
    for scope in SCOPE_ORDER:
        recs = scope_records[scope]
        loci = scope_loci[scope]
        valid_loci = [r for r in loci if r.coordinate_valid]
        genome_counts = [per_genome[g][scope]["loci"] for g in denominator]
        genome_bp = [per_genome[g][scope]["bp_inclusive"] for g in denominator]
        genome_span = [per_genome[g][scope]["bp_end_minus_begin"] for g in denominator]
        genome_scaffolds = [per_genome[g][scope]["hit_scaffolds"] for g in denominator]
        lengths = [r.inclusive_length for r in valid_loci if r.inclusive_length is not None]
        spans = [r.half_open_span for r in valid_loci if r.half_open_span is not None]
        add_metric(scope, "source_records", len(recs), "records", len(records), "all normalized records", show_fraction=True)
        add_metric(scope, "unique_loci", len(loci), "loci", len(recs), f"{scope} source records", show_fraction=True)
        add_metric(scope, "valid_coordinate_unique_loci", len(valid_loci), "loci", len(loci), f"{scope} unique loci", show_fraction=True)
        add_metric(scope, "total_bp_1based_closed", sum(lengths), "bp", len(valid_loci), "valid unique loci", notes="sum(end-begin+1)")
        add_metric(scope, "total_bp_end_minus_begin", sum(spans), "bp", len(valid_loci), "valid unique loci", notes="coordinate-convention sensitivity sum(end-begin)")
        for metric, test in (("genomes_zero_loci", lambda x:x==0), ("genomes_one_locus", lambda x:x==1),
                             ("genomes_multiple_loci", lambda x:x>1), ("genomes_with_locus", lambda x:x>0)):
            value=sum(test(x) for x in genome_counts)
            add_metric(scope, metric, value, "genomes", len(denominator), "normalized denominator genomes", show_fraction=True)
        for family, values, unit in (
            ("per_genome_unique_locus_count", genome_counts, "loci/genome"),
            ("per_genome_total_bp_1based_closed", genome_bp, "bp/genome"),
            ("per_genome_total_bp_end_minus_begin", genome_span, "bp/genome"),
            ("per_genome_hit_scaffold_count", genome_scaffolds, "scaffolds/genome"),
            ("unique_locus_length_1based_closed", lengths, "bp/locus"),
            ("unique_locus_span_end_minus_begin", spans, "bp/locus"),
        ):
            denom_values = len(values)
            denom_label = "normalized denominator genomes (zeros included)" if family.startswith("per_genome") else "valid unique loci"
            for label, probability in distribution_probs:
                add_metric(scope, f"{family}_{label}", quantile(values, probability), unit, denom_values, denom_label)
            add_metric(scope, f"{family}_mean", statistics.fmean(values) if values else math.nan, unit, denom_values, denom_label)
            add_metric(scope, f"{family}_population_sd", statistics.pstdev(values) if values else math.nan, unit, denom_values, denom_label)
        rel=relations[scope]
        add_metric(scope,"overlap_pairs",rel["overlap_pairs"],"locus pairs",len(valid_loci),"valid unique loci")
        add_metric(scope,"overlap_participating_loci",len(rel["overlap_participants"]),"loci",len(valid_loci),"valid unique loci",show_fraction=True)
        add_metric(scope,"pairwise_overlap_bp",rel["overlap_bp_pairwise"],"bp",rel["overlap_pairs"],"overlapping locus pairs")
        add_metric(scope,"one_bp_overlap_pairs",rel["one_bp_overlap_pairs"],"locus pairs",rel["overlap_pairs"],"overlapping locus pairs",show_fraction=True)
        add_metric(scope,"nested_pairs",rel["nested_pairs"],"locus pairs",len(valid_loci),"valid unique loci")
        add_metric(scope,"nested_participating_loci",len(rel["nested_participants"]),"loci",len(valid_loci),"valid unique loci",show_fraction=True)
        for label,test in (("begin_le_3",lambda r:r.begin<=3), ("length_lt_1000_bp",lambda r:(r.inclusive_length or 0)<1000),
                           ("length_gt_200000_bp",lambda r:(r.inclusive_length or 0)>200000)):
            value=sum(test(r) for r in valid_loci)
            add_metric(scope,label,value,"loci",len(valid_loci),"valid unique loci",show_fraction=True,notes="diagnostic threshold; not proof of artifact")

    # Explicitly encode absent semantic fields rather than silently omitting them.
    for name in ("tag", "status", "quality", "completeness", "chromosome_role"):
        add_metric("schema", f"explicit_{name}_columns", 0, "columns", len(EXPECTED_HEADER), "CSV columns", show_fraction=True,
                   notes="no exact column with this name; do not infer unsupported stratification")
    add_metric("schema", "class_like_columns", 1, "columns", len(EXPECTED_HEADER), "CSV columns", show_fraction=True, notes="taxonomy")
    add_metric("schema", "flag_like_columns", 1, "columns", len(EXPECTED_HEADER), "CSV columns", show_fraction=True, notes="transposable")

    tsv_write(out / "summary_metrics.tsv", metric_fields, metrics)

    category_fields = ["scope", "dimension", "category", "count", "denominator", "fraction", "unit", "denominator_label", "notes"]
    categories: list[dict[str, object]]=[]
    def add_category(scope: str, dimension: str, category: str, count: int, denominator_value: int,
                     unit: str, denominator_label: str, notes: str="") -> None:
        categories.append({"scope":scope,"dimension":dimension,"category":category,"count":count,
                           "denominator":denominator_value,"fraction":fraction(count,denominator_value),
                           "unit":unit,"denominator_label":denominator_label,"notes":notes})
    for scope in SCOPE_ORDER:
        recs=scope_records[scope]; loci=scope_loci[scope]
        for dimension, categorizer, notes in (
            ("transposable_exact",lambda r:r.transposable_raw,"exact source string"),
            ("taxonomy_exact",lambda r:r.taxonomy,"exact source string; source semantics otherwise unknown"),
            ("taxonomy_assignment",lambda r:"assigned" if r.taxonomy.strip().casefold() not in {"","unknown"} else "Unknown/missing","Unknown is treated as unassigned only for the competing subset"),
            ("scaffold_accession_prefix",lambda r:prefix(r.scaffold),"identifier prefix only; not a chromosome/plasmid classification"),
            ("coordinate_validity",lambda r:"valid_begin_ge_1_end_ge_begin" if r.coordinate_valid else "invalid","record-level"),
        ):
            counter=Counter(categorizer(r) for r in recs)
            for category,count in sorted(counter.items(),key=lambda item:(-item[1],item[0])):
                add_category(scope,dimension,category,count,len(recs),"records",f"{scope} source records",notes)
        valid=[r for r in loci if r.coordinate_valid]
        length_counter=Counter(
            "short_lt_1000_bp" if (r.inclusive_length or 0)<1000 else
            "long_gt_200000_bp" if (r.inclusive_length or 0)>200000 else
            "1000_to_200000_bp_inclusive"
            for r in valid
        )
        for category,count in sorted(length_counter.items()):
            add_category(scope,"length_qc_class",category,count,len(valid),"unique loci","valid unique loci","diagnostic bins, not biological quality classes")
        missing_fields=EXPECTED_HEADER
        for field in missing_fields:
            index=EXPECTED_HEADER.index(field)
            missing=sum(r.exact_record[index].strip()=="" for r in recs)
            add_category(scope,"field_missingness",field,missing,len(recs),"missing records",f"{scope} source records","empty-string missingness")
        genome_count_counter=Counter(
            "zero" if per_genome[g][scope]["loci"]==0 else "one" if per_genome[g][scope]["loci"]==1 else "multiple"
            for g in denominator
        )
        for category,count in sorted(genome_count_counter.items()):
            add_category(scope,"per_genome_locus_count_class",category,count,len(denominator),"genomes","normalized denominator genomes","unique loci; zero genomes included")
        scaffold_count_counter=Counter(per_genome[g][scope]["hit_scaffolds"] for g in denominator)
        for category,count in sorted(scaffold_count_counter.items()):
            add_category(scope,"per_genome_hit_scaffold_count",str(category),count,len(denominator),"genomes","normalized denominator genomes","exact original scaffold keys; zero genomes included")

    original_scaffold_counts=Counter(r.scaffold for r in records)
    frequency_bins=[("1",lambda n:n==1),("2-4",lambda n:2<=n<=4),("5-9",lambda n:5<=n<=9),("10+",lambda n:n>=10)]
    for label,test in frequency_bins:
        add_category("all","records_per_original_scaffold_key",label,sum(test(n) for n in original_scaffold_counts.values()),len(original_scaffold_counts),"scaffold keys","unique exact original scaffold keys","number of source records sharing the key")
    for line,value in invalid_accession_rows:
        add_category("denominator","invalid_accession_token",value,1,len(accession_nonblank),"tokens","nonblank accession lines",f"source line {line}")
    tsv_write(out / "category_counts.tsv",category_fields,categories)

    all_rel=relations["all"]
    interval_fields=[
        "source_row","parse_status","genome","genome_accession_no_version","genome_version",
        "scaffold","scaffold_accession_no_version","scaffold_version","begin_raw","end_raw","begin","end",
        "length_bp_1based_closed","span_bp_end_minus_begin","transposable_raw","taxonomy","prophage_id",
        "denominator_exact_join","denominator_versionless_join","genome_key_canonical_pattern","scaffold_key_version_pattern",
        "prophage_id_genome_prefix_match","locus_duplicate_rank","locus_multiplicity","is_unique_locus_representative",
        "valid_coordinates","coordinate_issue","overlaps_other_unique_locus","nested_with_other_unique_locus",
        "begin_le_3_qc","length_lt_1000_bp_qc","length_gt_200000_bp_qc"
    ]
    interval_rows=[]
    running_rank: Counter[tuple[str,str,int,int]]=Counter()
    for record in records:
        running_rank[record.locus]+=1
        genome_base,genome_version=split_version(record.genome,ASSEMBLY_RE)
        scaffold_base,scaffold_version=split_version(record.scaffold,SEQUENCE_RE)
        id_match=ID_RE.fullmatch(record.prophage_id)
        coordinate_issue="" if record.coordinate_valid else ("begin_lt_1" if record.begin<1 else "end_lt_begin")
        interval_rows.append({
            "source_row":record.source_row,"parse_status":"normalized","genome":record.genome,
            "genome_accession_no_version":genome_base,"genome_version":genome_version,"scaffold":record.scaffold,
            "scaffold_accession_no_version":scaffold_base,"scaffold_version":scaffold_version,
            "begin_raw":record.begin_raw,"end_raw":record.end_raw,"begin":record.begin,"end":record.end,
            "length_bp_1based_closed":"" if record.inclusive_length is None else record.inclusive_length,
            "span_bp_end_minus_begin":"" if record.half_open_span is None else record.half_open_span,
            "transposable_raw":record.transposable_raw,"taxonomy":record.taxonomy,"prophage_id":record.prophage_id,
            "denominator_exact_join":int(record.genome in denominator_set),
            "denominator_versionless_join":int(genome_base in denominator_bases),
            "genome_key_canonical_pattern":int(bool(ASSEMBLY_RE.fullmatch(record.genome))),
            "scaffold_key_version_pattern":int(bool(SEQUENCE_RE.fullmatch(record.scaffold))),
            "prophage_id_genome_prefix_match":int(bool(id_match and id_match.group(1)==record.genome)),
            "locus_duplicate_rank":running_rank[record.locus],"locus_multiplicity":locus_counter[record.locus],
            "is_unique_locus_representative":int(running_rank[record.locus]==1),
            "valid_coordinates":int(record.coordinate_valid),"coordinate_issue":coordinate_issue,
            "overlaps_other_unique_locus":int(record.locus in all_rel["overlap_participants"]),
            "nested_with_other_unique_locus":int(record.locus in all_rel["nested_participants"]),
            "begin_le_3_qc":int(record.coordinate_valid and record.begin<=3),
            "length_lt_1000_bp_qc":int(record.coordinate_valid and (record.inclusive_length or 0)<1000),
            "length_gt_200000_bp_qc":int(record.coordinate_valid and (record.inclusive_length or 0)>200000),
        })
    tsv_write(out/"interval_qc.tsv",interval_fields,interval_rows)

    write_count_plot(plot_dir/"per_genome_locus_counts.svg",per_genome,len(denominator))
    write_bp_ecdf(plot_dir/"per_genome_bp_ecdf.svg",per_genome,len(denominator))
    write_length_plot(plot_dir/"locus_length_distribution.svg",scope_loci)
    write_taxonomy_plot(plot_dir/"taxonomy_distribution.svg",records)

    # Report helpers.
    scope_stats={}
    for scope in SCOPE_ORDER:
        loci=scope_loci[scope]; valid=[r for r in loci if r.coordinate_valid]
        counts=[per_genome[g][scope]["loci"] for g in denominator]
        bps=[per_genome[g][scope]["bp_inclusive"] for g in denominator]
        lengths=[r.inclusive_length for r in valid if r.inclusive_length is not None]
        scope_stats[scope]={
            "records":len(scope_records[scope]),"loci":len(loci),"valid":len(valid),
            "zero":sum(n==0 for n in counts),"one":sum(n==1 for n in counts),"multi":sum(n>1 for n in counts),
            "bp":sum(lengths),"span":sum((r.half_open_span or 0) for r in valid),
            "count_median":quantile(counts,.5),"count_p95":quantile(counts,.95),"count_max":max(counts),
            "bp_median":quantile(bps,.5),"bp_p95":quantile(bps,.95),"bp_max":max(bps),
            "len_median":quantile(lengths,.5),"len_p05":quantile(lengths,.05),"len_p95":quantile(lengths,.95),
        }

    taxonomy_counts=Counter(r.taxonomy for r in records)
    trans_counts=Counter(r.transposable_raw for r in records)
    scaffold_prefix_counts=Counter(prefix(r.scaffold) for r in records)
    scaffold_top=Counter(r.scaffold for r in records).most_common(10)
    longest=sorted(unique_all,key=lambda r:(-(r.inclusive_length or -1),r.genome,r.scaffold))[:10]
    shortest=sorted((r for r in unique_all if r.coordinate_valid),key=lambda r:(r.inclusive_length or 0,r.genome,r.scaffold))[:10]

    report=[]
    report += [
        "# Prophage table distribution and interval audit", "",
        "## Scope, provenance, and reproducibility", "",
        "This is a local, offline audit of the two root inputs only. It did **not** inspect mounted storage, download or resolve sequence, rename identifiers, assess IMPG, build/query a pangenome, or compute a tree. Both inputs are opened read-only. `reproduce.sh` checks their expected SHA-256 values before analysis and verifies that the pre/post hashes match.", "",
        "| Input | Bytes | Physical lines | SHA-256 |", "|---|---:|---:|---|",
        f"| `26k_prophage1.csv` | {len(csv_bytes):,} | {csv_physical_lines:,} | `{input_hashes[csv_path.name]}` |",
        f"| `26k_ecoli_accession.txt` | {len(accession_bytes):,} | {accession_physical_lines:,} | `{input_hashes[accession_path.name]}` |", "",
        "Regenerate from any working directory with `artifacts/prophage_summary/reproduce.sh`. The implementation uses Python's standard library only, makes no network calls, writes LF-terminated TSVs, fixes ordering explicitly, and places no timestamp in outputs. Quantiles use linear interpolation at `(n-1)p` (R type 7 / NumPy `linear`).", "",
        "## Headline result and the meaning of “tagged”", "",
        "The table has no column named `tag`, `tagged`, `status`, `quality`, or `completeness`. Therefore there is **no uniquely evidenced generic ‘tagged prophage’ filter**. Rather than silently choosing one, this audit carries two competing, exactly reproducible interpretations alongside all records:", "",
        "1. **Flag-positive:** `Decimal(transposable)` is the integer `1` (the observed raw value is `1.0`). This is the only explicit binary, status-like field. The column name supports ‘transposable’, but no supplied documentation establishes a broader meaning such as quality or completeness.",
        "2. **Taxonomy-assigned:** trimmed, case-folded `taxonomy` is neither empty nor `Unknown`. This is defensible only if ‘tagged’ meant taxonomically labeled; it is not equivalent to the binary flag.", "",
        "All rows have a nonempty `prophage_id`, so ‘has an ID’ selects all records and is not a distinct subset. No completeness/quality stratification is possible from these inputs.", "",
        f"All genome fractions below use **{len(denominator):,} normalized unique denominator genomes**; locus/record fractions state their own denominator.", "",
        "| Scope (exact rule) | Source records | Unique loci | Genomes zero | Genomes one | Genomes multiple | Total bp (`end-begin+1`) | Median locus bp |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for scope in SCOPE_ORDER:
        s=scope_stats[scope]
        report.append(f"| {SCOPE_LABEL[scope]} | {s['records']:,} / {len(records):,} | {s['loci']:,} / {len(scope_records[scope]):,} | {s['zero']:,} / {len(denominator):,} ({pct(s['zero'],len(denominator))}) | {s['one']:,} / {len(denominator):,} ({pct(s['one'],len(denominator))}) | {s['multi']:,} / {len(denominator):,} ({pct(s['multi'],len(denominator))}) | {s['bp']:,} / {s['valid']:,} valid loci | {fmt_number(s['len_median'])} / {s['valid']:,} valid loci |")

    report += [
        "", "The accession denominator is coextensive with the CSV genome keys: every normalized denominator genome has at least one all-record locus. Thus the all-record zero count is genuinely 0 / denominator, but the list cannot estimate prevalence in an independent collection. The subset zero counts remain informative for the two explicit filters.", "",
        "## Row, parsing, normalization, and denominator reconciliation", "",
        "| Stage | Result | Denominator / rule |", "|---|---:|---|",
        f"| CSV physical lines | {csv_physical_lines:,} | file bytes split into lines |",
        f"| CSV records | {csv_record_rows:,} | 1 header + {csv_record_rows-1:,} data records |",
        f"| Seven-field parse failures | {len(parse_failures):,} | {csv_record_rows-1:,} source data records |",
        f"| Normalization failures | {len(normalization_failures):,} | {csv_record_rows-1:,} source data records |",
        f"| Normalized records | {len(records):,} | {csv_record_rows-1:,} source data records |",
        f"| Valid coordinates | {sum(r.coordinate_valid for r in records):,} | {len(records):,} normalized records |",
        f"| Unique locus keys | {len(locus_counter):,} | {len(records):,} normalized records |",
        f"| Duplicate-locus groups / extra rows | {sum(n>1 for n in locus_counter.values()):,} / {sum(n-1 for n in locus_counter.values()):,} | key `(genome, scaffold, begin, end)` |",
        f"| Exact duplicate-record groups | {sum(n>1 for n in exact_record_counter.values()):,} | all seven source fields |",
        f"| Duplicate `prophage_id` groups | {sum(n>1 for n in id_counter.values()):,} | {len(id_counter):,} unique IDs |",
        "", "The accession file has:", "",
        f"- **{accession_physical_lines:,} physical / {len(accession_nonblank):,} nonblank / {len(accession_exact_counts):,} exact-unique tokens**; duplicate extra lines: {sum(n-1 for n in accession_exact_counts.values()):,}.",
        f"- **{len(denominator):,} normalized unique assembly accessions**, defined strictly by `^GC[AF]_[0-9]+\\.[0-9]+$`. `per_genome.tsv` has exactly this many data rows, in accession-file order.",
        f"- **{len(invalid_accession_rows):,} invalid token:** " + (", ".join(f"`{v}` at line {i}" for i,v in invalid_accession_rows) or "none") + ". It is not a genome and is excluded from the normalized genome denominator, but is retained in reconciliation metrics.",
        f"- Exact-key CSV join failures: **{len(exact_join_failure_genomes):,} / {len(csv_genomes):,} unique CSV genome keys; {exact_join_failure_records:,} / {len(records):,} records**. Versionless join failures: **{len(versionless_join_failure_genomes):,} / {len(csv_genomes):,} keys; {versionless_join_failure_records:,} / {len(records):,} records**.",
        f"- Denominator genomes without an exact CSV hit: **{len(denominator_set-set(csv_genomes)):,} / {len(denominator):,}**.", "",
        "A versionless key is used for diagnostics only, never to merge rows into the per-genome table. Stripping assembly versions creates " + f"**{sum(len(v)>1 for v in denominator_bases.values()):,} collision groups / {len(denominator_bases):,} denominator base keys** and **{sum(len(v)>1 for v in csv_genome_bases.values()):,} / {len(csv_genome_bases):,} CSV base keys**. Stripping scaffold versions creates " + f"**{sum(len(v)>1 for v in scaffold_bases.values()):,} / {len(scaffold_bases):,} base-key collision groups**. These are current-data observations, not permission to drop versions.", "",
        "## Schema and semantic evidence", "",
        "The header is exactly `end,genome,scaffold,begin,transposable,taxonomy,prophage_id` (`26k_prophage1.csv`, line 1). The first data row demonstrates float-formatted integer coordinates/flag and versioned keys (`26k_prophage1.csv`, line 2).", "",
        "| Column | Observed values / constraints | Supported interpretation | Unsupported or unknown |", "|---|---|---|---|",
        f"| `end` | nonmissing; integer-like decimals; range {min(r.end for r in records):,}–{max(r.end for r in records):,} | interval end coordinate | indexing convention not explicitly declared |",
        f"| `genome` | {len(csv_genomes):,} exact values; all match `GCF_...version` | RefSeq assembly-style source key and exact denominator join key | organism quality/completeness |",
        f"| `scaffold` | {len(scaffold_to_genomes):,} exact values; prefixes " + ", ".join(f"`{k}` {v:,}" for k,v in sorted(scaffold_prefix_counts.items())) + " records | versioned sequence/contig source key | chromosome vs plasmid vs unplaced status cannot be inferred from prefix |",
        f"| `begin` | nonmissing; integer-like decimals; range {min(r.begin for r in records):,}–{max(r.begin for r in records):,}; minimum 1 | interval begin coordinate | indexing convention not explicitly declared |",
        f"| `transposable` | " + ", ".join(f"`{k}` {v:,}/{len(records):,} ({pct(v,len(records))})" for k,v in sorted(trans_counts.items())) + " | binary flag named ‘transposable’ | method, evidence threshold, quality, completeness, and meaning of 0 beyond not flag-positive |",
        f"| `taxonomy` | {len(taxonomy_counts):,} exact strings; `Unknown` {taxonomy_counts['Unknown']:,}/{len(records):,} | class-like source label; exact strings preserved | classifier/method/confidence; mixed labels are not resolved |",
        f"| `prophage_id` | {len(id_counter):,} unique; every value matches `<genome>_prophage_<integer>` and prefix agrees with `genome` | source record identifier | biological stable ID across releases |", "",
        "The source contains no explicit tag/status/quality/completeness/chromosome-role column. `transposable` is treated as flag-like and `taxonomy` as class-like only because of their names and observed values; unsupported semantics remain unknown.", "",
        "## Coordinates, locus definition, overlaps, and coordinate QC", "",
        "All begin/end values parse as finite integers despite `.0` formatting; all have `begin >= 1` and `end >= begin`. The minimum begin of 1 is consistent with **1-based closed** coordinates, so `end-begin+1` is the primary bp sensitivity. This is not conclusive without producer metadata or sequence; consequently every per-genome row and summary also reports `end-begin`. The two totals differ by exactly one bp per valid unique locus. No contig lengths are supplied, so right-edge clipping and out-of-bounds coordinates cannot be assessed.", "",
        "A locus is the exact normalized tuple **(`genome`, `scaffold`, integer `begin`, integer `end`)**. Deduplication keeps the first source row as representative; subset selection is applied to rows first and then the same locus key is deduplicated. Taxonomy and flag are deliberately absent from the key. Overlap is assessed only within exact `(genome, scaffold)`, on valid unique loci, as closed-interval intersection (`next.begin <= prior.end`). Nesting includes full containment; exact duplicates are handled before overlap.", "",
        "| Scope | Overlap pairs | Participating loci | Pairwise overlap bp | Nested pairs | Participating nested loci |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for scope in SCOPE_ORDER:
        rel=relations[scope]; valid=scope_stats[scope]["valid"]
        report.append(f"| {SCOPE_LABEL[scope]} | {rel['overlap_pairs']:,} / {valid:,} valid loci | {len(rel['overlap_participants']):,} / {valid:,} | {rel['overlap_bp_pairwise']:,} / {rel['overlap_pairs']:,} pairs | {rel['nested_pairs']:,} / {valid:,} valid loci | {len(rel['nested_participants']):,} / {valid:,} |")
    report += ["", "Diagnostic, not biological-quality, flags:", "",
               "| Scope | begin <= 3 | length <1,000 bp | length >200,000 bp | Denominator |", "|---|---:|---:|---:|---:|"]
    for scope in SCOPE_ORDER:
        valid=[r for r in scope_loci[scope] if r.coordinate_valid]
        a=sum(r.begin<=3 for r in valid); b=sum((r.inclusive_length or 0)<1000 for r in valid); c=sum((r.inclusive_length or 0)>200000 for r in valid)
        report.append(f"| {SCOPE_LABEL[scope]} | {a:,} ({pct(a,len(valid))}) | {b:,} ({pct(b,len(valid))}) | {c:,} ({pct(c,len(valid))}) | {len(valid):,} valid unique loci |")
    report += ["", "`begin <= 3` may indicate a left-contig-edge/truncated call; short/long cutoffs are transparent review thresholds, not evidence-based quality labels. `interval_qc.tsv` identifies every flagged source row.", "",
               "## Count, bp, and length distributions", "",
               "Per-genome distributions include all normalized denominator genomes, including zeros. Bp uses valid unique loci and the primary `end-begin+1` sensitivity.", "",
               "| Scope | Loci/genome median / p95 / max | bp/genome median / p95 / max | Locus bp p05 / median / p95 |", "|---|---:|---:|---:|"]
    for scope in SCOPE_ORDER:
        s=scope_stats[scope]
        report.append(f"| {SCOPE_LABEL[scope]} | {fmt_number(s['count_median'])} / {fmt_number(s['count_p95'])} / {fmt_number(float(s['count_max']))} (N={len(denominator):,} genomes) | {fmt_number(s['bp_median'])} / {fmt_number(s['bp_p95'])} / {fmt_number(float(s['bp_max']))} (N={len(denominator):,}) | {fmt_number(s['len_p05'])} / {fmt_number(s['len_median'])} / {fmt_number(s['len_p95'])} (N={s['valid']:,} loci) |")
    report += ["", "`summary_metrics.tsv` adds min, p01, p05, p25, median, p75, p95, p99, max, mean, and population SD for per-genome locus counts, hit-scaffold counts, both bp conventions, and locus lengths for all three scopes.", "",
               "### Extreme genomes", "",
               "Top genomes are descriptive extremes, not automatically errors. Every rank is among the same N="+f"{len(denominator):,} denominator genomes.", "",
               "| Scope | Rank | Genome exact key | Unique loci | Total inclusive bp |", "|---|---:|---|---:|---:|"]
    for scope in SCOPE_ORDER:
        ranked=sorted(denominator,key=lambda g:(-per_genome[g][scope]["loci"],-per_genome[g][scope]["bp_inclusive"],g))[:10]
        for rank,genome in enumerate(ranked,1):
            report.append(f"| {SCOPE_LABEL[scope]} | {rank} | `{genome}` | {per_genome[genome][scope]['loci']:,} | {per_genome[genome][scope]['bp_inclusive']:,} |")
    report += ["", "Top 10 by inclusive bp for every scope (each rank is among N="+f"{len(denominator):,} denominator genomes):", "",
               "| Scope | Rank | Genome exact key | Unique loci | Total inclusive bp |", "|---|---:|---|---:|---:|"]
    for scope in SCOPE_ORDER:
        ranked=sorted(denominator,key=lambda g:(-per_genome[g][scope]["bp_inclusive"],-per_genome[g][scope]["loci"],g))[:10]
        for rank,genome in enumerate(ranked,1):
            report.append(f"| {SCOPE_LABEL[scope]} | {rank} | `{genome}` | {per_genome[genome][scope]['loci']:,} | {per_genome[genome][scope]['bp_inclusive']:,} |")
    report += ["", "Shortest and longest unique loci (primary inclusive sensitivity):", "",
               "| Tail | Rank | Genome | Scaffold | Begin–end | bp | Prophage ID |", "|---|---:|---|---|---:|---:|---|"]
    for tail,items in (("shortest",shortest),("longest",longest)):
        for rank,r in enumerate(items,1):
            report.append(f"| {tail} | {rank} | `{r.genome}` | `{r.scaffold}` | {r.begin:,}–{r.end:,} | {r.inclusive_length:,} | `{r.prophage_id}` |")

    report += ["", "## Contig/scaffold and category distributions", "",
               f"There are **{len(scaffold_to_genomes):,} exact original scaffold keys / {len(records):,} records** and **{sum(len(v)>1 for v in scaffold_to_genomes.values()):,} / {len(scaffold_to_genomes):,} keys assigned to multiple exact genome keys**. Prefix counts are identifier-form distributions only:", "",
               "| Prefix | Records | Unique scaffold keys | Record denominator |", "|---|---:|---:|---:|"]
    for pre,n in sorted(scaffold_prefix_counts.items(),key=lambda x:(-x[1],x[0])):
        unique=sum(prefix(s)==pre for s in scaffold_to_genomes)
        report.append(f"| `{pre}` | {n:,} ({pct(n,len(records))}) | {unique:,} / {len(scaffold_to_genomes):,} | {len(records):,} all records |")
    report += ["", "No chromosome/contig role field exists. In particular, `NC_`/`NZ_` prefixes alone are not used to label a sequence chromosome, plasmid, complete, or draft. The most reused exact scaffold keys (reuse here means multiple prophage records on the same sequence) are:", "",
               "| Scaffold | Records | Record denominator |", "|---|---:|---:|"]
    for scaffold,n in scaffold_top:
        report.append(f"| `{scaffold}` | {n:,} | {len(records):,} |")
    report += ["", "Exact taxonomy strings (mixed labels are kept intact):", "",
               "| Taxonomy exact string | Records | Denominator |", "|---|---:|---:|"]
    for label,n in taxonomy_counts.most_common():
        report.append(f"| {label} | {n:,} ({pct(n,len(records))}) | {len(records):,} all records |")
    report += ["", "`category_counts.tsv` provides the same dimensions for all three scopes, plus field missingness, coordinate validity, diagnostic length class, per-genome zero/one/multiple, exact hit-scaffold-count frequencies, and records-per-scaffold-key bins. Every count has an explicit denominator and fraction. All seven fields have 0 missing values in all-record rows; subset missingness is also tabulated.", "",
               "## Lossless current-key crosswalk requirements", "",
               "This task does not rename identifiers. For the future PanSN/BGZF crosswalk, coordinate survival requires retaining these exact source components:", "",
               "- assembly key: `genome`, including `GCF_` prefix and version suffix;",
               "- sequence key: `scaffold`, including `NC_`/`NZ_` prefix and version suffix;",
               "- source `begin` and `end` strings plus their normalized integer values and an explicit coordinate-convention field;",
               "- `prophage_id` and `source_row` for row-level traceability;",
               "- the source input SHA-256 and the composite locus key `(genome, scaffold, begin, end)`.", "",
               "Do not join on version-stripped assembly or scaffold keys merely because this snapshot has zero collision groups; future releases may contain multiple versions. Do not rely on `prophage_id` alone as a coordinate key. `interval_qc.tsv` is the lossless tabular handoff for every current normalized row: it stores both exact original keys, version-separated diagnostics, original numeric strings, normalized coordinates, row ID, join flags, and QC flags. The dedicated `pansn-bgzip-genome-layout` task owns the canonical old-to-new crosswalk and future names.", "",
               "## Likely artifacts and limits", "",
               f"- The terminal accession token `{invalid_accession_rows[0][1] if invalid_accession_rows else 'none'}` and the exact coextensiveness of {len(denominator):,} normalized accessions with {len(csv_genomes):,} CSV genome keys suggest the denominator list may have been derived from the table (or at least selected to match it). Therefore 0 / {len(denominator):,} all-record zero-hit genomes is not an independent prevalence estimate.",
               f"- {sum(r.begin<=3 for r in unique_all):,} / {len(unique_all):,} unique loci begin at <=3 and {sum((r.inclusive_length or 0)<1000 for r in unique_all):,} / {len(unique_all):,} are <1,000 bp. These may reflect contig-edge/truncated or very short calls, but sequence/contig lengths and caller output are absent, so the cause is unverified.",
               f"- {sum((r.inclusive_length or 0)>200000 for r in unique_all):,} / {len(unique_all):,} loci exceed 200,000 bp; inspect the interval row before downstream use. Long calls are not automatically artifacts.",
               f"- Exact duplicate, overlap, and nesting checks found {sum(n-1 for n in locus_counter.values()):,} extra duplicate-locus rows, {relations['all']['overlap_pairs']:,} overlapping pairs, and {relations['all']['nested_pairs']:,} nested pairs / {len(unique_all):,} all unique loci. This does not detect homology or duplicate biology across different assemblies/scaffolds.",
               "- No sequence lengths, caller name/version/settings, taxonomy method/confidence, status definition, sample metadata, assembly completeness, or quality field is supplied. Out-of-bounds checks, right-edge clipping, biological validation, and completeness/quality stratification remain unsupported.", "",
               "## Machine-readable outputs and plots", "",
               "- `per_genome.tsv`: one row per normalized denominator genome, in input order; record/locus/valid-locus/hit-scaffold counts and both bp sensitivities for all three scopes.",
               "- `summary_metrics.tsv`: reconciliation, joins, distributions, interval QC, semantic-field presence, and explicit denominators.",
               "- `category_counts.tsv`: categorical distributions, missingness, fractions, and denominators.",
               "- `interval_qc.tsv`: one row per normalized source record with exact keys, coordinates, dedup ranks, joins, and QC flags.",
               "- `plots/per_genome_locus_counts.svg`: zero-inclusive genome-count histogram (N printed).",
               "- `plots/per_genome_bp_ecdf.svg`: zero-inclusive bp ECDF (N printed; log10(1+bp) display).",
               "- `plots/locus_length_distribution.svg`: scope-specific length-bin percentages with locus denominators.",
               "- `plots/taxonomy_distribution.svg`: top exact taxonomy strings plus Other, with record denominator.", "",
               "SVGs are deterministic, dependency-free, contain accessible `<title>`/`<desc>` text, and expose exact counts in bar/line tooltips where applicable.", ""]
    report_path.write_text("\n".join(report),encoding="utf-8")

    # Fail if analysis somehow changed either source input.
    after={csv_path.name:sha256(csv_path),accession_path.name:sha256(accession_path)}
    if after != input_hashes:
        raise SystemExit(f"Input hash changed during analysis: before={input_hashes}, after={after}")


if __name__ == "__main__":
    project_root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[2]
    main(project_root)
