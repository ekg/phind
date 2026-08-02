# Main history oversized-blob verification

Task: `verify-main-history` — verify `main` history stays free of oversized
blobs (>100 MB, GitHub GH001 hard limit) after the all-prophage-mash branch
lands, and fix if it does not.

Date: 2026-08-02 (UTC)

## Result: PASS — no action required

The all-prophage-mash work landed in `main` (commit `9c51cd4` "feat: all-prophage
MASH triangle + UPGMA tree", plus rescue checkpoint `102d904`) **without**
re-introducing the >100 MB blobs that previously blocked remote sync. Pushes to
`origin` succeed with no GH001 rejection.

## Checks performed

### 1. No blob > 100,000,000 bytes reachable from local `main`

```bash
git rev-list --objects main \
  | cut -d' ' -f1 \
  | git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize)' \
  | awk '$1=="blob" && $3 > 100000000 {print}'
```

Result: no output — **0 blobs > 100 MB**. Largest blob reachable from `main` is
37,928,119 bytes (`96d1e8c0...`), well under the limit.

### 2. Strip fix — not needed

No blob > 100 MB was found reachable from `main`, so no
`git-filter-repo --invert-paths` rewrite was required.

### 3. Push from a fresh worktree branch succeeds (no GH001)

This worktree branch (`wg/agent-11/verify-main-history`, created on top of the
landed all-prophage-mash `main`) was pushed to `origin` successfully — no GH001
rejection. The report you are reading is that branch's commit.

```bash
git push origin wg/agent-11/verify-main-history   # OK, no GH001
```

### 4. `origin/main` reflects local `main`

Local `main` was `102d904` while `origin/main` was `848ffcb` (2 commits behind,
fast-forward). Pushed:

```bash
git push origin main    # 848ffcb..102d904, clean fast-forward, no GH001
```

Verified after push:

```bash
git rev-parse origin/main main   # identical: 102d904bc191d15b433db156f751422e3236fa84
```

## Extra context

- `origin/wg/agent-7/all-prophage-mash` is `9c51cd4` (already rebased onto the
  rewritten main `848ffcb`) and contains **0 blobs > 100 MB**.
- The three stripped blobs from the original report
  (`prophage_homology_survey/full_prophages.fa` 3.26 GB,
  `pggb_analysis/cluster_1/cluster_1.alignments.paf` 678 MB,
  `cluster_1.mappings.paf` 76 MB) remain reachable **only** from local,
  non-remote refs: `refs/backup/agent-10/pre-fix-main`,
  `refs/heads/recovery/phind-control-plane-clean-20260801T164014Z`, and
  `refs/wg/candidates/*`, `refs/wg/finish-results/*`, `refs/wg/rescues/*`.
  None of these are pushed to `origin` and none are ancestors of `main`, so they
  cannot trigger GH001 on remote pushes. If the recovery/backup refs are ever
  deliberately pushed, they must be excluded first.

## Conclusion

`main` history is free of oversized blobs after all-prophage-mash landed. Remote
pushes (both `main` and a fresh worktree branch) succeed with no GH001
rejection. No fix was needed.
