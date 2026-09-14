# PHIND recovery manifest

Recovery timestamp: 20260801T164014Z

## Authoritative locations

- Fresh source checkout and fresh WG graph: `/home/erikg/phind`
- Preserved damaged checkout and all surviving worktrees: `/home/erikg/phind.recovery-20260801T164014Z`
- External recovery evidence: `/mnt/nvme3n1/erikg/phind-recovery-20260801T164014Z`
- Complete Git bundle: `/mnt/nvme3n1/erikg/phind-recovery-20260801T164014Z/git/phind-all-refs.bundle`
- Bundle SHA-256: `234c74cff8169034bfa5f977a4c211445672ee7de10e474b48664bb7066729fd`
- Task/agent/session map: `/mnt/nvme3n1/erikg/phind-recovery-20260801T164014Z/manifests/task-agent-session-map.tsv`
- Final chat captures: `tmux/chat-1-final.txt` and `tmux/chat-3-final.txt` under the external recovery directory
- Cached TUI captures: `tmux/tui-pid1210893.txt` and `prior-evidence/` under the external recovery directory

## Source state

- Fresh checkout `main`: `a5070e716fd4e06d5932c1aef8e77c829c79b1f4`
- Parent rescue commit: `2f991af1c4803b8debb2a63b27d17861a619946f`
- The only tree change in the recovery commit is deletion of the tracked `.wg` symlink.
- Fresh `main` is 20 commits ahead of GitHub `origin/main`.
- GitHub refused the recovery-branch push because local history contains files over 100 MB, including a 3.1 GB FASTA. Do not rewrite, force-push, or discard this history. Decide an LFS/artifact publication strategy separately.
- The old checkout remains 146 GB and was moved intact. Its 59 extant linked worktrees were repaired to point to the archived Git common directory.

## Recovered graph/scientific state

The old canonical `graph.jsonl` and configuration were destroyed in the control-plane incident and were not reconstructed as fictional history. Cached TUI evidence showed:

- `acquire-remaining-25k`: done
- `extract-full-132k`: done
- `full-prophage-homology`: done
- `full-132k-prophage`: done
- `focused-heatmap-connected`: done
- `pggb-per-cluster`: done
- `run-pggb-ancestral`: failed
- `research-impg-partition`: done
- `run-impg-partition`: failed
- `benchmark-find-fastest`: failed
- `run-all-wave`: failed
- `compare-ancestral-vs`: open
- `build-partition-stitching`: in progress at incident time; feature commit `2363ece2ae154779d1626b0a0138b0c978f371a0` is present in source history.

Latest recovered user intent from chat 3:

> Build a MASH triangle/tree from a FASTA containing all prophage elements, rather than full E. coli genomes, then map/highlight the ECOR elements so each accessible element can be inspected.

Parameter correction recorded in chat 3:

- pggb `p=90, l=2000` was considered inappropriate for divergent prophages.
- The intended phage-oriented retry used approximately `p=75, s=250, l=500, k=11, ani-diff=80`.
- Do not accept prior task terminal labels without verifying artifacts: historical evaluation incorrectly marked low-scoring/refused work as done.

## New WG runtime

The graph at `/home/erikg/phind/.wg` is deliberately fresh and empty.

- WG source/install: `0ef7802da66a215fab72d3f9dfb06bfb2e8dad02`
- Pi plugin: compat `0.2.0`, embedded cache install
- Exact route: `pi:openrouter:deepseek/deepseek-v4-flash-0731`
- Worker/chat/default/standard/premium reasoning: high
- Assigner/evaluator/FLIP/reviewer/fast reasoning: low
- Maximum agents: 8
- Auto-assign, auto-evaluate, and FLIP: enabled
- Automatic archive retention: disabled (`archive_retention_days = 0`)
- Direct low- and high-reasoning Pi canaries both returned the exact expected response.

## Recovery rules

1. Treat `/home/erikg/phind.recovery-20260801T164014Z` and `/mnt/nvme3n1/erikg/phind-recovery-20260801T164014Z` as read-only evidence.
2. Do not automatically resume old Pi sessions; inspect their JSONL copies read-only.
3. Do not bulk-clean or delete archived worktrees.
4. Do not claim old task terminal states as canonical without artifact verification.
5. Keep `.wg` outside Git and candidates.
6. Build only a small active graph representing current goals.
