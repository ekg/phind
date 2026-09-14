We are restarting the PHIND prophage-analysis project after a WG control-plane incident. Read `/home/erikg/phind/PHIND_RECOVERY_HANDOFF.md` first, then inspect the referenced recovery manifest and task/agent/session map read-only.

Important constraints:
- Do not modify `/home/erikg/phind.recovery-20260801T164014Z` or `/mnt/nvme3n1/erikg/phind-recovery-20260801T164014Z`.
- Do not resume old Pi sessions or bulk-clean old worktrees.
- Do not treat reconstructed status as exact historical truth; verify source commits and artifacts.
- Do not push or rewrite Git history: local `main` is 20 commits ahead and contains oversized scientific artifacts that GitHub rejected.
- Use the new empty WG graph only for present work.

Current scientific objective: construct a MASH triangle/tree from the FASTA of all prophage elements—not complete E. coli genomes—and map/highlight ECOR elements so each accessible prophage can be inspected. Also audit the partition-stitching implementation at commit `2363ece` and distinguish reusable outputs from failed IMPG/pggb attempts.

First response/action:
1. Summarize the recovered scientific state and identify exact artifact paths/commits supporting it.
2. Propose a minimal active WG graph (roughly 3–6 goal-bearing tasks) centered on artifact audit, the all-prophage MASH tree/triangle, ECOR mapping, and stitching validation.
3. Stop for human approval before publishing or launching computational work.
