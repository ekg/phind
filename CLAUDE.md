<!-- worksgood-managed-guide:v1:start -->
# WorksGood (`wg`) project guide

This file is the **layer-2** project guide for agents working in this
WorksGood task graph. It is NOT the universal chat-agent / worker-agent
contract — that is bundled inside the `wg` binary and emitted by:

```
wg agent-guide
```

Run `wg agent-guide` at session start (or read its output from a previous
session) to get the universal role contract: chat agent vs dispatcher vs worker
distinction, `## Validation` requirement, smoke-gate, cycle handling, git
hygiene, worktree isolation, "no built-in Task tool" rules, etc.

This file only covers things specific to this project. Add project-specific
build commands, test commands, architecture notes, and service recipes here.

**At the start of each session, run `wg quickstart` in your terminal to orient yourself.**
Use `wg service start` to dispatch work — do not manually claim tasks.

This guide is written to both `CLAUDE.md` and `AGENTS.md` and kept in
lock-step. The two files exist because Claude Code and Codex CLI look for
different filenames, but they should never drift in content. Any divergence is
a bug. Update both together.
<!-- worksgood-managed-guide:v1:end -->

## Long-running compute — MANDATORY process-tool pattern (workers AND chat)

**Why.** Provider stream timeouts have killed worker attempts on zai, lunaroute,
AND openrouter — always while an agent held a long foreground tool call or
busy-polled a background job (2026-09-12/13: ntm-v3-per attempts 1–3,
`agent-exit-nonzero`, "Provider failure detected after streaming: timeout";
historically ntm-v2-host on openrouter). The 15-hour ntm-v3-host run survived
only because it followed the detached pattern below.

**Rule.** For ANY command expected to run longer than ~5 minutes (aligners,
MASH triangles, batch drivers, Pharokka, downloads, index builds):

1. Start it with the **`process` tool** (not bash `&`, not `nohup`), with
   `notify.logMatches` on readiness/error patterns and `onFailure: turn`.
   Give it a specific `name` and check the process list first.
2. **End your turn immediately** after starting it. The notification wakes you
   on exit or on a matched log line — do not sleep, do not poll, do not hold
   the turn.
3. When woken: inspect output with SHORT commands (`tail`, counts), decide,
   and either end the turn again or take the next short step.
4. Keep every foreground bash call under ~5 minutes. If a driver loops for
   hours, chunk it (per-clade, per-batch) so each poll is short, or let the
   `process` tool hold the long loop and wake you per milestone.

Resume-first: check for existing partial outputs (NVMe `ntm/v3/clades/`,
driver logs, PAF counts) before launching anything — drivers must skip
completed work units.

The old `nohup ... &` + poll loop pattern still appears in historical recovery
notes; treat those as pre-process-tool workarounds. The `process` tool is the
-supported way to get wakeups.