#!/usr/bin/env bash
# Scheduled Stage-0 smoke re-check for the bounded NTM Logan pilot.
#
# Preregistered rule (PILOT_PLAN.md section 4, service degradation): pause
# >= 24 h after the 2026-08-17 NO-GO, then ONE bounded re-check of the
# already-submitted session (no new submission, no threshold change).
# Window opens 2026-08-18T21:00Z (>= 24 h after the 20:55:30Z submission).
set -u

TARGET_EPOCH=$(date -u -d "2026-08-18T21:00:30Z" +%s)
NOW=$(date -u +%s)
SLEEP=$((TARGET_EPOCH - NOW))

if [ "$SLEEP" -gt 0 ]; then
  echo "[$(date -u +%FT%TZ)] preregistered pause active; sleeping ${SLEEP}s until 2026-08-18T21:00:30Z"
  sleep "$SLEEP"
fi

echo "[$(date -u +%FT%TZ)] recheck window open — running bounded Stage-0 fetch (poll 6, min-interval 60)"
cd "$(dirname "$0")"
python3 run_pilot.py fetch --run-dir runs/stage0-smoke --poll 6 --min-interval 60
RC=$?
echo "[$(date -u +%FT%TZ)] fetch exit=${RC}"
echo "[$(date -u +%FT%TZ)] manifest status after fetch:"
python3 run_pilot.py status --run-dir runs/stage0-smoke
echo "[$(date -u +%FT%TZ)] DONE rc=${RC}"
