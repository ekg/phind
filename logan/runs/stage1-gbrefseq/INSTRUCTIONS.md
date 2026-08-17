# Logan Search dashboard submission instructions

Run: stage1-gbrefseq-t70
Generated: 2026-08-17T20:57:59.516+00:00

For each file in submissions/ (in this order, one at a time):

1. Open https://logan-search.org/dashboard
2. "Submit a query": upload (or paste) the single-sequence FASTA.
3. Group: GenBank_RefSeq   (recorded in the manifest row; do not change)
4. Threshold: 0.70   (recorded in the manifest row; do not change)
5. Email: leave empty (no notifications; results retained one month (docs.logan-search.org)).
6. Submit. The result page URL ends with the session id, e.g.
   https://logan-search.org/dashboard/kmviz-<uuid>  ->  session id
   kmviz-<uuid>. Copy it.
7. Record it WITHOUT any network request from this machine:
     python3 run_pilot.py record --run-dir . --bait <bait_id> \
         --session kmviz-<uuid>
8. Only after all submissions are recorded (or in bounded batches):
     python3 run_pilot.py fetch --run-dir .

Rate-limit policy: at most one dashboard submission per 60.0 s and
at most 40 HTTP requests per fetch invocation. Do not run
multiple fetch processes concurrently against the same run directory.
