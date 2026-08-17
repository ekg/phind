# Logan Search dashboard submission instructions

Run: smoke-ecoli-k12-150bp
Generated: 2026-08-17T16:53:29.569+00:00

For each file in submissions/ (in this order, one at a time):

1. Open https://logan-search.org/dashboard
2. "Submit a query": upload (or paste) the single-sequence FASTA.
3. Group: GenBank_RefSeq   (recorded in the manifest row; do not change)
4. Threshold: 0.50   (recorded in the manifest row; do not change)
5. Email: leave empty (no notifications; results retained one month (docs.logan-search.org)).
6. Submit. The result page URL ends with the session id, e.g.
   https://logan-search.org/dashboard/kmviz-<uuid>  ->  session id
   kmviz-<uuid>. Copy it.
7. Record it WITHOUT any network request from this machine:
     python3 run_pilot.py record --run-dir . --bait <bait_id> \
         --session kmviz-<uuid>
8. Only after all submissions are recorded (or in bounded batches):
     python3 run_pilot.py fetch --run-dir .

Rate-limit policy: at most one dashboard submission per 30.0 s and
at most 6 HTTP requests per fetch invocation. Do not run
multiple fetch processes concurrently against the same run directory.
