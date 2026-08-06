#!/usr/bin/env bash
# Full 26k host mash tree build.
set -e
MASH=/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/run-host-structure-1000/host-structure-1000-v1-run/tool-env/bin/mash
SRC=/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/26k/canonical_objects
BASE=$(cd "$(dirname "$0")" && pwd)
WORK=$BASE/data
IDS=$BASE/ids.txt
DIST=$BASE/host_mash.dist
MSH=$BASE/hosts.msh

echo "== Step 1: build chunk fastas + ids.txt ($(date))"
python3 $BASE/scripts/build_chunks.py --src $SRC --workdir $WORK --ids $IDS --chunk-size 5000

echo "== Step 2: sketch + pairwise triangle + merge ($(date))"
python3 $BASE/scripts/host_pipeline.py --workdir $WORK --out $DIST --ids $IDS \
  --msh $MSH --chunk-size 5000 --threads 32 --procs 8 --mash $MASH

echo "== Step 3: verify triangle ($(date))"
python3 $BASE/scripts/verify_triangle.py --ids $IDS --dist $DIST --msh $MSH \
  --src $SRC --mash $MASH --out $BASE/triangle_verify.json

echo "== Step 4: UPGMA tree ($(date))"
python3 $BASE/scripts/build_tree.py --ids $IDS --dist $DIST --out $BASE/host_tree

echo "== Step 5: validate tree ($(date))"
python3 $BASE/scripts/validate_tree.py --ids $IDS --dist $DIST \
  --nwk $BASE/host_tree.nwk --out $BASE/tree_verify.json

echo "== DONE ($(date)) =="
