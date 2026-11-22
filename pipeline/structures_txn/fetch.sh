#!/usr/bin/env bash
# Fetch the mmCIF files this directory builds from. Re-run only if dl/ is missing.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p dl && cd dl
for id in 6GMH 5XJC 4PJO; do
  l=$(echo "$id" | tr 'A-Z' 'a-z')
  [ -f "$l.cif" ] && { echo "have $l.cif"; continue; }
  nice -n 19 curl -fsS --max-time 300 -o "$l.cif.gz" "https://files.rcsb.org/download/$id.cif.gz"
  gunzip -f "$l.cif.gz"
  echo "fetched $l.cif"
done
# Alternatives kept for reference (see txn_meta.json rejected_alternatives):
#   9EHZ human Pol II EC 2.6 A (short scaffold), 5FLM pig Pol II EC 3.4 A, 8XSO human 2.7 A
