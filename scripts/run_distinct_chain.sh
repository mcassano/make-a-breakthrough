#!/bin/sh
# Exhaustive distinct-distance searches (onion order): validation against
# known OEIS A271490 values, then the open case n = 18, k = 15.
set -e
cd "$(dirname "$0")/.."
BIN=./target/release/distinct
LOG=results/distinct
for nk in "15 13" "16 14" "18 15" "17 14"; do
    set -- $nk
    echo "=== n=$1 k=$2 started $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> $LOG/chain.log
    $BIN $1 $2 --onion 2> $LOG/n$1_k$2.progress | tee -a $LOG/chain.log
done
echo "=== chain finished $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> $LOG/chain.log
