#!/bin/bash

set -euo pipefail

mkdir -p ../../results/mrbayes/test_run

cp ../../results/mrbayes/test_run_mrbayes.nex ../../results/mrbayes/test_run/

cd ../../results/mrbayes/test_run

echo "Starting MrBayes test run at $(date)"
time mb test_run_mrbayes.nex > test_run_mrbayes.log 2>&1
echo "Finished MrBayes test run at $(date)"