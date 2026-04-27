#!/bin/bash

set -euo pipefail

mkdir -p ../../results/mrbayes/full_run

cp ../../results/mrbayes/full_run_mrbayes.nex ../../results/mrbayes/full_run/

cd ../../results/mrbayes/full_run

echo "Starting MrBayes full run at $(date)"
time mb full_run_mrbayes.nex > _full_run_mrbayes.log 2>&1
echo "Finished MrBayes full run at $(date)"