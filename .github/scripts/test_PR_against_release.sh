#!/bin/bash
set -eo pipefail

# unzip testdata
gzip -cd .github/data/testset.fa.gz > testset.fa

# Run PR code
aln2type --csv_N  \
--gb NC_045512.2.gb \
--no_gzip_json \
--no_call_deletion \
new_pr/json \
new_pr/csv \
new_pr/summary.csv \
MN908947.3 \
testset.fa \
variant_definitions-current/variant_yaml/*.yml

# run tests against previous previous_release to compare outputs 
git clone https://github.com/connor-lab/aln2type.git previous_release 
cd previous_release
python -m venv venv && source venv/bin/activate && python -m pip install .

aln2type --csv_N  \
--gb ../NC_045512.2.gb \
--no_gzip_json \
--no_call_deletion \
../main/json \
../main/csv \
../main/summary.csv \
MN908947.3 \
../testset.fa \
../variant_definitions-current/variant_yaml/*.yml

cd ..

if ! git diff --stat --no-index new_pr main > diffs.txt ; then
  echo "test failed: differences found between PR and previous release" >> artifacts/artifact.log
  echo see diffs.txt >> artifacts/artifact.log 
  cp diffs.txt artifacts/
  tar -cvf artifacts/outputs.tgz new_pr main
  exit 1
else
  tar -cvf artifacts/outputs.tgz new_pr main
  echo no differences found between PR and previous release >> artifacts/artifact.log
fi


