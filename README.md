# aln2type

Use [variant_definitions](https://github.com/phe-genomics/variant_definitions) to type an MSA of SARS-CoV-2 sequences.

MSA must include MN908947.3.

Recreate the test data with:

`git clone https://github.com/connor-lab/aln2type`

`cd aln2type`

`pip install .`

`git clone https://github.com/phe-genomics/variant_definitions`

`aln2type testdata/outputs/json_out testdata/outputs/output.csv MN908947.3 testdata/alignment/test.aln variant_definitions/variant_yaml/*.yml`
