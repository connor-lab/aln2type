# aln2type

Use [variant_definitions](https://github.com/phe-genomics/variant_definitions) to type an MSA of SARS-CoV-2 sequences, handling insertions.

MSA must include MN908947.3 to correctly label variant positions.

## Help
```
usage: aln2type [-h] [--csv_N] [--no_gzip_json]
                json_outdir sample_csv_outdir summary_csv_outfile ref_name msa
                typing_yaml [typing_yaml ...]

positional arguments:
  json_outdir          Output directory for typing JSON
  sample_csv_outdir    Output directory for sample variant CSVs
  summary_csv_outfile  Output summary CSV file
  ref_name             Name of reference sequence in MSA
  msa                  Path to MSA
  typing_yaml          Path to Variant definition YAML files

optional arguments:
  -h, --help           show this help message and exit
  --csv_N              Include Ns in sample variant CSVs
  --no_gzip_json       Don't gzip typing JSON files
```

## Install

`git clone https://github.com/connor-lab/aln2type`

`cd aln2type`

`pip install .`

`git clone https://github.com/phe-genomics/variant_definitions`


## Run
`aln2type sample_json_out sample_csv_out typing_summary.csv MN908947.3 sars_cov_2.aln variant_definitions/variant_yaml/*.yml`
