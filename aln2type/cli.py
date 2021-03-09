import argparse
from .aln2type import go

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("json_outdir", help="Output directory for typing JSON")
    parser.add_argument("sample_csv_outdir", help="Output directory for sample variant CSVs")
    parser.add_argument("summary_csv_outfile", help="Output summary CSV file")
    parser.add_argument("ref_name", help="Name of reference sequence in MSA")
    parser.add_argument("msa", help="Path to MSA")
    parser.add_argument("typing_yaml", nargs='+', help="Path to Variant definition YAML files")
    args = parser.parse_args()

    return args

def main():
    go(parse_args())
