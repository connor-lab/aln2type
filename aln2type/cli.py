import argparse
from .aln2type import go

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("oj", help="Output directory for typing JSON")
    parser.add_argument("ocsv", help="Output summary CSV file")
    parser.add_argument("refname", help="Name of reference sequence")
    parser.add_argument("msa", help="Path to MSA")
    parser.add_argument("typing_yaml", nargs='+', help="Path to Variant definition YAML files")
    args = parser.parse_args()

    return args

def main():
    go(parse_args())
