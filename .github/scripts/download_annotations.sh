#!/usr/bin/env bash

set -eo pipefail

# Download annotations
curl --silent -L  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_045512.2&rettype=gb&retmode=txt" > NC_045512.2.gb