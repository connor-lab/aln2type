#!/usr/bin/env bash

set -eo pipefail

# Get current definitions release
current_gh_release=$(curl --silent "https://api.github.com/repos/phe-genomics/variant_definitions/releases/latest" | grep -Po '"tag_name": "\K.*?(?=")')

curl -L --silent https://github.com/phe-genomics/variant_definitions/archive/${current_gh_release}.tar.gz | tar -xz
ln -sf $(pwd)/variant_definitions-${current_gh_release} $(pwd)/variant_definitions-current
echo -e "Current definitions release: ${current_gh_release}\n\n" > artifacts/artifact.log

