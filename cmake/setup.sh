#!/bin/sh

echo "Error: This script has ceased to exist at '$(realpath --no-symlinks ${BASH_SOURCE[0]})'."
echo "       Please use the version at '$(realpath --no-symlinks $(dirname ${BASH_SOURCE[0]})/bin/thisepic.sh)'."
return 1
