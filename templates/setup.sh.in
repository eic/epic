#!/bin/sh

## Error if not the right name (this script is sourced, hence $1)
if [[ "$(basename ${BASH_SOURCE[0]})" != "setup.sh" ]]; then
        echo "Error: This script has ceased to exist at '$(realpath --no-symlinks ${BASH_SOURCE[0]})'."
        echo "       Please use the version at '$(realpath --no-symlinks $(dirname ${BASH_SOURCE[0]})/setup.sh)'."
        return 1
fi
