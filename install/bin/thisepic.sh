#!/bin/sh

export DETECTOR=epic
export DETECTOR_PATH=/w/eic-scshelf2104/users/anatochii/epic/install/share/epic
export DETECTOR_CONFIG=${1:-epic}
export DETECTOR_VERSION=26.02.0-11-g048c5469b

## Warn is not the right name (this script is sourced, hence $1)
if [[ "$(basename ${BASH_SOURCE[0]})" != "thisepic.sh" ]]; then
        echo "Warning: This script will cease to exist at '$(realpath --no-symlinks ${BASH_SOURCE[0]})'."
        echo "         Please use the version at '$(realpath --no-symlinks $(dirname ${BASH_SOURCE[0]})/bin/thisepic.sh)'."
fi

## Export detector libraries
if [[ "$(uname -s)" = "Darwin" ]] || [[ "$OSTYPE" == "darwin"* ]]; then
        export DYLD_LIBRARY_PATH="/w/eic-scshelf2104/users/anatochii/epic/install/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
else
        export LD_LIBRARY_PATH="/w/eic-scshelf2104/users/anatochii/epic/install/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi
