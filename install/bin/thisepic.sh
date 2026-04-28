#!/bin/sh

export DETECTOR=epic
export DETECTOR_PATH=/home/igorko/EIC/epic_github/epic/install/share/epic
export DETECTOR_CONFIG="${1:-epic}"
export DETECTOR_VERSION=26.04.0-15-g308efb360-dirty

## Export detector libraries
case "$(uname -s)" in
    Darwin*)
        if [ -n "$DYLD_LIBRARY_PATH" ]; then
            export DYLD_LIBRARY_PATH="/home/igorko/EIC/epic_github/epic/install/lib:$DYLD_LIBRARY_PATH"
        else
            export DYLD_LIBRARY_PATH="/home/igorko/EIC/epic_github/epic/install/lib"
        fi
        ;;
    *)
        if [ -n "$LD_LIBRARY_PATH" ]; then
            export LD_LIBRARY_PATH="/home/igorko/EIC/epic_github/epic/install/lib:$LD_LIBRARY_PATH"
        else
            export LD_LIBRARY_PATH="/home/igorko/EIC/epic_github/epic/install/lib"
        fi
        ;;
esac
