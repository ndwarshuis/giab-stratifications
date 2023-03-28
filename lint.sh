#! /usr/bin/env bash

# lint all scripts

# usage: lint.sh [ENV_BASE_PREFIX]
#
# where ENV_BASE_PREFIX is an optional location to store each environment
# (meant to be used by the CI/CD pipeline) and should correspond to what was
# used to create the environments using `setup_dev.sh`

anyfail=0 # global var (ew...) to track if we have any lint failures

_run_test () {
    "$@"
    local status=$?
    if (( status != 0 )) || (( anyfail == 1 )); then
        anyfail=1
    fi
}

_mypy () {
    # run mypy without the cache to ensure a clean (albeit slow) lint
    mypy --no-incremental --cache-dir=/dev/null "$1"
}

_conda_activate () {
    # use either a local env or a user-installed env
    base="$1"
    if [ -z "$prefix" ]; then
        conda activate giab-strats-"$base"
    else
        conda activate "./$prefix-$base"
    fi
}

eval "$(${CONDA_EXE} shell.bash hook 2> /dev/null)"

python_root="workflow/scripts/python"

if [ -n "$1" ]; then
    prefix="$1"
fi

# test bedtools scripts
_conda_activate bedtools

echo "Testing common module"
echo ""

# test the common dir since flake8 doesn't follow imports
_run_test flake8 "$python_root/common"
_run_test mypy "$python_root/common"

for mod in "$python_root/bedtools"/*; do
    if [[ "$mod" != *.yml ]]; then
        echo "Testing scripts for bedtools module: $mod"
        echo ""
    
        _run_test flake8 "$mod"
        _run_test _mypy "$mod"
    
        echo ""
    fi
done


if (( anyfail == 1 )); then
   exit 1
fi
