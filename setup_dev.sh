#! /bin/bash

# Install all development envs for python code. Assumes mamba is installed

# usage: setup_dev.sh [ENV_BASE_PREFIX]
#
# where ENV_BASE_PREFIX is an optional location to store each environment
# (meant to be used by the CI/CD pipeline)

_update () {
    env_file="$1"
    shift
    mamba env update -f "$env_file" "$@" --prune
}

for i in workflow/scripts/{python,rmarkdown}/*; do
    base="$(basename "$i")"
    env_file="$i/env.yml"
    if [ "$base" != "common" ]; then
        if [ -z "$1" ]; then
            _update "$env_file" -n "giab-strats-$base"
        else
            local_path="$1-$base"
            if [ ! -e "$local_path" ]; then
                _update "$env_file" -p "$local_path"
            fi
        fi
        
    fi
done

