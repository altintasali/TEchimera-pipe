#!/usr/bin/env bash
# activate-env.sh -- activate the pre-built techimera-pipe conda
# environment. SOURCE this file (it sets environment variables in your current
# shell and must run in it):
#
#   source workflow/scripts/activate-env.sh [PREFIX]
#
# PREFIX defaults to $HOME/software/techimera-pipe-env, matching the
# default install location used by install-env.sh. Pass a different PREFIX when
# you installed elsewhere with `./workflow/scripts/install-env.sh -o PREFIX`.
#
# On activation it prints a one-line confirmation (only when interactive) and
# verifies the environment looks complete; a missing or broken install is
# reported with a reinstall hint instead of half-activating.
_is_exec=0
if [[ -n "${BASH_SOURCE[0]:-}" ]]; then
    [[ "${BASH_SOURCE[0]}" == "$0" ]] && _is_exec=1    # bash: $0 == script path when run directly
elif [[ "$0" == *"activate-env.sh" ]]; then
    _is_exec=1                                          # zsh: no BASH_SOURCE; $0 is the script when run directly
fi
if [[ $_is_exec -eq 1 ]]; then
    echo "error: source this file instead of running it:" >&2
    echo "  source workflow/scripts/activate-env.sh [PREFIX]" >&2
    exit 1
fi

prefix="${1:-$HOME/software/techimera-pipe-env}"
if [[ ! -f "$prefix/bin/activate" ]]; then
    echo "error: no environment found at $prefix" >&2
    echo "       install it first with: ./workflow/scripts/install-env.sh [-o $prefix]" >&2
    return 1
fi
if [[ ! -x "$prefix/bin/python" || ! -x "$prefix/bin/snakemake" ]]; then
    echo "error: environment at $prefix looks incomplete (missing bin/python or bin/snakemake)." >&2
    echo "       Reinstall it with: ./workflow/scripts/install-env.sh -o $prefix" >&2
    return 1
fi
source "$prefix/bin/activate"

if [[ -t 2 ]]; then
    echo "Activated techimera-pipe environment at $prefix" >&2
    snakemake_version=$("$prefix/bin/python" -c "import snakemake; print(snakemake.__version__)" 2>/dev/null || true)
    [[ -n "$snakemake_version" ]] && echo "  snakemake $snakemake_version" >&2
fi
unset _is_exec snakemake_version 2>/dev/null || true
