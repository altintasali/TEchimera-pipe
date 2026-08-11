#!/usr/bin/env bash
# install-env.sh -- download and install the pre-built conda environment for
# techimera-pipe, published as a GitHub Release asset alongside each
# version tag. The tarball is built for Linux x86_64 (see
# .github/workflows/release-env.yml); there is no need to install conda.
#
# Usage:
#   ./workflow/scripts/install-env.sh [-o PREFIX] [-r VERSION] [-f]
#
# Options:
#   -o PREFIX   install into PREFIX (default: $HOME/software/techimera-pipe-env)
#   -r VERSION  release tag to fetch, e.g. v0.1.0 (default: latest release)
#   -f          overwrite even if PREFIX doesn't look like a previous env
#               install, and skip the platform check
#   -h          show this help
#
# An existing environment at PREFIX is replaced automatically with a warning
# message; -f is only needed when PREFIX holds unrelated files.
#
# Afterwards, either activate the environment in your current shell:
#   source workflow/scripts/activate-env.sh [PREFIX]
# or just prepend its bin directory to your PATH (handy on SLURM compute nodes):
#   export PATH="$PREFIX/bin:$PATH"
set -euo pipefail

repo="altintasali/TEchimera-pipe"
stem="techimera-pipe"
default_prefix="$HOME/software/${stem}-env"

usage() {
    sed -n '2,23p' "$0"
    exit "${1:-0}"
}

prefix="$default_prefix"
version=""
force=0
while getopts ":o:r:fh" opt; do
    case $opt in
        o) prefix=$OPTARG ;;
        r) version=$OPTARG ;;
        f) force=1 ;;
        h) usage ;;
        :) echo "option -$OPTARG requires an argument" >&2; usage 1 ;;
        \?) echo "invalid option: -$OPTARG" >&2; usage 1 ;;
    esac
done

if [[ ! -e "$prefix" ]]; then
    mkdir -p "$prefix"
fi

# Platform check first: never touch a prefix on an OS/arch the tarball can't
# support -- wiping an existing env and then aborting would leave nothing.
if [[ $force -eq 0 ]]; then
    os=$(uname -s)
    arch=$(uname -m)
    if [[ "$os" != "Linux" || "$arch" != "x86_64" ]]; then
        echo "error: the pre-built environment is Linux x86_64 only (got $os/$arch);" >&2
        echo "       create the env locally instead with: conda env create -f workflow/environment.yaml" >&2
        exit 1
    fi
fi

# Overwriting an existing install. A prefix that looks like a previous
# environment (conda-pack ships bin/conda-unpack, conda ships bin/activate) is
# replaced automatically with a warning; a prefix holding unrelated files is
# left alone unless -f is given, so a stray -o can't destroy data by accident.
if [[ -n "$(ls -A "$prefix")" ]]; then
    if [[ $force -eq 0 && ! -f "$prefix/bin/conda-unpack" && ! -f "$prefix/bin/activate" ]]; then
        echo "error: $prefix is not empty and does not look like a previously-installed" >&2
        echo "       techimera-pipe environment; pick another prefix with -o," >&2
        echo "       or re-run with -f to overwrite it anyway." >&2
        exit 1
    fi
    echo "Warning: replacing the existing environment at $prefix." >&2
    rm -rf "$prefix"/* "$prefix"/.[!.]*
fi

# Resolve the release to fetch. This deliberately avoids the GitHub API (which
# rate-limits anonymous requests); it uses the releases Atom feed for the
# latest tag, and SHA256SUMS for the asset list. RNASEQ_INSTALL_ORIGIN
# overrides the host for tests and mirrors.
origin="${RNASEQ_INSTALL_ORIGIN:-https://github.com/$repo}"
if [[ -z "$version" ]]; then
    version=$(curl -fsSL --retry 3 "$origin/releases.atom" 2>/dev/null \
        | grep -oE 'releases/tag/[^"<]+' \
        | head -n1 \
        | sed -E 's#.*/##')
fi
if [[ -z "$version" ]]; then
    echo "error: could not determine a release to download" >&2
    exit 1
fi

releases_dl="$origin/releases/download/$version"
work_dir=$(mktemp -d)
trap 'rm -rf "$work_dir"' EXIT

# SHA256SUMS is shipped with every env release and lists the exact files to
# fetch (the single tarball, or the split .part.* chunks), so we read the file
# list from it instead of querying the API.
echo "Downloading $stem env assets from release $version ..."
if ! curl -fsSL --progress-bar --retry 3 -o "$work_dir/SHA256SUMS" "$releases_dl/SHA256SUMS"; then
    echo "error: release $version has no $stem env assets (no SHA256SUMS)" >&2
    exit 1
fi
(
    cd "$work_dir"
    while IFS= read -r fname; do
        [[ -n "$fname" ]] || continue
        echo "  $releases_dl/$fname"
        curl -fsSL --progress-bar --retry 3 -o "$fname" "$releases_dl/$fname"
    done < <(awk '{print $2}' SHA256SUMS)
)

# Verify checksums (or at least download SHA256SUMS to confirm they line up).
if command -v sha256sum >/dev/null 2>&1; then
    cksum="sha256sum"
elif command -v shasum >/dev/null 2>&1; then
    cksum="shasum -a 256"
else
    echo "error: neither sha256sum nor shasum found on PATH" >&2
    exit 1
fi
(cd "$work_dir" && $cksum -c SHA256SUMS)

# Recombine split parts, if the release shipped them.
tarball="$work_dir/${stem}-${version}-env.tar.gz"
if [[ ! -f "$tarball" ]]; then
    if compgen -G "$tarball.part.*" >/dev/null; then
        cat "$tarball".part.* > "$tarball"
    else
        echo "error: no env tarball (or split parts) found in the release" >&2
        exit 1
    fi
fi

echo "Extracting into $prefix ..."
tar -xzf "$tarball" -C "$prefix"

echo "Relocating hard-coded prefixes (conda-unpack) ..."
"$prefix/bin/conda-unpack"

echo "Sanity check ..."
"$prefix/bin/python" -c "import snakemake; print('snakemake', snakemake.__version__)"

echo
echo "Done. Activate it in your shell with:"
echo "  source workflow/scripts/activate-env.sh"
if [[ "$prefix" != "$default_prefix" ]]; then
    echo "  # or, since you installed with a custom -o prefix:"
    echo "  source workflow/scripts/activate-env.sh \"$prefix\""
fi
echo "or add its bin directory to your PATH (e.g. once in ~/.bashrc or on SLURM compute nodes):"
echo "  export PATH=\"$prefix/bin:\$PATH\""
echo "Then run the workflow as usual:"
echo "  snakemake --configfile config/test.yaml --cores 4   # smoke test"
echo
echo "For a shared cluster install, extract to shared storage once (e.g. -o /shared/software/${stem}-env)"
echo "so all nodes see it; no conda or container runtime is needed."
