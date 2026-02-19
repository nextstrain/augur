conda-ish() {
    for condaish in micromamba mamba conda; do
        if type -a "$condaish" >/dev/null 2>&1; then
            "$condaish" "$@"
            return $?
        fi
    done
    echo "No conda-ish program (micromamba, mamba, conda) found" >&2
    return 1
}
