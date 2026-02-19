# NOTE: This script only works when called from the augur repo root.
base="$(pwd)"
env="$base/.env"

source "$base/devel/common.sh"

conda-ish activate "$env"
