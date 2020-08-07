
function errorFound {
  echo -e "\nTest Script Failed at Line $1"
  exit 2
}
trap 'errorFound $LINENO' ERR
shopt -s nullglob # glob results are empty if no matches

echo -e "Running all tests\n-----------------\n\n"

cd $(dirname "$BASH_SOURCE")

if [ -d ./auspice ]; then
    rm -rf auspice
fi
mkdir auspice


for snakefile in ./*/Snakefile; do
    echo -e "\nRunning ${snakefile} (quietly)\n"
    pushd $(dirname "${snakefile}") >/dev/null
    snakemake --cores 1 --quiet clean 1>/dev/null
    snakemake --cores 1 --quiet 1>/dev/null
    cp auspice/*.json ../auspice
    popd >/dev/null
done

echo -e "\nAll tests passed. You can view the results by running"
echo -e "auspice view --datasetDir $(pwd)/auspice"

