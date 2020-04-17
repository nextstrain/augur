
function errorFound {
  echo -e "\nTest Script Failed at Line $1"
  exit 2
}
trap 'errorFound $LINENO' ERR
shopt -s nullglob # glob results are empty if no matches

echo -e "\nRunning all tests\n-----------------\n\n"

cd $(dirname "$BASH_SOURCE")

# exclude_test_dirs=(export filter import mask refine traits translate tree validate various_export_settings zika runner.sh)
exclude_test_dirs=(traits runner.sh)

for test_dir in *; do
    if [[ ! " ${exclude_test_dirs[@]} " =~ " ${test_dir} " ]]; then
      echo -e "\nRunning tests in ${test_dir}...\n"
      cram -v $test_dir
    else
      echo -e "\nSkipping all tests in ${test_dir}.\n"
    fi
done

echo -e "\nAll tests passed."
