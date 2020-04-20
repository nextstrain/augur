
function errorFound {
  echo -e "\nTest Script Failed at Line $1"
  exit 2
}
trap 'errorFound $LINENO' ERR
shopt -s nullglob # glob results are empty if no matches

echo -e "\nRunning all tests\n-----------------\n\n"

cd $(dirname "$BASH_SOURCE")

# exclude_test_dirs=(refine tree)
exclude_test_dirs=(export filter traits runner.sh)

for test_dir in *; do
    if [[ ! " ${exclude_test_dirs[@]} " =~ " ${test_dir} " ]]; then
      echo -e "\nRunning tests in ${test_dir}...\n"
      cram $test_dir --verbose --shell=/bin/bash
    else
      echo -e "\nSkipping all tests in ${test_dir}.\n"
    fi
done

echo -e "\nAll tests passed."
