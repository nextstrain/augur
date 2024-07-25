#!/bin/sh

# Default to running all tests.
partial_test=0

# If user requests a subset of tests for pytest, note this preference to avoid
# running other tests.
case "$@" in
    *-k*) partial_test=1 ;;
esac

if [ "$partial_test" = 1 ]; then
    # skip test coverage when running a subset of the test suite
    coverage_arg='--no-cov'
else
    coverage_arg=''
fi

echo "Running unit tests and doctests with pytest"
python3 -m pytest $coverage_arg "$@"

# Only run functional tests if we are not running a subset of tests for pytest.
if [ "$partial_test" = 0 ]; then
    echo "Running functional tests with cram"
    cram tests/
else
    echo "Skipping functional tests when running a subset of unit tests"
fi

exit $?
