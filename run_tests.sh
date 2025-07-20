#!/bin/sh

# Default to running all tests without coverage.
partial_test=0
run_coverage=0

# Check if user explicitly requests coverage
case "$@" in
    *--cov*) run_coverage=1 ;;
esac

# If user requests a subset of tests for pytest, note this preference to avoid
# running other tests.
case "$@" in
    *-k*) partial_test=1 ;;
esac

if [ "$run_coverage" = 1 ]; then
    # Remove --cov from args since pytest will handle coverage via --cov=augur
    echo "Running tests with coverage enabled, remove --cov to disable"
    filtered_args=$(echo "$*" | sed 's/--cov//g')
    coverage_arg='--cov=augur'

    # Set env variable COVERAGE_CORE to sysmon for coverage
    # This reduces overhead of coverage collection
    export COVERAGE_CORE=sysmon
else
    # Default to no coverage
    echo "Running tests without coverage, use --cov to enable"
    filtered_args="$*"
    coverage_arg='--no-cov'
fi

echo "Running unit tests and doctests with pytest"
python3 -m pytest $coverage_arg $filtered_args

# Only run functional tests if we are not running a subset of tests for pytest.
if [ "$partial_test" = 0 ]; then
    echo "Running functional tests with cram"
    if [ "$run_coverage" = 1 ]; then
        # Use coverage for cram tests when coverage is requested
        AUGUR="coverage run -a bin/augur" cram tests/
    else
        # Use regular augur binary when coverage is not requested
        cram tests/
    fi
else
    echo "Skipping functional tests when running a subset of unit tests"
fi

exit $?