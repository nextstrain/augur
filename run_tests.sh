#!/bin/sh

# Default to running all tests without coverage.
run_coverage=0

# Check if user explicitly requests coverage
case "$@" in
    *--cov*) run_coverage=1 ;;
esac

if [ "$run_coverage" = 1 ]; then
    # Remove --cov from args since pytest will handle coverage via --cov=augur
    echo "Running tests with coverage enabled, remove --cov to disable"
    filtered_args=$(echo "$*" | sed 's/--cov//g')
    coverage_arg='--cov=augur'

    # Set env variable COVERAGE_CORE to sysmon for coverage
    # This reduces overhead of coverage collection
    # But only if Python version >3.14
    if python3 -c 'import sys; print(sys.version_info >= (3, 14))' | grep -q True; then
        echo "Using sysmon for coverage to reduce overhead"
        export COVERAGE_CORE=sysmon
    else
        echo "Using default coverage collection method"
    fi
else
    # Default to no coverage
    echo "Running tests without coverage, use --cov to enable"
    filtered_args="$*"
    coverage_arg='--no-cov'
fi

echo "Running unit tests, doctests,and functional tests with pytest"
python3 -m pytest $coverage_arg $filtered_args

exit $?
