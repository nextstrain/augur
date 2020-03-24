#!/bin/sh

case "$@" in
	*-k*) partial_test=1 ;;
esac

if [ "$partial_test" = 1 ]; then
	# skip test coverage when running a subset of the test suite
	coverage_arg='--no-cov'
else
	coverage_arg=''
fi

python3 -m pytest -c pytest.python3.ini $coverage_arg "$@"

exit $?
