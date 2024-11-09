#!/bin/bash
# Default values
QUIET=0
VERBOSE=0
KEEP_TMPDIR=0
SHELL_PATH="/bin/bash"
SHELL_OPTS=""
INDENT=2
TEST_DIR="tests"
PRIORITY_PATTERNS="merge|tree"  # Slow tests that should start first
# Get default number of processors
NCPU=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
PARALLEL_JOBS=$NCPU  # Default to CPU count if not specified

# Function to show usage
usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] [TESTS...]
Wrapper script to run cram tests in parallel.
Options:
  -h, --help          show this help message and exit
  -q, --quiet         don't print diffs
  -v, --verbose       show filenames and test status
  -j, --jobs=N        number of tests to run in parallel (default: number of CPUs)
  --keep-tmpdir       keep temporary directories
  --shell=PATH        shell to use for running tests (default: /bin/bash)
  --shell-opts=OPTS   arguments to invoke shell with
  --indent=NUM        number of spaces to use for indentation (default: 2)
EOF
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            ;;
        -q|--quiet)
            QUIET=1
            shift
            ;;
        -v|--verbose)
            VERBOSE=1
            shift
            ;;
        -j|--jobs=*)
            if [[ "$1" == *"="* ]]; then
                PARALLEL_JOBS="${1#*=}"
            else
                shift
                PARALLEL_JOBS="$1"
            fi
            if ! [[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]]; then
                echo "Error: -j/--jobs requires a numeric argument"
                exit 1
            fi
            if [ "$PARALLEL_JOBS" -lt 1 ]; then
                echo "Error: Number of jobs must be at least 1"
                exit 1
            fi
            shift
            ;;
        -j*)  # Handle -j4 format
            PARALLEL_JOBS="${1#-j}"
            if ! [[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]]; then
                echo "Error: -j requires a numeric argument"
                exit 1
            fi
            if [ "$PARALLEL_JOBS" -lt 1 ]; then
                echo "Error: Number of jobs must be at least 1"
                exit 1
            fi
            shift
            ;;
        --keep-tmpdir)
            KEEP_TMPDIR=1
            shift
            ;;
        --shell=*)
            SHELL_PATH="${1#*=}"
            shift
            ;;
        --shell-opts=*)
            SHELL_OPTS="${1#*=}"
            shift
            ;;
        --indent=*)
            INDENT="${1#*=}"
            shift
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            ;;
        *)
            TEST_DIR="$1"
            shift
            ;;
    esac
done

# Validate test directory
if [ ! -d "$TEST_DIR" ]; then
    echo "Error: Directory $TEST_DIR does not exist"
    exit 1
fi

# Build cram command with options
CRAM_CMD="cram"
[ $QUIET -eq 1 ] && CRAM_CMD="$CRAM_CMD -q"
[ $VERBOSE -eq 1 ] && CRAM_CMD="$CRAM_CMD -v"
[ $KEEP_TMPDIR -eq 1 ] && CRAM_CMD="$CRAM_CMD --keep-tmpdir"
[ -n "$SHELL_PATH" ] && CRAM_CMD="$CRAM_CMD --shell=$SHELL_PATH"
[ -n "$SHELL_OPTS" ] && CRAM_CMD="$CRAM_CMD --shell-opts=$SHELL_OPTS"
[ -n "$INDENT" ] && CRAM_CMD="$CRAM_CMD --indent=$INDENT"

# Find all .t files and order with priority files first
TEST_FILES=$(find "$TEST_DIR" -name "*.t" | sort)
PRIORITY_FILES=$(echo "$TEST_FILES" | grep -E "$PRIORITY_PATTERNS" || true)
OTHER_FILES=$(echo "$TEST_FILES" | grep -vE "$PRIORITY_PATTERNS" || true)
ALL_FILES=$(echo -e "${PRIORITY_FILES}\n${OTHER_FILES}")

# Create a lock directory for synchronized output
LOCK_DIR="/tmp/cram_lock_$$"
mkdir "$LOCK_DIR"
trap 'rm -rf "$LOCK_DIR"' EXIT

# Create files to store temporary directories and test results
TMPDIR_FILE="$LOCK_DIR/tmpdirs"
RESULTS_FILE="$LOCK_DIR/results"
touch "$TMPDIR_FILE" "$RESULTS_FILE"

# Export variables for subshell
export CRAM_CMD LOCK_DIR TMPDIR_FILE RESULTS_FILE KEEP_TMPDIR

# Run tests in parallel with specified number of jobs
echo "$ALL_FILES" | xargs -P "$PARALLEL_JOBS" -I {} sh -c '
    file="$1"
    lock_dir="$LOCK_DIR"
    output=$($CRAM_CMD "$file" 2>&1)
    status=$?
    
    while ! mkdir "$lock_dir/lock" 2>/dev/null; do
        sleep 0.05
    done
    
    # Record test status
    echo "$status" >> "$RESULTS_FILE"
    
    if [ $status -eq 0 ]; then
        printf "."
    else
        echo "$output" | grep -v "^# Ran "
    fi
    
    # If --keep-tmpdir is enabled, capture the temporary directory
    if [ $KEEP_TMPDIR -eq 1 ]; then
        echo "$output" | grep "Kept temporary directory:" >> "$TMPDIR_FILE"
    fi
    
    rm -rf "$lock_dir/lock"
' - {}
echo

# Display captured temporary directories if --keep-tmpdir was used
if [ $KEEP_TMPDIR -eq 1 ]; then
    sort -u "$TMPDIR_FILE" | while read -r line; do
        echo "$line"
    done
fi

# Check if any tests failed, i.e. non-zero exit status
# Count the number of failed tests
# Count the number of passed tests
failed_tests=$(grep -c -v "^0$" "$RESULTS_FILE")
passed_tests=$(grep -c "^0$" "$RESULTS_FILE")

# Print summary of test results
echo "Passed $passed_tests, failed $failed_tests tests"
if [ "$failed_tests" -gt 0 ]; then
    exit 1
fi

exit 0
