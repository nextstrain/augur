#!/usr/bin/env python3
"""Run cram tests in parallel using a worker pool."""
import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# These tests were identified as being particularly slow <https://github.com/nextstrain/augur/issues/1994>
# We can re-check these over time as we speed up individual tests.
# We run these tests first to improve parallel efficiency.
SLOW_TESTS = [
    "tests/functional/merge/cram/merge-metadata.t",
    "tests/functional/tree/cram/iqtree-more-threads.t",
    "tests/functional/subsample/cram/proximal-subsampling.t",
    "tests/functional/measurements_export.t",
    "tests/functional/curate/cram/metadata-input.t",
    "tests/functional/export_v2/cram/metadata-columns.t",
    "tests/functional/curate/cram/titlecase.t",
    "tests/functional/tree/cram/iqtree-override-args.t",
    "tests/functional/subsample/cram/proximal-subsampling-errors.t",
    "tests/functional/merge/cram/merge-metadata-and-sequences.t",
]


def run_test(test_file, cram_args):
    start = time.monotonic()
    result = subprocess.run(
        ["cram", *cram_args, str(test_file)],
        capture_output=True,
    )
    elapsed = time.monotonic() - start
    return test_file, result.returncode, elapsed, result.stdout, result.stderr


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        usage="%(prog)s [OPTIONS] [TESTS] [-- CRAM_ARGS...]",
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=os.cpu_count(),
        help="number of parallel workers (default: all available (%(default)s))",
    )
    parser.add_argument(
        "tests", nargs="*", default=["tests/"],
        help="files or directories to find .t files in (default: tests/)",
    )

    argv = sys.argv[1:]
    if "--" in argv:
        split = argv.index("--")
        args = parser.parse_args(argv[:split])
        cram_args = argv[split + 1:]
    else:
        args = parser.parse_args(argv)
        cram_args = []

    all_tests = []
    for path in map(Path, args.tests):
        if path.is_file():
            all_tests.append(path)
        elif path.is_dir():
            all_tests.extend(path.rglob("*.t"))
        else:
            parser.error(f"not a file or directory: {path}")
    all_tests = sorted(set(all_tests))
    if not all_tests:
        parser.error(f"no .t files found in {' '.join(args.tests)}")

    slow_set = [Path(p) for p in SLOW_TESTS]
    slow_tests = [t for t in slow_set if t in all_tests]
    rest_tests = [t for t in all_tests if t not in slow_set]
    test_files = slow_tests + rest_tests

    cram_cmd = " ".join(["cram", *cram_args])
    print(f"Running {len(test_files)} tests with {args.jobs} workers")
    print(f"  ({len(slow_tests)} slow tests scheduled first)")
    print(f"cram invocation: {cram_cmd} <test>\n")

    results = []
    passed = failed = 0
    wall_start = time.monotonic()

    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(run_test, t, cram_args): t for t in test_files
        }
        for future in as_completed(futures):
            test_file, rc, elapsed, stdout, stderr = future.result()
            results.append((elapsed, rc, test_file))
            status = "PASS" if rc == 0 else "FAIL"
            if rc == 0:
                passed += 1
            else:
                failed += 1
            print(f"  {status}  {elapsed:6.1f}s  {test_file}")
            if rc != 0:
                if stdout:
                    print(stdout.decode(errors="replace"))
                if stderr:
                    print(stderr.decode(errors="replace"))

    wall_elapsed = time.monotonic() - wall_start
    total_cpu = sum(e for e, _, _ in results)

    print(f"\n{'='*60}")
    print(f"Passed: {passed}  Failed: {failed}  Total: {len(results)}")
    print(f"Wall time:  {wall_elapsed:.1f}s")
    print(f"Total CPU:  {total_cpu:.1f}s")
    print(f"Speedup:    {total_cpu / wall_elapsed:.1f}x")

    failures = [(e, rc, t) for e, rc, t in results if rc != 0]
    if failures:
        print(f"\n{'='*60}")
        print(f"Failing tests ({len(failures)}):\n")
        for elapsed, rc, test_file in sorted(failures, key=lambda x: x[2]):
            print(f"  {elapsed:6.1f}s  exit={rc}  {test_file}")

    print(f"\n{'='*60}")
    print("Slowest tests:\n")
    results.sort(reverse=True)
    for elapsed, rc, test_file in results[:20]:
        status = "PASS" if rc == 0 else "FAIL"
        print(f"  {elapsed:6.1f}s  {status}  {test_file}")

    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
