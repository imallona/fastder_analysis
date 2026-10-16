#!/usr/bin/env python3
"""Wait until the one-minute load is below a limit; exit 1 on timeout."""

import argparse
import os
import sys
import time


def wait_until_quiet(below, timeout_s, read_load=lambda: os.getloadavg()[0],
                     sleep=time.sleep, interval_s=10):
    """True once the load is below the limit, False if timeout_s passes first."""
    waited = 0
    while read_load() >= below:
        if waited >= timeout_s:
            return False
        sleep(interval_s)
        waited += interval_s
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--below", type=float, required=True)
    parser.add_argument("--timeout", type=int, default=600, help="seconds")
    args = parser.parse_args()
    if wait_until_quiet(args.below, args.timeout):
        return 0
    print(f"wait_quiet: load {os.getloadavg()[0]:.2f} stayed at or above {args.below} "
          f"for {args.timeout} s; timed pass not started", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
