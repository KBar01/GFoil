#!/bin/bash
# Usage:
#   bash tests/run_regression.sh --create-golden
#   bash tests/run_regression.sh --test
#   bash tests/run_regression.sh --build --test
set -e
python3 tests/regression_test.py "$@"
