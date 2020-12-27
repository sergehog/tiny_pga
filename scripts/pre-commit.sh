#!/usr/bin/env bash

set -e

echo "RUN CLANG-FORMAT"
./scripts/pre-commit-format.sh

echo "CHECK UNIT-TESTS"
./scripts/run-tests.sh

# $? stores exit value of the last command
if [ $? -ne 0 ]; then
 echo "!!!!!TESTS MUST PASS BEFORE COMMIT!!!!!"
 exit 1
fi
