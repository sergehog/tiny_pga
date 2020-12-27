#!/usr/bin/env bash

echo "RUN CLANG-FORMAT"
./scripts/run-format.sh

echo "CHECK UNIT-TESTS"
./scripts/run-tests.sh

# $? stores exit value of the last command
if [ $? -ne 0 ]; then
 echo "!!!!!TESTS MUST PASS BEFORE COMMIT!!!!!"
 exit 1
fi
