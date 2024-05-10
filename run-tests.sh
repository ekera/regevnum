#!/bin/sh

export PYTHONPATH=$PYTHONPATH:$PWD
mkdir -p tests/tmp

sage tests/test-logarithm-finding.sage
sage tests/test-factoring.sage
sage tests/test-order-finding.sage
