#!/usr/bin/env bash

set -o pipefail # get errors even when piped
set -o nounset # prevent using undeclared variables
set -o errexit # exit on command fail; allow failure: || true

cc -Weverything -O3 -framework Accelerate profile.c consensus_contour.c -o profile
./profile
