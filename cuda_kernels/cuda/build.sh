#!/usr/bin/env bash
set -euo pipefail

/usr/local/cuda/bin/nvcc \
  -arch=sm_87 \
  -O3 \
  --compiler-options '-fPIC' \
  -shared horb_accel.cu \
  -o libhorb_cuda.so

ls -lh libhorb_cuda.so
