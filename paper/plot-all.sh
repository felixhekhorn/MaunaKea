#!/usr/bin/bash
# Iterate all groups. Pass `pto`, `alphas`, or `mass` as argument.

# charm
python plot.py 1.51 3 --$1 --short-range
python plot.py 1.4 3 --$1 --short-range
python plot.py 1.3 3 --$1 --pdf_set CT18NNLO_rescaled_NF3 --short-range

# bottom
python plot.py 4.92 4 --$1 --short-range
python plot.py 4.75 4 --$1 --short-range
python plot.py 4.75 4 --$1 --pdf_set CT18NNLO_rescaled_NF4 --short-range
