#!/usr/bin/bash
# Iterate all PDF groups and quarks.
# 1. argument: `pto`, `lumi`, `alphas`, or `mass`
# 2. argument: `+` (normal mode) or `-` (data mode)

# charm
python plot.py 1.51 ${2}3 --$1 --short-range
python plot.py 1.40 ${2}3 --$1 --short-range
python plot.py 1.30 ${2}3 --$1 --pdf_set CT18NNLO_rescaled_NF3 --short-range

# bottom
python plot.py 4.92 ${2}4 --$1 --short-range
python plot.py 4.75 ${2}4 --$1 --short-range
python plot.py 4.75 ${2}4 --$1 --pdf_set CT18NNLO_rescaled_NF4 --short-range
