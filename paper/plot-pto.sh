#!/usr/bin/bash

# charm
python plot.py 1.51 3 --pto --short-range
python plot.py 1.4 3 --pto --short-range
python plot.py 1.3 3 --pto --pdf_set CT18NNLO_rescaled_NF3 --short-range

# bottom
python plot.py 4.92 4 --pto --short-range
python plot.py 4.75 4 --pto --short-range
python plot.py 4.75 4 --pto --pdf_set CT18NNLO_rescaled_NF4 --short-range
