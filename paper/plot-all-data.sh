#!/usr/bin/bash

# charm
python plot.py 0 -3 --pdf --mass
python plot.py 1.51 -3 --pto

# bottom
python plot.py 0 -4 --pdf --mass
python plot.py 4.92 -4 --pto
