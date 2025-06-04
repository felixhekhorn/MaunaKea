#!/usr/bin/bash

NPROCESS=6

# charm
python run.py 1.20 -3 0 -c -m --process $NPROCESS
python run.py 1.25 -3 0 -c -m --process $NPROCESS
python run.py 1.30 -3 0 -c -m --process $NPROCESS
python run.py 1.35 -3 0 -c -m --process $NPROCESS
python run.py 1.40 -3 0 -c -m --process $NPROCESS
python run.py 1.45 -3 0 -c -m --process $NPROCESS
python run.py 1.50 -3 0 -c -m --process $NPROCESS
python run.py 1.51 -3 0 -c -m --process $NPROCESS
python run.py 1.55 -3 0 -c -m --process $NPROCESS

# bottom
python run.py 4.00 -4 0 -c -m --process $NPROCESS
python run.py 4.25 -4 0 -c -m --process $NPROCESS
python run.py 4.50 -4 0 -c -m --process $NPROCESS
python run.py 4.75 -4 0 -c -m --process $NPROCESS
python run.py 4.92 -4 0 -c -m --process $NPROCESS
python run.py 5.00 -4 0 -c -m --process $NPROCESS
python run.py 5.25 -4 0 -c -m --process $NPROCESS
python run.py 5.50 -4 0 -c -m --process $NPROCESS
