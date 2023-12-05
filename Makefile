CXX = c++
CXXFLAGS = -std=c++11 -Wall -Wextra -O3
PINEAPPL_DEPS != pkg-config --cflags --libs pineappl_capi
LHAPDF_DEPS != pkg-config --cflags --libs lhapdf
#GSL_DEPS != gsl-config --cflags --libs
DVEGAS_DEPS != pkg-config --cflags --libs dvegas
# config.h Integration.hpp $(DVEGAS_DEPS) 
MaunaKea: main.cpp PineAPPL.hpp 
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) $(LHAPDF_DEPS) -o $@

PHONY: clean

clean:
	rm -f MaunaKea
