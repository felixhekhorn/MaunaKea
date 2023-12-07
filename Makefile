CXX = c++
CXXFLAGS = -std=c++11 -Wall -Wextra -O3
PINEAPPL_DEPS != pkg-config --cflags --libs pineappl_capi
LHAPDF_DEPS != pkg-config --cflags --libs lhapdf
#GSL_DEPS != gsl-config --cflags --libs
DVEGAS_DEPS != pkg-config --cflags --libs dvegas

MaunaKea: main.cpp PineAPPL.hpp config.h FO.hpp Integration.hpp Kernel.hpp MaunaKea.hpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) $(LHAPDF_DEPS) $(DVEGAS_DEPS) -o $@

PHONY: clean

clean:
	rm -f MaunaKea
