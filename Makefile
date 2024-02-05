CXX = c++
CXXFLAGS = -std=c++11 -Wall -Wextra -O3
PINEAPPL_DEPS != pkg-config --cflags --libs pineappl_capi
LHAPDF_DEPS != pkg-config --cflags --libs lhapdf
GSL_DEPS != gsl-config --cflags --libs
DVEGAS_DEPS != pkg-config --cflags --libs dvegas
PY_DEPS != python3.10-config --includes --ldflags
PYBIND_DEPS != python3 -m pybind11 --includes

MaunaKea: main.cpp PineAPPL.hpp config.h FO.hpp Integration.hpp Kernel.hpp MaunaKea.hpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) $(LHAPDF_DEPS) $(GSL_DEPS) $(DVEGAS_DEPS) -o $@

dE: 1612.05582/dE.cpp PineAPPL.hpp config.h FO.hpp Integration.hpp Kernel.hpp MaunaKea.hpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) $(LHAPDF_DEPS) $(GSL_DEPS) $(DVEGAS_DEPS) -o 1612.05582/dE

ccbar: ccbar/ccbar.cpp PineAPPL.hpp config.h FO.hpp Integration.hpp Kernel.hpp MaunaKea.hpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) $(LHAPDF_DEPS) $(GSL_DEPS) $(DVEGAS_DEPS) -o ccbar/ccbar

MaunaKeaPy: Python.cpp PineAPPL.hpp config.h FO.hpp Integration.hpp Kernel.hpp MaunaKea.hpp
	$(CXX) $(CXXFLAGS) -shared -fPIC $< $(PINEAPPL_DEPS) $(LHAPDF_DEPS) $(GSL_DEPS) $(DVEGAS_DEPS) $(PY_DEPS) $(PYBIND_DEPS) -o MaunaKea`python3.10-config --extension-suffix`

PHONY: clean

clean:
	rm -f MaunaKea
	rm -f 1612.05582/dE
	rm -f ccbar/ccbar
