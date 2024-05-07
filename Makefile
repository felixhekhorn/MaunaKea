CXX = c++
CXXFLAGS = -std=c++14 -Wall -Wextra -O3 -IMaunaKea/
PINEAPPL_DEPS != pkg-config --cflags --libs pineappl_capi
LHAPDF_DEPS != pkg-config --cflags --libs lhapdf
GSL_DEPS != gsl-config --cflags --libs
DVEGAS_DEPS != pkg-config --cflags --libs dvegas
PY310_DEPS != python3.10-config --includes --ldflags
PYBIND_DEPS != python3 -m pybind11 --includes

FILES != ls MaunaKea/*.h*

py310: Python.cpp $(FILES)
	$(CXX) $(CXXFLAGS) -shared -fPIC $< $(PINEAPPL_DEPS) $(LHAPDF_DEPS) $(GSL_DEPS) $(DVEGAS_DEPS) $(PY310_DEPS) $(PYBIND_DEPS) -o MaunaKea`python3.10-config --extension-suffix`

PHONY: clean

clean:
	rm -f MaunaKea`python3.10-config --extension-suffix`
