CXX = c++
CXXFLAGS = -Wall -Wextra -O2
HATHORPATH = ./../../../Hathor-2.0
HATHOR_DEPS = -I$(HATHORPATH)/include -L $(HATHORPATH)/lib -lHathor -lff -lgfortran -lm
LHAPDF_DEPS != pkg-config --cflags --libs lhapdf

run: run.cpp
	$(CXX) $(CXXFLAGS) $< $(HATHOR_DEPS) $(LHAPDF_DEPS) -o run.exe

PHONY: clean

clean:
	rm -f run
