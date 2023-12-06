# MaunaKea

MaunaKea is a patched version of [top++](https://www.precision.hep.phy.cam.ac.uk/top-plus-plus/) by
M. Czakon and A. Mitov. It makes the number of light flavors dynamic and adds a
[PineAPPL](https://github.com/NNPDF/pineappl) interface - this is why it can also be called
"Hawaiian top++", but then the top of Hawai'i is [Mauna Kea](https://en.wikipedia.org/wiki/Mauna_Kea).

## External dependencies
- [LHAPDF](https://lhapdf.hepforge.org/)
- [Dvegas](https://dvegas.hepforge.org/) >= 2.0.3 with patched [`dvegas.h`](https://github.com/felixhekhorn/LeProHQ/blob/main/Patches/dvegas.h.patch) and [`dvegas.cpp`](https://github.com/felixhekhorn/LeProHQ/blob/main/Patches/dvegas.cpp.patch)
- [PineAPPL](https://github.com/NNPDF/pineappl)
