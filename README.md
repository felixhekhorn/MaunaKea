# MaunaKea

MaunaKea is a patched version of [top++](https://www.precision.hep.phy.cam.ac.uk/top-plus-plus/) by
M. Czakon and A. Mitov. It makes the number of light flavors dynamic and adds a
[PineAPPL](https://github.com/NNPDF/pineappl) interface - this is why it can also be called
"Hawaiian top++", but then the top of Hawai'i is [Mauna Kea](https://en.wikipedia.org/wiki/Mauna_Kea).

## External dependencies
- [LHAPDF](https://lhapdf.hepforge.org/) >= 6.4
- [Dvegas](https://dvegas.hepforge.org/) >= 2.0.3 with patched [`dvegas.h`](https://github.com/felixhekhorn/LeProHQ/blob/main/Patches/dvegas.h.patch) and [`dvegas.cpp`](https://github.com/felixhekhorn/LeProHQ/blob/main/Patches/dvegas.cpp.patch)
- [gsl](https://www.gnu.org/software/gsl/) >= 2.7
- [PineAPPL](https://github.com/NNPDF/pineappl) >= 0.6
- [pybind11](https://pybind11.readthedocs.io/en/latest/index.html) >= 2 for Python interface

## Citation policy
When using the code please cite
- this repository [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14185847.svg)](https://doi.org/10.5281/zenodo.14185847)
- the original paper of top++ [![arXiv](https://img.shields.io/badge/arXiv-1112.5675-b31b1b?labelColor=222222)](https://arxiv.org/abs/1112.5675) (+ if necessary the papers for the coefficient functions, see their website)
