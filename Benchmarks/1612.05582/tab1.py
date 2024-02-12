"""Generate Table 1."""
import argparse
import pathlib

import lhapdf
import numpy as np
import pineappl

import MaunaKea

# pico to milli
PB2MB = 1e-9

MASSES = {3: pow(1.67, 2), 4: pow(4.66, 2)}


def labels(nf: int, sqrt_s: int) -> str:
    """Build label and grid path."""
    q = "c" if nf == 3 else "b"
    lab = f"{q}{q}bar"
    return lab, f"tab1-{lab}-{sqrt_s}TeV.pineappl.lz4"


def compute(nf: int, sqrt_s: int) -> None:
    """Compute grid."""
    m2: float = MASSES[nf]
    # init object
    mk = MaunaKea.MaunaKea(m2, nf, MaunaKea.ORDER_ALL, MaunaKea.LUMI_ALL)
    mk.intCfg.calls = 50000
    mk.setHadronicS(sqrt_s * sqrt_s * 1e6)  # TeV to GeV
    mk.setPDF(f"ABMP16_{nf}_nnlo", 0)
    mk.setCentralScaleRatio(2.0)
    # fill the grid
    mk.run()
    int_out = mk.getIntegrationOutput()
    print(f"sigma_tot = {int_out.result:e} +- {int_out.error:e} [pb]\n")
    # save
    mk.write(labels(nf, sqrt_s)[1])


def plot(nf: int, sqrt_s: int) -> None:
    """Plot table."""
    # prepare objects
    lab, path = labels(nf, sqrt_s)
    grid_path = pathlib.Path(path)
    lhapdf.setVerbosity(0)
    pdfs = lhapdf.mkPDFs(f"ABMP16_{nf}_nnlo")
    grid = pineappl.grid.Grid.read(grid_path)
    # compute central value
    central_pdf = pdfs[0]
    central = (
        grid.convolute_with_one(2212, central_pdf.xfxQ2, central_pdf.alphasQ2)[0]
        * PB2MB
    )
    # compute SV
    xis = []
    for xif in (0.5, 1.0, 2.0):
        for xir in (0.5, 1.0, 2.0):
            if xif / xir >= 4.0 or xir / xif >= 4.0:
                continue
            xis.append((xir, xif))
    # xis = [(0.5*xi0, 0.5*xi0), (2.*xi0, 2.*xi0)]
    sv_vals = (
        grid.convolute_with_one(2212, central_pdf.xfxQ2, central_pdf.alphasQ2, xi=xis)
        * PB2MB
    )
    sv_err = np.max(np.abs(sv_vals - central))
    # compute PDF uncert
    pdf_vals = []
    for pdf in pdfs[1:]:
        pdf_vals.append(
            grid.convolute_with_one(2212, pdf.xfxQ2, pdf.alphasQ2)[0] * PB2MB
        )
    pdf_err = np.sqrt(np.sum(np.power(pdf_vals - central, 2)))
    # print
    if nf == 4 and sqrt_s == 14:
        print(
            f"sigma({lab}+X, sqrt(S)={sqrt_s}TeV) = {central:.2f} ± {sv_err:.2f}_sc ± {pdf_err:.2f}_pdf"
        )
    else:
        print(
            f"sigma({lab}+X, sqrt(S)={sqrt_s}TeV) = {central:.1f} ± {sv_err:.1f}_sc ± {pdf_err:.1f}_pdf"
        )


def main() -> None:
    """CLI entry point"""
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", help="compute or plot")
    parser.add_argument("nf", help="number of light flavors")
    parser.add_argument("sqrtS", help="square root of c.o.m. energy [TeV]")
    args = parser.parse_args()

    nf = int(args.nf)
    sqrt_s = int(args.sqrtS)
    mode = args.mode.strip().lower()
    if mode == "compute":
        compute(nf, sqrt_s)
    elif mode == "plot":
        plot(nf, sqrt_s)
    else:
        raise ValueError(f"mode has to be compute or plot, but was {mode}")


if __name__ == "__main__":
    main()
