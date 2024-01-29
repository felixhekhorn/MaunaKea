import argparse
import pathlib

import lhapdf
import numpy as np
import pineappl

# pico to milli
PB2MB = 1e-9


def main() -> None:
    """CLI entry point"""
    parser = argparse.ArgumentParser()
    parser.add_argument("nf", help="number of light flavors")
    parser.add_argument("sqrtS", help="square root of c.o.m. energy")
    args = parser.parse_args()

    nf = int(args.nf)
    sqrt_s = int(args.sqrtS)

    # prepare objects
    xi0 = 2.0
    q = "c" if nf == 3 else "b"
    lab = f"{q}{q}bar"
    grid_path = pathlib.Path(f"MaunaKea-{lab}-{sqrt_s}TeV.pineappl.lz4")
    lhapdf.setVerbosity(0)
    pdfs = lhapdf.mkPDFs(f"ABMP16_{nf}_nnlo")
    grid = pineappl.grid.Grid.read(grid_path)
    # compute central value
    central_pdf = pdfs[0]
    central = (
        grid.convolute_with_one(
            2212, central_pdf.xfxQ2, central_pdf.alphasQ2, xi=[(xi0, xi0)]
        )[0]
        * PB2MB
    )
    # compute SV
    xis = []
    for xif in (0.5, 1.0, 2.0):
        for xir in (0.5, 1.0, 2.0):
            if xif / xir >= 4.0 or xir / xif >= 4.0:
                continue
            xis.append((xir * xi0, xif * xi0))
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
            grid.convolute_with_one(2212, pdf.xfxQ2, pdf.alphasQ2, xi=[(xi0, xi0)])[0]
            * PB2MB
        )
    pdf_err = np.sqrt(np.sum(np.power(pdf_vals - central, 2)))
    # print
    if nf == 4 and sqrt_s == 14:
        print(
            f"sigma({lab}+X, sqrt(S)={sqrt_s}TeV) = {central:.2f}±{sv_err:.2f}_sc±{pdf_err:.2f}_pdf"
        )
    else:
        print(
            f"sigma({lab}+X, sqrt(S)={sqrt_s}TeV) = {central:.1f}±{sv_err:.1f}_sc±{pdf_err:.1f}_pdf"
        )


if __name__ == "__main__":
    main()
