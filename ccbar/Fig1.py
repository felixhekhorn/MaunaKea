import argparse
import pathlib

import lhapdf
import matplotlib.pyplot as plt
import numpy as np
import pineappl

# pico to µ
PB2MUB = 1e-6


def merge(ndata: int) -> None:
    """Merge grids"""
    # merge grids according to their c.o.m. energy
    base = None
    for j in range(ndata):
        grid_path = pathlib.Path(f"MaunaKea-ccbar-fig1-{j}.pineappl.lz4")
        grid = pineappl.grid.Grid.read(grid_path)
        # rebin by S_h
        sqrt_s_h = np.sqrt(float(grid.key_values()["MaunaKea/hadronicS"]))
        print(f"Merging j={j} with √S_h={sqrt_s_h} GeV")
        br = pineappl.bin.BinRemapper([1.0], [(sqrt_s_h, sqrt_s_h)])
        grid.set_remapper(br)
        # start or continue?
        if base is None:
            base = grid
        else:
            base.merge(grid)
    # done!
    merged_path = pathlib.Path("MaunaKea-ccbar-fig1.pineappl.lz4")
    print(f"Write combined grid to {merged_path}")
    base.write_lz4(merged_path)


def plot() -> None:
    """Plot grid"""
    # prepare objects
    grid_path = pathlib.Path("MaunaKea-ccbar-fig1.pineappl.lz4")
    lhapdf.setVerbosity(0)
    pdfs = lhapdf.mkPDFs("NNPDF40_nlo_pch_as_01180_nf_3")
    grid = pineappl.grid.Grid.read(grid_path)
    sqrt_s = grid.bin_left(0)
    # compute central value
    central_pdf = pdfs[0]
    central = (
        grid.convolute_with_one(2212, central_pdf.xfxQ2, central_pdf.alphasQ2) * PB2MUB
    )
    # # compute SV
    # xis = []
    # for xif in (0.5, 1.0, 2.0):
    #     for xir in (0.5, 1.0, 2.0):
    #         if xif / xir >= 4.0 or xir / xif >= 4.0:
    #             continue
    #         xis.append((xir * xi0, xif * xi0))
    # # xis = [(0.5*xi0, 0.5*xi0), (2.*xi0, 2.*xi0)]
    # sv_vals = (
    #     grid.convolute_with_one(2212, central_pdf.xfxQ2, central_pdf.alphasQ2, xi=xis)
    #     * PB2MUB
    # ).reshape(len(sqrt_s), 7)
    # sv_errs = np.array([sv - c for sv, c in zip(sv_vals, central)])
    # sv_err = np.sum(np.power(sv_errs, 2) / 7.0, axis=1)
    # # compute PDF uncert
    # pdf_vals = []
    # for pdf in pdfs[1:]:
    #     pdf_vals.append(
    #         grid.convolute_with_one(2212, pdf.xfxQ2, pdf.alphasQ2, xi=[(xi0, xi0)])
    #         * PB2MUB
    #         - central
    #     )
    # pdf_err = np.sum(np.power(pdf_vals, 2), axis=0)
    # # combine uncert
    # tot_err = np.sqrt(sv_err + pdf_err)

    # plot
    fig, ax = plt.subplots(1, 1)
    # ax.fill_between(sqrt_s, central - tot_err, central + tot_err, alpha=0.5)
    ax.plot(sqrt_s, central, color="red")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sqrt{s}$ (GeV)")
    # ax.set_xlim(35.0, 1e5)
    ax.set_yscale("log")
    ax.set_ylabel(r"$\sigma_{c\bar{c}}$ (µb)")
    # ax.set_ylim(1e-5, 5e5)
    ax.tick_params(
        "both",
        which="both",
        direction="in",
        bottom=True,
        top=True,
        left=True,
        right=True,
    )
    fig.tight_layout()
    fig.savefig("ccbar-fig1.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="sub-command help", required=True)
    parser_plot = subparsers.add_parser("plot", help="plot grid")
    parser_plot.set_defaults(mode="plot")
    parser_merge = subparsers.add_parser("merge", help="merge grids")
    parser_merge.add_argument("ndata", help="number of points")
    parser_merge.set_defaults(mode="merge")
    args = parser.parse_args()
    if args.mode == "merge":
        merge(int(args.ndata))
    elif args.mode == "plot":
        plot()
