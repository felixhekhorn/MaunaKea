"""Generate Figure 1."""
import argparse
import pathlib

import lhapdf
import matplotlib.pyplot as plt
import numpy as np
import pineappl

import MaunaKea

# pico to µ
PB2MUB = 1e-6

MASSES = {3: pow(1.67, 2), 4: pow(4.66, 2)}
NAME = {3: "c", 4: "b"}


def labels(nf: int, j: int) -> str:
    """Build label and grid path."""
    q = NAME[nf]
    lab = f"{q}{q}bar"
    return lab, f"fig1-{lab}-{j}.pineappl.lz4", f"fig1-{lab}.pineappl.lz4"


def compute(nl: int, ndata: int) -> None:
    """Compute grids."""
    m2: float = MASSES[nl]
    Sh_min: float = 35.0**2
    Sh_max: float = 100e3**2
    for j in range(ndata):
        logS_h: float = np.log(Sh_min) + (np.log(Sh_max) - np.log(Sh_min)) * j / (
            ndata - 1
        )
        S_h: float = np.exp(logS_h)
        print(f"j = {j:d}, sqrt(S) = {S_h:e}")
        # init object
        mk = MaunaKea.MaunaKea(m2, nl, MaunaKea.ORDER_ALL, MaunaKea.LUMI_ALL)
        mk.intCfg.calls = 50000
        mk.setHadronicS(S_h)
        mk.setPDF(f"ABMP16_{nl}_nnlo", 0)
        mk.setCentralScaleRatio(2.0)
        # fill the grid
        mk.run()
        int_out = mk.getIntegrationOutput()
        print(f"sigma_tot = {int_out.result:e} +- {int_out.error:e} [pb]\n")
        # save
        mk.write(labels(nl, j)[1])


def merge(nf: int, ndata: int) -> None:
    """Merge grids."""
    # merge grids according to their c.o.m. energy
    base = None
    for j in range(ndata):
        lab, pathj, path = labels(nf, j)
        grid_path = pathlib.Path(pathj)
        grid = pineappl.grid.Grid.read(grid_path)
        # rebin by S_h
        s_h = float(grid.key_values()["MaunaKea/hadronicS"])
        print(f"Merging {lab} j={j} with S_h={s_h} GeV^2")
        br = pineappl.bin.BinRemapper([1.0], [(s_h, s_h)])
        grid.set_remapper(br)
        # start or continue?
        if base is None:
            base = grid
        else:
            base.merge(grid)
    # done!
    merged_path = pathlib.Path(path)
    print(f"Write combined grid to {merged_path}")
    base.write_lz4(merged_path)


def plot(nf: int) -> None:
    """Plot grid."""
    # prepare objects
    lab, _pathj, path = labels(nf, 0)
    grid_path = pathlib.Path(path)
    lhapdf.setVerbosity(0)
    pdfs = lhapdf.mkPDFs(f"ABMP16_{nf}_nnlo")
    grid = pineappl.grid.Grid.read(grid_path)
    sqrt_s = np.sqrt(grid.bin_left(0))
    # compute central value
    central_pdf = pdfs[0]
    central = (
        grid.convolute_with_one(2212, central_pdf.xfxQ2, central_pdf.alphasQ2) * PB2MUB
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
        * PB2MUB
    ).reshape(len(sqrt_s), 7)
    sv_errs = np.array([sv - c for sv, c in zip(sv_vals, central)])
    sv_err = np.sum(np.power(sv_errs, 2) / 7.0, axis=1)
    # compute PDF uncert
    pdf_vals = []
    for pdf in pdfs[1:]:
        pdf_vals.append(
            grid.convolute_with_one(2212, pdf.xfxQ2, pdf.alphasQ2) * PB2MUB - central
        )
    pdf_err = np.sum(np.power(pdf_vals, 2), axis=0)
    # combine uncert
    tot_err = np.sqrt(sv_err + pdf_err)

    # plot
    fig, ax = plt.subplots(1, 1)
    ax.fill_between(sqrt_s, central - tot_err, central + tot_err, alpha=0.5)
    ax.plot(sqrt_s, central, color="red")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sqrt{s}$ (GeV)")
    ax.set_xlim(35.0, 1e5)
    ax.set_yscale("log")
    q = NAME[nf]
    ax.set_ylabel(f"$\\sigma_{{{q}\\bar{{{q}}}}}$ (µb)")
    ax.set_ylim(1e-5, 5e5)
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
    fig.savefig(f"fig1-{lab}.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="sub-command help", required=True)
    parser_compute = subparsers.add_parser("compute", help="compute grids")
    parser_compute.add_argument("nf", help="number of light flavors")
    parser_compute.add_argument("ndata", help="number of points")
    parser_compute.set_defaults(mode="compute")
    parser_merge = subparsers.add_parser("merge", help="merge grids")
    parser_merge.add_argument("nf", help="number of light flavors")
    parser_merge.add_argument("ndata", help="number of points")
    parser_merge.set_defaults(mode="merge")
    parser_plot = subparsers.add_parser("plot", help="plot grid")
    parser_plot.add_argument("nf", help="number of light flavors")
    parser_plot.set_defaults(mode="plot")
    args = parser.parse_args()
    mode = args.mode.strip().lower()
    if mode == "compute":
        compute(int(args.nf), int(args.ndata))
    elif mode == "merge":
        merge(int(args.nf), int(args.ndata))
    elif mode == "plot":
        plot(int(args.nf))
