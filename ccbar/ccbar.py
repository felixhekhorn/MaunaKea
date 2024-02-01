import argparse
import pathlib

import lhapdf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pineappl

TARGET_NAME = "MaunaKea-ccbar-fig1.pineappl.lz4"


def merge(ndata: int) -> None:
    """Merge grids."""
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
    # update metadata
    base.scale(1e-6)  # pico to µ
    base.set_key_value("x1_label", "sqrtS_h")
    base.set_key_value("x1_label_tex", r"\sqrt{S_h}")
    base.set_key_value("x1_unit", "GeV")
    base.set_key_value("y_label", "sigma_tot")
    base.set_key_value("y_label_tex", r"\sigma_{tot}")
    base.set_key_value("y_unit", "µb")
    # save back to disk
    merged_path = pathlib.Path(TARGET_NAME)
    print(f"Write combined grid to {merged_path}")
    base.write_lz4(merged_path)


def extract_sv_by_order(
    grid: pineappl.grid.Grid, central_pdf: lhapdf.PDF, pto: int
) -> pd.DataFrame:
    """Extract a given PTO together with its SV."""
    # compute central value
    order_mask = pineappl.grid.Order.create_mask(grid.orders(), 2 - 1 + pto, 0, True)
    central = grid.convolute_with_one(
        2212, central_pdf.xfxQ2, central_pdf.alphasQ2, order_mask=order_mask
    )
    # compute SV
    xis = []
    for xif in (0.5, 1.0, 2.0):
        for xir in (0.5, 1.0, 2.0):
            if xif / xir >= 4.0 or xir / xif >= 4.0:
                continue
            xis.append((xir, xif))
    # xis = [(0.5, 0.5), (2.0, 2.0)]
    sv_vals = (
        grid.convolute_with_one(
            2212, central_pdf.xfxQ2, central_pdf.alphasQ2, xi=xis, order_mask=order_mask
        )
    ).reshape(len(central), 7)
    df = pd.DataFrame()
    df["sqrt_s"] = grid.bin_left(0)
    df["central"] = central
    df["sv_min"] = np.min(sv_vals, axis=1)
    df["sv_max"] = np.max(sv_vals, axis=1)
    return df


def extract_lumis_by_order(
    grid: pineappl.grid.Grid, central_pdf: lhapdf.PDF, pto: int
) -> pd.DataFrame:
    """Extract lumi fraction for a given PTO."""
    order_mask = pineappl.grid.Order.create_mask(grid.orders(), 2 - 1 + pto, 0, True)
    full = grid.convolute_with_one(
        2212, central_pdf.xfxQ2, central_pdf.alphasQ2, order_mask=order_mask
    )
    df = pd.DataFrame()
    df["sqrt_s"] = grid.bin_left(0)
    df["total"] = full
    for lu, lab in enumerate(("gg", "qqbar", "gq")):
        lumi_mask = [False] * 6
        lumi_mask[lu] = True
        lumi = grid.convolute_with_one(
            2212,
            central_pdf.xfxQ2,
            central_pdf.alphasQ2,
            lumi_mask=lumi_mask,
            order_mask=order_mask,
        )
        df[lab] = lumi
    return df


# central_pdf = pdfs[0]
# central = grid.convolute_with_one(2212, central_pdf.xfxQ2, central_pdf.alphasQ2)
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


def plot_pto() -> None:
    """Plot convergence of PTO."""
    # prepare objects
    grid_path = pathlib.Path(TARGET_NAME)
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF("NNPDF40_nlo_pch_as_01180_nf_3", 0)
    grid = pineappl.grid.Grid.read(grid_path)
    # prepare data
    dfs = {}
    for k in range(2 + 1):
        df = extract_sv_by_order(grid, central_pdf, k)
        df.to_csv(f"ccbar-pto-{k}.csv")
        dfs[k] = df

    # plot bare
    fig, axs = plt.subplots(2, 1, height_ratios=[1, 0.35], sharex=True)
    for k, lab in [(0, "LO"), (1, "NLO"), (2, "NNLO")]:
        df = dfs[k]
        axs[0].fill_between(df["sqrt_s"], df["sv_min"], df["sv_max"], alpha=0.4)
        axs[0].plot(df["sqrt_s"], df["central"], label=lab)
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_ylabel(r"$\sigma_{c\bar{c}}$ [µb]")
    axs[0].tick_params(
        "both",
        which="both",
        direction="in",
        bottom=True,
        top=True,
        left=True,
        right=True,
    )
    axs[0].legend()
    # plot K-factor
    axs[1].plot([])  # add empty plot to align colors
    for k1, k2, lab in [(1, 0, "NLO/LO"), (2, 1, "NNLO/NLO")]:
        axs[1].plot(
            dfs[k1]["sqrt_s"], dfs[k1]["central"] / dfs[k2]["central"], label=lab
        )
    axs[1].set_xlabel(r"$\sqrt{s}$ [GeV]")
    axs[1].set_ylabel(r"K factor")
    axs[1].tick_params(
        "both",
        which="both",
        direction="in",
        bottom=True,
        top=True,
        left=True,
        right=True,
    )
    axs[1].legend()
    fig.tight_layout()
    fig.savefig("ccbar-pto.pdf")


def plot_lumi() -> None:
    """Plot lumi separation."""
    # prepare objects
    grid_path = pathlib.Path(TARGET_NAME)
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF("NNPDF40_nlo_pch_as_01180_nf_3", 0)
    grid = pineappl.grid.Grid.read(grid_path)
    # prepare data
    dfs = {}
    for k in range(2 + 1):
        df = extract_lumis_by_order(grid, central_pdf, k)
        df.to_csv(f"ccbar-lumi-{k}.csv")
        dfs[k] = df

    # plot
    fig, ax = plt.subplots(1, 1)
    for pto, lab_pto, ls in [
        (0, "LO", "dotted"),
        (1, "NLO", "dashed"),
        (2, "NNLO", "solid"),
    ]:
        df = dfs[pto]
        for lu, (raw_lab, plt_lab) in enumerate(
            [("gg", "$gg$"), ("qqbar", r"$q\bar{q}$"), ("gq", "$gq$")]
        ):
            if raw_lab == "gg" or pto == 2:
                lab = f"{plt_lab} ({lab_pto})"
            else:
                lab = None
            ax.plot(
                df["sqrt_s"],
                df[raw_lab] / df["total"] * 100,
                label=lab,
                color=f"C{lu}",
                linestyle=ls,
            )
    ax.set_xscale("log")
    ax.set_ylabel(r"$\sigma^{ij}/\sigma^{tot}$ [%]")
    ax.tick_params(
        "both",
        which="both",
        direction="in",
        bottom=True,
        top=True,
        left=True,
        right=True,
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig("ccbar-lumi.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="sub-command help", required=True)
    parser_merge = subparsers.add_parser("merge", help="Merge grids.")
    parser_merge.add_argument("ndata", help="number of points")
    parser_merge.set_defaults(mode="merge")
    parser_plot = subparsers.add_parser("plot-pto", help="Plot convergence of PTO.")
    parser_plot.set_defaults(mode="plot-pto")
    parser_plot = subparsers.add_parser("plot-lumi", help="Plot lumi separation.")
    parser_plot.set_defaults(mode="plot-lumi")
    args = parser.parse_args()
    if args.mode == "merge":
        merge(int(args.ndata))
    elif args.mode == "plot-pto":
        plot_pto()
    elif args.mode == "plot-lumi":
        plot_lumi()
