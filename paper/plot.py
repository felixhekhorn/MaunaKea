"""Plotting routines."""

import argparse
import pathlib
from collections.abc import Callable, Collection, Mapping

import lhapdf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pineappl
from scipy.integrate import quad

from run import LABELS, MASSES, PDFS, grid_path

TEX_LABELS = {3: r"c\bar{c}", 4: r"b\bar{b}"}


PDF_SET_NAMES = {
    3: [
        "NNPDF40_nnlo_pch_as_01180_nf_3",
        "MSHT20nnlo_nf3",
        "ABMP16_3_nnlo",
        "CT18NNLO_NF3",
    ],
    4: [
        "NNPDF40_nnlo_as_01180_nf_4",
        "MSHT20nnlo_nf4",
        "ABMP16_4_nnlo",
        "CT18NNLO_NF4",
    ],
}


def extract_sv_by_order(
    grid: pineappl.grid.Grid, central_pdf: lhapdf.PDF, pto_: int
) -> pd.DataFrame:
    """Extract a given PTO together with its SV."""
    # compute central value
    order_mask = pineappl.grid.Order.create_mask(grid.orders(), 2 - 1 + pto_, 0, True)
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


def load_pto(nl: int) -> Mapping[int, pd.DataFrame]:
    """Load PTO data."""
    # prepare objects
    grid_path_ = pathlib.Path(grid_path(nl))
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF(PDFS[nl], 0)
    grid = pineappl.grid.Grid.read(grid_path_)
    # prepare data
    dfs = {}
    for k in range(2 + 1):
        df = extract_sv_by_order(grid, central_pdf, k)
        df.to_csv(f"data/{LABELS[nl]}-pto-{k}.csv")
        dfs[k] = df
    return dfs


def extract_lumis_by_order(
    grid: pineappl.grid.Grid, central_pdf: lhapdf.PDF, pto_: int
) -> pd.DataFrame:
    """Extract lumi fraction for a given PTO."""
    order_mask = pineappl.grid.Order.create_mask(grid.orders(), 2 - 1 + pto_, 0, True)
    full = grid.convolute_with_one(
        2212, central_pdf.xfxQ2, central_pdf.alphasQ2, order_mask=order_mask
    )
    df = pd.DataFrame()
    df["sqrt_s"] = grid.bin_left(0)
    df["total"] = full
    for lu, lab in enumerate(("gg", "qqbar", "gq")):
        lumi_mask = [False] * 6
        lumi_mask[lu] = True
        df[lab] = grid.convolute_with_one(
            2212,
            central_pdf.xfxQ2,
            central_pdf.alphasQ2,
            lumi_mask=lumi_mask,
            order_mask=order_mask,
        )
    return df


def load_lumi(nl: int, pdf_set: str) -> Mapping[int, pd.DataFrame]:
    """Load lumi data."""
    # prepare objects
    grid_path_ = pathlib.Path(grid_path(nl))
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF(pdf_set, 0)
    grid = pineappl.grid.Grid.read(grid_path_)
    # prepare data
    dfs = {}
    for k in range(2 + 1):
        df = extract_lumis_by_order(grid, central_pdf, k)
        df.to_csv(f"data/{LABELS[nl]}-lumi-{pdf_set}-{k}.csv")
        dfs[k] = df
    return dfs


def load_pdf(
    nl: int, pdf_sets: Collection[str], f: Callable
) -> Mapping[str, pd.DataFrame]:
    """Load PDF data."""
    # prepare objects
    grid_path_ = pathlib.Path(grid_path(nl))
    grid = pineappl.grid.Grid.read(grid_path_)
    lhapdf.setVerbosity(0)
    dfs = {}
    for pdf_set_name in pdf_sets:
        pdf_set = lhapdf.getPDFSet(pdf_set_name)
        pdfs = pdf_set.mkPDFs()
        # compute PDF uncert
        pdf_vals = []
        for pdf_ in pdfs:
            pdf_vals.append(f(grid, pdf_))
        pdf_vals = np.array(pdf_vals)
        df = pd.DataFrame()
        df["sqrt_s"] = grid.bin_left(0)
        # add exception for for CT18NNLO_NF*
        if len(pdfs) == 1:
            df["central"] = pdf_vals[0]
            df["pdf_minus"] = pdf_vals[0]
            df["pdf_plus"] = pdf_vals[0]
            dfs[pdf_set_name] = df
            continue
        # else use lhapdf to do the math
        pdf_errs = []
        for obs in pdf_vals.T:
            pdf_errs.append(pdf_set.uncertainty(obs))
        df["central"] = list(map(lambda u: u.central, pdf_errs))
        df["pdf_minus"] = list(map(lambda u: u.central - u.errminus, pdf_errs))
        df["pdf_plus"] = list(map(lambda u: u.central + u.errplus, pdf_errs))
        df.to_csv(f"data/{LABELS[nl]}-pdf-{pdf_set_name}.csv")
        dfs[pdf_set_name] = df

    return dfs


def pdf_raw(
    nl: int, f: Callable, ylabel: str, suffix: str, update_ax0: Callable = None
) -> None:
    """Plot PDF dependence."""
    # prepare data
    pdf_set_names = PDF_SET_NAMES[nl]
    dfs = load_pdf(nl, pdf_set_names, f)

    fig, axs = plt.subplots(2, 1, height_ratios=[1, 0.5], sharex=True)
    # plot nominal x_min
    for j, pdf_set_name in enumerate(pdf_set_names):
        xmin = lhapdf.getPDFSet(pdf_set_name).get_entry("XMin")
        m2 = MASSES[nl]
        for ax in axs:
            ax.axvline(
                x=_xmin2sqrts(m2, float(xmin)), c=f"C{j}", linestyle="--", alpha=0.3
            )
    # plot data
    for pdf_set, df in dfs.items():
        axs[0].fill_between(df["sqrt_s"], df["pdf_minus"], df["pdf_plus"], alpha=0.4)
        axs[0].plot(df["sqrt_s"], df["central"], label=pdf_set)
        axs[0].set_xlim(df["sqrt_s"].min(), df["sqrt_s"].max())
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_ylabel(ylabel)
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
    if update_ax0:
        update_ax0(axs[0])
    # rel. size
    norm = dfs[pdf_set_names[0]]["central"]
    for _, df in dfs.items():
        axs[1].fill_between(
            df["sqrt_s"], df["pdf_plus"] / norm, df["pdf_minus"] / norm, alpha=0.4
        )
        axs[1].plot(df["sqrt_s"], df["central"] / norm)
    axs[1].set_ylim(-0.5, 2)
    axs[1].set_xlabel(r"$\sqrt{s}$ [GeV]")
    axs[1].set_ylabel(r"rel. PDF unc.")
    axs[1].tick_params(
        "both",
        which="both",
        direction="in",
        bottom=True,
        top=True,
        left=True,
        right=True,
    )
    fig.tight_layout()
    fig.savefig(f"plots/{LABELS[nl]}-{suffix}.pdf")


def pdf_obs(nl: int) -> None:
    """Plot PDF dependence."""
    pdf_raw(
        nl,
        lambda grid, pdf_: grid.convolute_with_one(2212, pdf_.xfxQ2, pdf_.alphasQ2),
        f"$\\sigma_{{{TEX_LABELS[nl]}}}$ [µb]",
        "pdf",
    )


def _sqrts2xmin(m2: float, sqrt_s: float) -> float:
    r"""Convert :math:`\sqrt s` to :math:`x_{min}`."""
    return 4.0 * m2 / np.power(sqrt_s, 2)


def _xmin2sqrts(m2: float, xmin: float) -> float:
    r"""Convert :math:`x_{min}` to :math:`\sqrt s`."""
    return np.sqrt(4.0 * m2 / xmin)


def add_xmin(m2: float, ax0) -> None:
    """Add x_min as second x-axis on top."""
    ax0.tick_params(
        "both",
        which="both",
        direction="in",
        bottom=True,
        top=False,
        left=True,
        right=True,
    )
    secax = ax0.secondary_xaxis(
        "top",
        functions=(
            lambda sqrt_s: _sqrts2xmin(m2, sqrt_s),
            lambda xmin: _xmin2sqrts(m2, xmin),
        ),
    )
    secax.set_xlabel(r"$x_{min}$")
    secax.tick_params("x", which="both", direction="in")


def pdf_gluon(nl: int) -> None:
    """Plot gluon(x_min) dependence."""

    def extract(grid, pdf_):
        res = []
        for b in range(len(grid.bin_left(0))):
            sg = grid.subgrid(0, b, 0)
            res.append(pdf_.xfxQ2(21, np.min(sg.x1_grid()), sg.mu2_grid()[0].fac))
        return res

    m2 = MASSES[nl]

    pdf_raw(nl, extract, r"$xg(x_{min})$", "gluon", lambda ax0: add_xmin(m2, ax0))


def pdf_gg(nl: int) -> None:
    """Plot gg(x_min) dependence."""

    def lumi_ker(z: float, pdf_, x, mu2):
        return 1.0 / z * pdf_.xfxQ2(21, z, mu2) * pdf_.xfxQ2(21, x / z, mu2)

    def pdf_lumi(pdf_, x, mu2):
        i, _e = quad(lumi_ker, x, 1.0, args=(pdf_, x, mu2), limit=100)
        return i

    def extract(grid, pdf_):
        res = []
        for b in range(len(grid.bin_left(0))):
            sg = grid.subgrid(0, b, 0)
            res.append(pdf_lumi(pdf_, np.min(sg.x1_grid()), sg.mu2_grid()[0].fac))
        return res

    m2 = MASSES[nl]

    pdf_raw(nl, extract, r"$x L_{gg}(x_{min})$", "gg", lambda ax0: add_xmin(m2, ax0))


def pto(nl: int) -> None:
    """Plot convergence of PTO."""
    # prepare data
    dfs = load_pto(nl)

    # plot bare
    fig, axs = plt.subplots(2, 1, height_ratios=[1, 0.35], sharex=True)
    for k, lab in [(0, "LO"), (1, "NLO"), (2, "NNLO")]:
        df = dfs[k]
        axs[0].fill_between(df["sqrt_s"], df["sv_min"], df["sv_max"], alpha=0.4)
        axs[0].plot(df["sqrt_s"], df["central"], label=lab)
        axs[0].set_xlim(df["sqrt_s"].min(), df["sqrt_s"].max())
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_ylabel(f"$\\sigma_{{{TEX_LABELS[nl]}}}$ [µb]")
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
    axs[1].set_ylim(1.0, 4.0)
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
    fig.savefig(f"plots/{LABELS[nl]}-pto.pdf")


def lumi(nl: int, pdf_set: str) -> None:
    """Plot lumi separation."""
    # prepare data
    dfs = load_lumi(nl, pdf_set)

    # plot
    fig, ax = plt.subplots(1, 1)
    for pto_, lab_pto, ls in [
        (0, "LO", "dotted"),
        (1, "NLO", "dashed"),
        (2, "NNLO", "solid"),
    ]:
        df = dfs[pto_]
        for lu, (raw_lab, plt_lab) in enumerate(
            [("gg", "$gg$"), ("qqbar", r"$q\bar{q}$"), ("gq", "$gq$")]
        ):
            if raw_lab == "gg" or pto_ == 2:
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
        ax.set_xlim(df["sqrt_s"].min(), df["sqrt_s"].max())
    ax.set_ylim(-25, 107)
    ax.set_xscale("log")
    ax.set_ylabel(r"$\sigma^{ij}/\sigma^{tot}$ [%]")
    ax.set_xlabel(r"$\sqrt{s}$ [GeV]")
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
    fig.savefig(f"plots/{LABELS[nl]}-lumi-{pdf_set}.pdf")


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("nl", help="Number of light flavors")
    h_pto = "Plot convergence with PTO"
    parser.add_argument("--pto", help=h_pto, action="store_true")
    h_lumi = "Plot lumi separation"
    parser.add_argument("--lumi", help=h_lumi, action="store_true")
    h_pdf = "Plot PDF dependence"
    parser.add_argument("--pdf", help=h_pdf, action="store_true")
    h_gluon = "Plot gluon(x_min)"
    parser.add_argument("--gluon", help=h_gluon, action="store_true")
    h_gg = "Plot gg(x_min)"
    parser.add_argument("--gg", help=h_gg, action="store_true")
    parser.add_argument("--all", help="Plot everything", action="store_true")
    parser.add_argument(
        "--iterate-pdf-sets", help="Iterate on PDF sets.", action="store_true"
    )
    parser.add_argument("--pdf-set", help="PDF set used for plots.")
    args = parser.parse_args()
    nl_ = int(args.nl)
    pdf_sets = [PDFS[nl_]]
    if args.iterate_pdf_sets:
        pdf_sets = PDF_SET_NAMES[nl_]
    if args.pdf_set:
        pdf_sets = [args.pdf_set]
    if args.pto or args.all:
        print(h_pto)
        pto(nl_)
    if args.lumi or args.all:
        print(h_lumi)
        for ps in pdf_sets:
            lumi(nl_, ps)
    if args.pdf or args.all:
        print(h_pdf)
        pdf_obs(nl_)
    if args.gluon or args.all:
        print(h_gluon)
        pdf_gluon(nl_)
    if args.gg or args.all:
        print(h_gg)
        pdf_gg(nl_)


if __name__ == "__main__":
    main()
