"""Plotting routines."""

import argparse
import pathlib
from collections.abc import Callable, Collection, Mapping
from dataclasses import dataclass
from typing import Tuple

import lhapdf
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pineappl
from scipy.integrate import quad

from run import LABELS, MSHT20_MBRANGE, MSHT20_MCRANGE, PDFS, grid_path

TEX_LABELS = {3: r"c\bar{c}", 4: r"b\bar{b}"}


@dataclass(frozen=True)
class Extrapolation:
    """Extrapolation configuration"""

    x: float
    const: bool

    def masked_xfxQ2(self, pdf) -> Callable[[int, float, float], float]:
        """The masked xfxQ2 function."""

        def xfxQ2(pid: int, x: float, Q2: float):
            xmin = self.x if self.x >= 0.0 else float(pdf.set().get_entry("XMin"))
            if x < xmin:
                return pdf.xfxQ2(pid, xmin, Q2) if self.const else 0.0
            return pdf.xfxQ2(pid, x, Q2)

        return xfxQ2

    @property
    def suffix(self) -> str:
        """File name suffix."""
        if self.x == 0.0:
            return ""
        x = f"{self.x:f}" if self.x >= 0.0 else "n"
        return f"-ex{x}" + ("_c" if self.const else "")


# Set the default color cycle
mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
    color=["#001D66", "#B01C64", "#FFBF33", "#3D2E85", "#FF6000", "#FFF180"]
)


def m2str(m2: float) -> str:
    """Sanitize `.` to `p` for overleaf."""
    return f"{m2:.2f}".replace(".", "p")


def extract_sv_by_order(
    grid: pineappl.grid.Grid, central_pdf: lhapdf.PDF, extra: Extrapolation, pto_: int
) -> pd.DataFrame:
    """Extract a given PTO together with its SV."""
    # compute central value
    order_mask = pineappl.grid.Order.create_mask(grid.orders(), 2 - 1 + pto_, 0, True)
    central = grid.convolute_with_one(
        2212,
        extra.masked_xfxQ2(central_pdf),
        central_pdf.alphasQ2,
        order_mask=order_mask,
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
            2212,
            extra.masked_xfxQ2(central_pdf),
            central_pdf.alphasQ2,
            xi=xis,
            order_mask=order_mask,
        )
    ).reshape(len(central), 7)
    df = pd.DataFrame()
    df["sqrt_s"] = grid.bin_left(0)
    df["central"] = central
    df["sv_min"] = np.min(sv_vals, axis=1)
    df["sv_max"] = np.max(sv_vals, axis=1)
    return df


def load_pto(
    m2: float, nl: int, pdf: str, extra: Extrapolation
) -> Mapping[int, pd.DataFrame]:
    """Load PTO data."""
    # prepare objects
    grid_path_ = pathlib.Path(grid_path(m2, nl))
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF(pdf)
    grid = pineappl.grid.Grid.read(grid_path_)
    # prepare data
    dfs = {}
    for k in range(2 + 1):
        df = extract_sv_by_order(grid, central_pdf, extra, k)
        df.to_csv(f"data/{LABELS[nl]}-{m2:.2f}-{pdf}-pto-{k}{extra.suffix}.csv")
        dfs[k] = df
    return dfs


def load_xmean_pto(
    m2: float, nl: int, pdf: str, extra: Extrapolation
) -> Mapping[int, pd.DataFrame]:
    """Load PTO data."""
    # prepare objects
    grid_path_ = pathlib.Path(grid_path(m2, nl))
    lhapdf.setVerbosity(0)
    pdf_ = lhapdf.mkPDF(pdf)
    grid = pineappl.grid.Grid.read(grid_path_)
    # prepare data
    dfs = {}
    for k in range(2 + 1):
        order_mask = pineappl.grid.Order.create_mask(grid.orders(), 2 - 1 + k, 0, True)
        # <x> = <(x1+x2)/2>
        x1 = grid.convolute_with_two(
            2212,
            lambda pid, x, q2: x * extra.masked_xfxQ2(pdf_)(pid, x, q2),
            2212,
            extra.masked_xfxQ2(pdf_),
            pdf_.alphasQ2,
            order_mask=order_mask,
        )
        x2 = grid.convolute_with_two(
            2212,
            extra.masked_xfxQ2(pdf_),
            2212,
            lambda pid, x, q2: x * extra.masked_xfxQ2(pdf_)(pid, x, q2),
            pdf_.alphasQ2,
            order_mask=order_mask,
        )
        sig = grid.convolute_with_one(
            2212, extra.masked_xfxQ2(pdf_), pdf_.alphasQ2, order_mask=order_mask
        )
        mean = (x1 + x2) / (2.0 * sig)
        # Δx^2 = <(x-<x>)^2> = <x^2> - <x>^2
        x1x1 = grid.convolute_with_two(
            2212,
            lambda pid, x, q2: x * x * extra.masked_xfxQ2(pdf_)(pid, x, q2),
            2212,
            extra.masked_xfxQ2(pdf_),
            pdf_.alphasQ2,
            order_mask=order_mask,
        )
        x1x2 = grid.convolute_with_one(
            2212,
            lambda pid, x, q2: x * extra.masked_xfxQ2(pdf_)(pid, x, q2),
            pdf_.alphasQ2,
            order_mask=order_mask,
        )
        x2x2 = grid.convolute_with_two(
            2212,
            extra.masked_xfxQ2(pdf_),
            2212,
            lambda pid, x, q2: x * x * extra.masked_xfxQ2(pdf_)(pid, x, q2),
            pdf_.alphasQ2,
            order_mask=order_mask,
        )
        unc2 = (x1x1 + 2.0 * x1x2 + x2x2) / sig - mean**2.0
        df = pd.DataFrame()
        df["sqrt_s"] = grid.bin_left(0)
        df["mean"] = mean
        df["unc"] = np.sqrt(unc2)
        df.to_csv(f"data/{LABELS[nl]}-{m2:.2f}-{pdf}-xmean-pto-{k}{extra.suffix}.csv")
        dfs[k] = df
    return dfs


def load_extra(
    m2: float, nl: int, pdf: str, extras: Collection[Extrapolation]
) -> Collection[pd.DataFrame]:
    """Load extrapolation data."""
    dfs = []
    for extra in extras:

        def conv(grid, pdf_, extra=extra):
            return grid.convolute_with_one(
                2212, extra.masked_xfxQ2(pdf_), pdf_.alphasQ2
            )

        pdfs = load_pdf(m2, nl, [pdf], f"extra{extra.suffix}", conv)
        dfs.append(pdfs[pdf])
    return dfs


def extract_lumis_by_order(
    grid: pineappl.grid.Grid, central_pdf: lhapdf.PDF, extra: Extrapolation, pto_: int
) -> pd.DataFrame:
    """Extract lumi fraction for a given PTO."""
    order_mask = pineappl.grid.Order.create_mask(grid.orders(), 2 - 1 + pto_, 0, True)
    full = grid.convolute_with_one(
        2212,
        extra.masked_xfxQ2(central_pdf),
        central_pdf.alphasQ2,
        order_mask=order_mask,
    )
    df = pd.DataFrame()
    df["sqrt_s"] = grid.bin_left(0)
    df["total"] = full
    for lu, lab in enumerate(("gg", "qqbar", "gq")):
        lumi_mask = [False] * 6
        lumi_mask[lu] = True
        df[lab] = grid.convolute_with_one(
            2212,
            extra.masked_xfxQ2(central_pdf),
            central_pdf.alphasQ2,
            lumi_mask=lumi_mask,
            order_mask=order_mask,
        )
    return df


def load_lumi(
    m2: float, nl: int, pdf: str, extra: Extrapolation
) -> Mapping[int, pd.DataFrame]:
    """Load lumi data."""
    # prepare objects
    grid_path_ = pathlib.Path(grid_path(m2, nl))
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF(pdf)
    grid = pineappl.grid.Grid.read(grid_path_)
    # prepare data
    dfs = {}
    for k in range(2 + 1):
        df = extract_lumis_by_order(grid, central_pdf, extra, k)
        df.to_csv(
            f"data/{LABELS[nl]}-{m2:.2f}-{pdf.replace('/','__')}-lumi-{k}{extra.suffix}.csv"
        )
        dfs[k] = df
    return dfs


def load_pdf(
    m2: float, nl: int, pdf_sets: Collection[str], suffix: str, f: Callable
) -> Mapping[str, pd.DataFrame]:
    """Load PDF data from a grid."""
    # prepare objects
    grid_path_ = pathlib.Path(grid_path(m2, nl))
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
        # add exception for CT18NNLO_NF*
        if len(pdfs) == 1:
            df["central"] = pdf_vals[0]
            df["pdf_minus"] = pdf_vals[0]
            df["pdf_plus"] = pdf_vals[0]
        else:  # use lhapdf to do the uncertainty math
            pdf_errs = []
            # make exception for MSHT20nnlo_mXrange_nfY
            if pdf_set.errorType == "unknown":
                df["central"] = pdf_vals.mean(axis=0)
                df["pdf_minus"] = pdf_vals.min(axis=0)
                df["pdf_plus"] = pdf_vals.max(axis=0)
            else:
                for obs in pdf_vals.T:
                    # pdf_errs.append(pdf_set.uncertainty(obs,alternative=("NNPDF" in pdf_set_name)))
                    pdf_errs.append(pdf_set.uncertainty(obs))
                df["central"] = list(map(lambda u: u.central, pdf_errs))
                df["pdf_minus"] = list(map(lambda u: u.central - u.errminus, pdf_errs))
                df["pdf_plus"] = list(map(lambda u: u.central + u.errplus, pdf_errs))
        df.to_csv(f"data/{LABELS[nl]}-{m2:.2f}-{pdf_set_name}-{suffix}.csv")
        dfs[pdf_set_name] = df

    return dfs


def load_mass(
    ms: Collection[float], nl: int, pdf_set: str, extra: Extrapolation
) -> pd.DataFrame:
    """Load PDF data from a grid."""
    # prepare objects
    lhapdf.setVerbosity(0)
    df = pd.DataFrame()
    pdf_vals = []
    m_central = (MSHT20_MCRANGE if abs(nl) == 3 else MSHT20_MBRANGE)[0]
    for j, m in enumerate(ms):
        grid_path_ = pathlib.Path(grid_path(m, nl))
        grid = pineappl.grid.Grid.read(grid_path_)
        pdf = lhapdf.mkPDF(pdf_set, j)
        df["sqrt_s"] = grid.bin_left(0)
        # Fix µ to the central choice
        xi = m_central / m
        vals = grid.convolute_with_one(
            2212, extra.masked_xfxQ2(pdf), pdf.alphasQ2, xi=[(xi, xi)]
        )
        df[f"{j}"] = vals
        pdf_vals.append(vals)
    pdf_vals = np.array(pdf_vals)
    df["min"] = pdf_vals.min(0)
    df["max"] = pdf_vals.max(0)
    df.to_csv(f"data/{LABELS[nl]}-{pdf_set}{extra.suffix}.csv")
    return df


def pdf_raw(
    m2: float,
    nl: int,
    suffix: str,
    f: Callable,
    ylabel: str,
) -> Tuple[mpl.figure.Figure, str]:
    """Plot PDF dependence."""
    # load data
    elems = to_elems(m2, abs(nl))
    mass_label = f"{m2:.2f}-" if m2 > 0.0 else ""

    # collect datapoints
    dfs = {}
    for actual_m2_, pdfs in elems:
        for pdf, df in load_pdf(actual_m2_, nl, pdfs, suffix, f).items():
            # relabel if necessary
            new_label = pdf if m2 > 0.0 else f"{pdf} $m_h={actual_m2_:.2f}$ GeV"
            dfs[new_label] = df
    output = f"{LABELS[nl]}-{mass_label}{suffix}.pdf"

    fig, axs = plt.subplots(2, 1, height_ratios=[1, 0.5], sharex=True)
    # plot nominal x_min
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
    if m2 != 0.0:
        add_xmin(m2, axs[0])
    # rel. size
    norm = list(dfs.values())[0]["central"]
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
    path = f"plots/{output}"
    fig.savefig(path)
    return fig, path


def to_elems(m2: float, nl: int) -> Collection[Tuple[float, Collection[str]]]:
    """Collect default settings."""
    elems = []
    if m2 > 0.0:  # use a fixed mass
        if nl == 3:
            elems = [
                (
                    m2,
                    [
                        "NNPDF40_nnlo_pch_as_01180_nf_3",
                        # "NNPDF40_nnlo_as_01180",
                        # "NNPDF31_nnlo_as_0118",
                        # "NNPDF30_nnlo_as_0118",
                        # "NNPDF31_nnlo_pch_as_0118",
                        "MSHT20nnlo_nf3",
                        "CT18NNLO_rescaled_NF3",
                        # "MSHT20nnlo_mcrange_nf3"
                    ],
                )
            ]
        elif nl == 4:
            elems = [
                (
                    m2,
                    [
                        "NNPDF40_nnlo_as_01180_nf_4",
                        # "NNPDF31_nnlo_as_0118_nf_4"
                        "MSHT20nnlo_nf4",
                        "CT18NNLO_rescaled_NF4",
                        # "MSHT20nnlo_mbrange_nf4"
                    ],
                )
            ]
    else:  # use "dynamic" mass
        if nl == 3:
            elems = [
                (1.51, ["NNPDF40_nnlo_pch_as_01180_nf_3"]),
                (1.4, ["MSHT20nnlo_nf3"]),
                (1.3, ["CT18NNLO_rescaled_NF3"]),
            ]
        elif nl == 4:
            elems = [
                (4.92, ["NNPDF40_nnlo_as_01180_nf_4"]),
                (4.75, ["MSHT20nnlo_nf4", "CT18NNLO_rescaled_NF4"]),
            ]
    return elems


def pdf_obs(m2: float, nl: int, extra: Extrapolation) -> None:
    """Plot PDF dependence."""

    # plot the actual predictions
    def conv(grid, pdf_):
        return grid.convolute_with_one(2212, extra.masked_xfxQ2(pdf_), pdf_.alphasQ2)

    fig, path = pdf_raw(
        m2,
        nl,
        "pdf" + extra.suffix,
        conv,
        f"$\\sigma_{{{TEX_LABELS[abs(nl)]}}}$ [µb]",
    )
    if nl == 3:
        fig.axes[0].set_ylim(10, 2e5)
    fig.savefig(path)


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


def pdf_gluon(m2: float, nl: int, extra: Extrapolation) -> None:
    """Plot gluon(x_min) dependence."""

    def extract(grid, pdf_):
        res = []
        for b in range(len(grid.bin_left(0))):
            sg = grid.subgrid(0, b, 0)
            res.append(
                extra.masked_xfxQ2(pdf_)(21, np.min(sg.x1_grid()), sg.mu2_grid()[0].fac)
            )
        return res

    pdf_raw(
        m2,
        nl,
        "gluon" + extra.suffix,
        extract,
        r"$xg(x_{min})$",
    )


def pdf_xmean(m2: float, nl: int, extra: Extrapolation) -> None:
    """Plot x_mean PDF dependence."""

    def conv(grid, pdf_):
        x1 = grid.convolute_with_two(
            2212,
            lambda pid, x, q2: x * extra.masked_xfxQ2(pdf_)(pid, x, q2),
            2212,
            extra.masked_xfxQ2(pdf_),
            pdf_.alphasQ2,
        )
        x2 = grid.convolute_with_two(
            2212,
            extra.masked_xfxQ2(pdf_),
            2212,
            lambda pid, x, q2: x * extra.masked_xfxQ2(pdf_)(pid, x, q2),
            pdf_.alphasQ2,
        )
        sig = grid.convolute_with_one(2212, extra.masked_xfxQ2(pdf_), pdf_.alphasQ2)
        return (x1 + x2) / (2.0 * sig)

    pdf_raw(
        m2,
        nl,
        "xmean" + extra.suffix,
        conv,
        r"$\langle x\rangle$",
    )


def pdf_gg(m2: float, nl: int, extra: Extrapolation) -> None:
    """Plot gg(x_min) dependence."""

    def lumi_ker(z: float, pdf_, x, mu2):
        return (
            1.0
            / z
            * extra.masked_xfxQ2(pdf_)(21, z, mu2)
            * extra.masked_xfxQ2(pdf_)(21, x / z, mu2)
        )

    def pdf_lumi(pdf_, x, mu2):
        i, _e = quad(lumi_ker, x, 1.0, args=(pdf_, x, mu2), limit=100)
        return i

    def extract(grid, pdf_):
        res = []
        for b in range(len(grid.bin_left(0))):
            sg = grid.subgrid(0, b, 0)
            res.append(pdf_lumi(pdf_, np.min(sg.x1_grid()), sg.mu2_grid()[0].fac))
        return res

    pdf_raw(
        m2,
        nl,
        "gg" + extra.suffix,
        extract,
        r"$x L_{gg}(x_{min})$",
    )


# restrict sqrt(s) if requested
SHORT_RANGE_MAX: float = 20e3


def pto(m2: float, nl: int, pdf: str, extra: Extrapolation, short_range: bool) -> None:
    """Plot convergence of PTO."""
    # prepare data
    dfs = load_pto(m2, nl, pdf, extra)

    # plot bare
    fig, axs = plt.subplots(3, 1, height_ratios=[1, 0.35, 0.35], sharex=True)
    axP = axs[0]
    for k, lab in [(0, "LO"), (1, "NLO"), (2, "NNLO")]:
        df = dfs[k]
        axP.fill_between(df["sqrt_s"], df["sv_min"], df["sv_max"], alpha=0.4)
        axP.plot(df["sqrt_s"], df["central"], label=lab)
        axP.set_xlim(
            df["sqrt_s"].min(), SHORT_RANGE_MAX if short_range else df["sqrt_s"].max()
        )
    axP.set_xscale("log")
    axP.set_yscale("log")
    axP.set_ylabel(f"$\\sigma_{{{TEX_LABELS[abs(nl)]}}}$ [µb]")
    add_xmin(m2, axP)
    axP.legend()
    # plot uncertainty
    axU = axs[1]
    for k, lab in [(0, "LO"), (1, "NLO"), (2, "NNLO")]:
        df = dfs[k]
        axU.plot(
            df["sqrt_s"],
            np.max([df["central"] - df["sv_min"], df["sv_max"] - df["central"]], 0)
            / df["central"],
            label=lab,
        )
    if nl == 3:  # There seems to be a spike somewhere
        axU.set_ylim(0.5, 3.0)
    axU.set_ylabel("rel. Unc.")
    # plot K-factor
    axK = axs[2]
    axK.plot([])  # add empty plot to align colors
    for k1, k2, lab in [(1, 0, "NLO/LO"), (2, 1, "NNLO/NLO")]:
        axK.plot(dfs[k1]["sqrt_s"], dfs[k1]["central"] / dfs[k2]["central"], label=lab)
    axK.set_ylim(1.0, 3.1)
    axK.set_xlabel(r"$\sqrt{s}$ [GeV]")
    axK.set_ylabel(r"K factor")
    axK.legend()
    for ax in axs:
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
    sr_suffix = "-sr" if short_range else ""
    fig.savefig(
        f"plots/{LABELS[nl]}-{m2str(m2)}-{pdf}-pto{extra.suffix}{sr_suffix}.pdf"
    )


def xmean_pto(m2: float, nl: int, pdf: str, extra: Extrapolation) -> None:
    """Plot x_mean PTO dependence."""
    # prepare data
    dfs = load_xmean_pto(m2, nl, pdf, extra)

    # plot bare
    fig, axs = plt.subplots(2, 1, height_ratios=[1, 0.35], sharex=True)
    for k, lab in [(0, "LO"), (1, "NLO"), (2, "NNLO")]:
        df = dfs[k]
        axs[0].fill_between(
            df["sqrt_s"], df["mean"] - df["unc"], df["mean"] + df["unc"], alpha=0.4
        )
        axs[0].plot(df["sqrt_s"], df["mean"], label=lab)
        axs[0].set_xlim(df["sqrt_s"].min(), df["sqrt_s"].max())
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_ylabel(f"$\\langle x \\rangle_{{{TEX_LABELS[abs(nl)]}}}$")
    axs[0].tick_params(
        "both",
        which="both",
        direction="in",
        bottom=True,
        top=True,
        left=True,
        right=True,
    )
    add_xmin(m2, axs[0])
    axs[0].legend()
    # plot K-factor
    axs[1].fill_between([], [])  # add empty plot to align colors
    axs[1].plot([])
    for k1, k2, lab in [(1, 0, "NLO/LO"), (2, 1, "NNLO/NLO")]:
        df1, df2 = dfs[k1], dfs[k2]
        norm = df2["mean"]
        axs[1].fill_between(
            df1["sqrt_s"],
            (df1["mean"] - df1["unc"]) / norm,
            (df1["mean"] + df1["unc"]) / norm,
            alpha=0.4,
        )
        axs[1].plot(df1["sqrt_s"], df1["mean"] / norm, label=lab)
    axs[1].set_ylim(0, 3.0)
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
    fig.savefig(f"plots/{LABELS[nl]}-{m2:.2f}-{pdf}-xmean-pto{extra.suffix}.pdf")


def extra_dep(m2: float, nl: int, pdf: str) -> None:
    """Plot extrapolation dependency."""
    # prepare data
    extras = [
        Extrapolation(0, False),
        Extrapolation(-1, False),
        Extrapolation(1e-5, False),
    ]
    dfs = load_extra(m2, nl, pdf, extras)

    # plot xs
    fig, axs = plt.subplots(2, 1, height_ratios=[1, 0.35], sharex=True)
    for extra, df in zip(extras, dfs):
        axs[0].fill_between(df["sqrt_s"], df["pdf_minus"], df["pdf_plus"], alpha=0.4)
        lab = "LHAPDF"
        if extra.x != 0.0:
            xmin = r"x_{min}^*" if extra.x < 0.0 else f"{extra.x}"
            if np.isclose(extra.x, 1e-4):
                xmin = r"10^{-4}"
            elif np.isclose(extra.x, 1e-5):
                xmin = r"10^{-5}"
            lab = f"$f(x<{xmin})$=" + ("0" if not extra.const else f"f({xmin})")
        axs[0].plot(df["sqrt_s"], df["central"], label=lab)
        axs[0].set_xlim(df["sqrt_s"].min(), df["sqrt_s"].max())
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_ylabel(f"$\\sigma_{{{TEX_LABELS[abs(nl)]}}}$ [µb]")
    axs[0].tick_params(
        "both",
        which="both",
        direction="in",
        bottom=True,
        top=True,
        left=True,
        right=True,
    )
    add_xmin(m2, axs[0])
    axs[0].legend()
    # plot rel. PDF uncertainty
    norm = dfs[0]["central"]
    for df in dfs:
        axs[1].fill_between(
            df["sqrt_s"], df["pdf_plus"] / norm, df["pdf_minus"] / norm, alpha=0.4
        )
        axs[1].plot(df["sqrt_s"], df["central"] / norm)
    axs[1].set_xlabel(r"$\sqrt{s}$ [GeV]")
    axs[1].set_ylabel(r"rel. PDF unc.")
    axs[1].set_ylim(0.5, 1.5)
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
    fig.savefig(f"plots/{LABELS[nl]}-{m2:.2f}-{pdf}-extra.pdf")


def lumi(m2: float, nl: int, pdf: str, extra: Extrapolation, short_range: bool) -> None:
    """Plot lumi separation."""
    # prepare data
    dfs = load_lumi(m2, nl, pdf, extra)

    # plot
    fig, ax = plt.subplots(1, 1)
    for pto_, ls in [(0, "dotted"), (1, "dashed"), (2, "solid")]:
        df = dfs[pto_]
        for lu, raw_lab in enumerate(["gg", "qqbar", "gq"]):
            ax.plot(
                df["sqrt_s"],
                df[raw_lab] / df["total"] * 100,
                label=" ",
                color=f"C{lu}",
                linestyle=ls,
            )
        ax.set_xlim(
            df["sqrt_s"].min(), SHORT_RANGE_MAX if short_range else df["sqrt_s"].max()
        )
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
    add_xmin(m2, ax)
    # Fiddle with legend! Thanks https://stackoverflow.com/a/44072076
    h, l = ax.get_legend_handles_labels()
    empty = [plt.plot([], marker="", ls="")[0]]
    handles = np.array(h).reshape(3, 3)
    labels = np.array(l, dtype="U30").reshape(3, 3)
    handles = np.insert(handles, 0, empty * 3, 0)
    labels = np.insert(labels, 0, ["$gg$", r"$q\bar{q}$", "$gq$"], 0)
    handles = np.insert(handles, 0, empty * 4, 1)
    labels = np.insert(labels, 0, [" ", "LO", "NLO", "NNLO"], 1)
    leg = ax.legend(handles.flatten(), labels.flatten(), ncol=4)
    for vpack in leg._legend_handle_box.get_children():
        for hpack in vpack.get_children():
            hpack.get_children()[0].set_width(0)
    for vpack in leg._legend_handle_box.get_children()[:1]:
        for hpack in vpack.get_children():
            hpack.get_children()[0].set_width(0)
    fig.tight_layout()
    sr_suffix = "-sr" if short_range else ""
    fig.savefig(
        f"plots/{LABELS[nl]}-{m2str(m2)}-{pdf.replace('/','__')}-lumi{extra.suffix}{sr_suffix}.pdf"
    )


def mass(nl: int, extra: Extrapolation, short_range: bool) -> None:
    """Plot mass dependency."""
    # prepare data
    if abs(nl) == 3:
        ms = MSHT20_MCRANGE
        pdf = "MSHT20nnlo_mcrange_nf3"
    else:
        ms = MSHT20_MBRANGE
        pdf = "MSHT20nnlo_mbrange_nf4"
    df = load_mass(ms, nl, pdf, extra)
    label = f"{pdf}\n$m_{TEX_LABELS[abs(nl)][0]}={min(ms):.2f} - {max(ms):.2f}$ GeV\n$\\mu=2\\cdot{ms[0]:.2f}$ GeV"

    # plot bare
    fig, axs = plt.subplots(2, 1, height_ratios=[1, 0.35], sharex=True)
    axP = axs[0]
    axP.fill_between(df["sqrt_s"], df["min"], df["max"], alpha=0.4)
    axP.plot(df["sqrt_s"], df["0"], label=label)
    axP.set_xlim(
        df["sqrt_s"].min(), SHORT_RANGE_MAX if short_range else df["sqrt_s"].max()
    )
    axP.set_xscale("log")
    axP.set_yscale("log")
    axP.set_ylabel(f"$\\sigma_{{{TEX_LABELS[abs(nl)]}}}$ [µb]")
    axP.tick_params(
        "both",
        which="both",
        direction="in",
        bottom=True,
        top=True,
        left=True,
        right=True,
    )
    axP.legend(loc="lower right")
    # plot rel. uncertainty
    axU = axs[1]
    axU.fill_between(df["sqrt_s"], df["min"] / df["0"], df["max"] / df["0"], alpha=0.4)
    axU.plot(df["sqrt_s"], np.ones(len(df["0"])))
    axU.set_xlabel(r"$\sqrt{s}$ [GeV]")
    axU.set_ylabel(r"rel. Unc.")
    axU.set_ylim(1.0 / 2.1, 2.1)
    axU.set_yscale("log")
    axU.set_yticks(np.geomspace(0.5, 2.0, 3), ["0.5", "1.0", "2.0"])
    axU.set_yticks(np.geomspace(0.5, 2.0, 5), [], minor=True)
    axU.tick_params(
        "both",
        which="both",
        direction="in",
        bottom=True,
        top=True,
        left=True,
        right=True,
    )
    fig.tight_layout()
    sr_suffix = "-sr" if short_range else ""
    fig.savefig(f"plots/{LABELS[nl]}-mass{extra.suffix}{sr_suffix}.pdf")


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("m2", type=float, help="Mass of heavy quark")
    parser.add_argument("nl", type=int, help="Number of light flavors")
    h_pto = "Plot convergence with PTO"
    parser.add_argument("--pto", help=h_pto, action="store_true")
    h_lumi = "Plot lumi separation"
    parser.add_argument("--lumi", help=h_lumi, action="store_true")
    h_pdf = "Plot PDF dependence"
    parser.add_argument("--pdf", help=h_pdf, action="store_true")
    h_gluon = "Plot gluon(x_min)"
    parser.add_argument("--gluon", help=h_gluon, action="store_true")
    h_xmean = "Plot x_mean PDF dependence"
    parser.add_argument("--xmean", help=h_xmean, action="store_true")
    h_xmean_pto = "Plot x_mean PTO dependence"
    parser.add_argument("--xmean-pto", help=h_xmean_pto, action="store_true")
    h_gg = "Plot gg(x_min)"
    parser.add_argument("--gg", help=h_gg, action="store_true")
    h_mass = "Plot mass dependency"
    parser.add_argument("--mass", help=h_mass, action="store_true")
    h_extra = "Plot extrapolation dependency"
    parser.add_argument("--extra", help=h_extra, action="store_true")
    parser.add_argument("--all", help="Plot everything", action="store_true")
    # configuration parameter
    parser.add_argument("--pdf_set", help="PDF used for plots")
    parser.add_argument(
        "--extrapolate-xmin", type=float, default=0.0, help="Extrapolation region"
    )
    parser.add_argument(
        "--extrapolate-const",
        help="Extrapolate with f(x_extra) instead of 0",
        action="store_true",
    )
    parser.add_argument("--short-range", help="Restrict abscissa", action="store_true")
    args = parser.parse_args()
    m2_: float = float(args.m2)
    nl_: int = int(args.nl)
    extra_ = Extrapolation(float(args.extrapolate_xmin), bool(args.extrapolate_const))
    # multi PDF plots
    if args.pdf or args.all:
        print(h_pdf)
        pdf_obs(m2_, nl_, extra_)
    if args.gluon or args.all:
        print(h_gluon)
        pdf_gluon(m2_, nl_, extra_)
    if args.xmean or args.all:
        print(h_xmean)
        pdf_xmean(m2_, nl_, extra_)
    if args.gg or args.all:
        print(h_gg)
        pdf_gg(m2_, nl_, extra_)
    if args.mass or args.all:
        print(h_mass)
        mass(nl_, extra_, args.short_range)
    # single PDF plots
    if m2_ > 0:
        pdf = PDFS[abs(nl_)][f"{m2_:.2f}"]
    else:
        m2_, pdf = to_elems(m2_, abs(nl_))[0]
        pdf = pdf[0]
    if args.pdf_set:
        pdf = args.pdf_set
    if args.pto or args.all:
        print(h_pto)
        pto(m2_, nl_, pdf, extra_, args.short_range)
    if args.xmean_pto or args.all:
        print(h_xmean_pto)
        xmean_pto(m2_, nl_, pdf, extra_)
    if args.lumi or args.all:
        print(h_lumi)
        lumi(m2_, nl_, pdf, extra_, args.short_range)
    if args.extra or args.all:
        print(h_extra)
        extra_dep(m2_, nl_, pdf)


if __name__ == "__main__":
    main()
