"""Compute grids."""

import argparse
import os
import pathlib
import time
from dataclasses import dataclass
from typing import Optional

import lhapdf
import numpy as np
import pineappl

import MaunaKea


@dataclass(frozen=True)
class ExpConfig:
    """Experiment config"""

    no: str
    """Reference number."""
    sqrts: float
    """Center-of-mass energy in GeV."""
    correction: Optional[int] = None
    """Correction to compute."""


A: dict[str, int] = {"S": 32, "Au": 197, "Pb": 208}
"""Map element names to their closest integer A."""

DATA: dict[int, list[ExpConfig]] = {
    3: [
        ExpConfig("55", 0.0114e3),
        ExpConfig("48", 0.0217e3),
        ExpConfig("48", 0.026e3),
        ExpConfig("48", 0.0274e3),
        ExpConfig("48", 0.0289e3),
        ExpConfig("48", 0.0389e3),
        ExpConfig("47", 0.0416e3),
        ExpConfig("46", 0.0685e3),
        ExpConfig("45,55", 0.1104e3),
        ExpConfig("37,41,42,44", 0.2e3),
        ExpConfig("34", 1.96e3),
        ExpConfig("19", 2.76e3),
        ExpConfig("31", 5e3),
        ExpConfig("20,32,92", 5e3, A["Pb"]),
        ExpConfig("94,97,98", 5.02e3),
        ExpConfig("15,17,21,23", 7e3),
        ExpConfig("93,95", 8.16e3, A["Pb"]),
        ExpConfig("11", 13e3),
    ],
    4: [
        # ExpConfig("89", 0.0387e3, A["S"]), !no EPPS PDFs for S!
        ExpConfig("90", 0.0387e3, A["Au"]),
        ExpConfig("88", 0.0416e3),
        ExpConfig("43,84", 0.2e3),
        ExpConfig("82", 0.51e3),
        ExpConfig("80", 0.63e3, -1),
        ExpConfig("77", 1.96e3),
        ExpConfig("76", 2.76e3),
        ExpConfig("75", 5.02e3, A["Pb"]),
        ExpConfig("15,16,56", 7e3),
        ExpConfig("56", 13e3),
    ],
}
"""Experiment configurations for c and b."""

LABELS = {3: "ccbar", 4: "bbbar", -3: "data-ccbar", -4: "data-bbbar"}
"""File prefixes."""

SH_MIN: float = 20.0**2
"""Minimum energy for plotting."""

SH_MAX: float = 400e3**2
"""Maximum energy for plotting."""

CALLS: int = 50000
"""MC calls for full run."""

# setup default PDFs
PDFS = {
    3: {
        "1.51": "NNPDF40_nnlo_pch_as_01180_nf_3",
        "1.40": "MSHT20nnlo_nf3",
        "1.30": "CT18NNLO_NF3",
    },
    4: {
        "4.92": "NNPDF40_nnlo_as_01180_nf_4",
        "4.75": "MSHT20nnlo_nf4",  # or "CT18NNLO_NF4",
    },
}
"""Default PDF choices."""

MSHT20_MCRANGE = [1.4, 1.2, 1.25, 1.3, 1.35, 1.45, 1.5, 1.55]
"""Charm mass range in MSHT20."""
for j_, msht_mc in enumerate(MSHT20_MCRANGE[1:]):
    PDFS[3][f"{msht_mc:.2f}"] = f"MSHT20nnlo_mcrange_nf3/{j_+1}"

MSHT20_MBRANGE = [4.75, 4.0, 4.25, 4.5, 5.0, 5.25, 5.5]
"""Bottom mass range in MSHT20."""
for j_, msht_mb in enumerate(MSHT20_MBRANGE[1:]):
    PDFS[4][f"{msht_mb:.2f}"] = f"MSHT20nnlo_mbrange_nf4/{j_+1}"


def subgrid_path(mass: float, nf: int, j: int) -> str:
    """Return grid path for single configuration."""
    return f"subgrids/{LABELS[nf]}-{mass:.2f}-{j}.pineappl.lz4"


def grid_path(mass: float, nf: int) -> str:
    """Return combined grid path."""
    return f"{LABELS[nf]}-{mass:.2f}.pineappl.lz4"


def compute(
    ndata: int, mass: float, nl: int, pdf: str, processes: int = 1, quick: bool = False
) -> None:
    """Compute grids."""
    # determine energy range
    if nl < 0:
        sh_range = [c.sqrts**2.0 for c in DATA[abs(nl)]]
        ndata = len(sh_range)
    else:
        sh_range = np.geomspace(SH_MIN, SH_MAX, ndata)
    # parallelize if requested
    if processes <= 0:
        processes = max(os.cpu_count() + processes, 1)
    start = time.perf_counter()
    print(
        f"Computing with mass={mass} GeV, nl={abs(nl)}, pdf={pdf} using {processes} threads"
    )
    args = zip(
        range(ndata),
        [ndata] * ndata,
        [mass] * ndata,
        [nl] * ndata,
        sh_range,
        [pdf] * ndata,
        [quick] * ndata,
    )
    if processes > 1:
        from multiprocessing import Pool  # pylint: disable=import-outside-toplevel

        with Pool(processes) as p:
            p.starmap(
                compute_subgrid,
                args,
            )
    else:
        for arg in args:
            compute_subgrid(*arg)
    delta = time.perf_counter() - start
    print("---")
    print(f"computed {ndata} grids in {delta/60:.2f} min")


def compute_subgrid(
    j: int, ndata: int, mass: float, nl: int, sh: float, pdf: str, quick: bool
) -> None:
    """Compute a subgrid."""
    print(f"j = {j:d}/{ndata:d}, S = {sh:e} GeV^2")
    lhapdf.setVerbosity(0)
    start = time.perf_counter()
    # init object
    mk = MaunaKea.MaunaKea(mass * mass, abs(nl), MaunaKea.ORDER_ALL, MaunaKea.LUMI_ALL)
    mk.intCfg.calls = 5000 if quick else CALLS
    mk.setHadronicS(sh)
    mk.setPDF(pdf)
    mk.setCentralScaleRatio(2.0)
    # fill the grid
    mk.run()
    int_out = mk.getIntegrationOutput()
    delta = time.perf_counter() - start
    print(f"sigma_tot = {int_out.result:e} +- {int_out.error:e} [pb]")
    print(f"took {delta:.2f} s and chi2iter={int_out.chi2iter}")
    # save
    mk.write(subgrid_path(mass, nl, j))


def merge(ndata: int, mass: float, nl: int) -> None:
    """Merge grids."""
    # merge grids according to their c.o.m. energy
    base = None
    if nl < 0:
        ndata = len(DATA[abs(nl)])
    for j in range(ndata):
        subgrid_path_ = pathlib.Path(subgrid_path(mass, nl, j))
        grid = pineappl.grid.Grid.read(subgrid_path_)
        # rebin by S_h
        sqrt_sh = np.sqrt(float(grid.key_values()["MaunaKea::hadronicS"]))
        print(f"Merging j={j} with √S_h={sqrt_sh} GeV")
        br = pineappl.bin.BinRemapper([1.0], [(sqrt_sh, sqrt_sh)])
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
    merged_path = pathlib.Path(grid_path(mass, nl))
    print(f"Write combined grid to {merged_path}")
    base.write_lz4(merged_path)


def cli() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("mass", help="Mass of heavy quark")
    parser.add_argument("nl", help="Number of light flavors")
    parser.add_argument("ndata", help="Number of points")
    parser.add_argument("-c", "--compute", help="Compute grids", action="store_true")
    parser.add_argument("-m", "--merge", help="Merge grids", action="store_true")
    parser.add_argument("--pdf", help="PDF set used for computing")
    parser.add_argument(
        "--processes", default=1, help="Number of parallel threads for computing"
    )
    parser.add_argument("--quick", help="Use low statistics", action="store_true")
    # prepare args
    args = parser.parse_args()
    mass_: float = float(args.mass)
    nl_: int = int(args.nl)
    # determine PDF
    pdf_ = None
    if args.pdf:
        pdf_ = args.pdf
    else:
        pdf_ = PDFS[abs(nl_)].get(f"{mass_:.2f}")
    if pdf_ is None or len(pdf_.strip()) <= 0:
        raise ValueError("No PDF set given!")
    # do something
    if args.compute:
        compute(
            int(args.ndata), mass_, nl_, pdf_, int(args.processes), bool(args.quick)
        )
    if args.merge:
        merge(int(args.ndata), mass_, nl_)


if __name__ == "__main__":
    cli()
