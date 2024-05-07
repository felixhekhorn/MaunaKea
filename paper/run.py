"""Compute grids."""

import argparse
import os
import pathlib
import time
from multiprocessing import Pool

import lhapdf
import numpy as np
import pineappl

import MaunaKea

LABELS = {3: "ccbar", 4: "bbbar"}
SH_MIN: float = 20.0**2
SH_MAX: float = 400e3**2
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
for j_, m2_ in enumerate([1.2, 1.25, 1.3, 1.35, 1.45, 1.5, 1.55]):
    PDFS[3][f"{m2_:.2f}"] = f"MSHT20nnlo_mcrange_nf3/{j_+1}"
for j_, m2_ in enumerate([4, 4.25, 4.5, 5, 5.25, 5.5]):
    PDFS[4][f"{m2_:.2f}"] = f"MSHT20nnlo_mbrange_nf4/{j_+1}"


def subgrid_path(m2: float, nf: int, j: int) -> str:
    """Return grid path for single configuration."""
    return f"subgrids/{LABELS[nf]}-{m2:.2f}-{j}.pineappl.lz4"


def grid_path(m2: float, nf: int) -> str:
    """Return combined grid path."""
    return f"{LABELS[nf]}-{m2:.2f}.pineappl.lz4"


def compute(ndata: int, m2: float, nl: int, pdf: str, processes: int = -1) -> None:
    """Compute grids."""
    if processes <= 0:
        processes = max(os.cpu_count() + processes, 1)
    start = time.perf_counter()
    print(f"Computing with m2={m2}, nl={nl}, pdf={pdf} using {processes} threads")
    with Pool(processes) as p:
        p.starmap(
            compute_subgrid,
            zip(
                range(ndata),
                [ndata] * ndata,
                [m2] * ndata,
                [nl] * ndata,
                np.geomspace(SH_MIN, SH_MAX, ndata),
                [pdf] * ndata,
            ),
        )
    delta = time.perf_counter() - start
    print("---")
    print(f"computed {ndata} grids in {delta/60:.2f} min")


def compute_subgrid(
    j: int, ndata: int, m2: float, nl: int, sh: float, pdf: str
) -> None:
    """Compute a subgrid."""
    print(f"j = {j:d}/{ndata:d}, S = {sh:e} GeV^2")
    lhapdf.setVerbosity(0)
    start = time.perf_counter()
    # init object
    mk = MaunaKea.MaunaKea(m2, nl, MaunaKea.ORDER_ALL, MaunaKea.LUMI_ALL)
    mk.intCfg.calls = 50000
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
    mk.write(subgrid_path(m2, nl, j))


def merge(ndata: int, m2: float, nl: int) -> None:
    """Merge grids."""
    # merge grids according to their c.o.m. energy
    base = None
    for j in range(ndata):
        subgrid_path_ = pathlib.Path(subgrid_path(m2, nl, j))
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
    merged_path = pathlib.Path(grid_path(m2, nl))
    print(f"Write combined grid to {merged_path}")
    base.write_lz4(merged_path)


def cli() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("m2", help="Mass of heavy quark")
    parser.add_argument("nl", help="Number of light flavors")
    parser.add_argument("ndata", help="Number of points")
    parser.add_argument("-c", "--compute", help="Compute grids", action="store_true")
    parser.add_argument("-m", "--merge", help="Merge grids", action="store_true")
    parser.add_argument("--pdf", help="PDF set used for computing")
    parser.add_argument(
        "--processes", default=-1, help="Number of parallel threads for computing"
    )
    # prepare args
    args = parser.parse_args()
    m2_: float = float(args.m2)
    nl_: int = int(args.nl)
    # determine PDF
    pdf_ = None
    if args.pdf:
        pdf_ = args.pdf
    else:
        pdf_ = PDFS[nl_].get(f"{m2_:.2f}")
    if pdf_ is None or len(pdf_.strip()) <= 0:
        raise ValueError("No PDF set given!")
    # do something
    if args.compute:
        compute(int(args.ndata), m2_, nl_, pdf_, int(args.processes))
    if args.merge:
        merge(int(args.ndata), m2_, nl_)


if __name__ == "__main__":
    cli()
