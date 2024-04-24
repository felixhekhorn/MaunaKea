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
MASSES = {3: 1.51**2, 4: 4.92**2}
SH_MIN: float = 20.0**2
SH_MAX: float = 400e3**2
PDFS = {3: "NNPDF40_nnlo_pch_as_01180_nf_3", 4: "NNPDF40_nnlo_as_01180_nf_4"}


def sub_grid_path(nf: int, j: int) -> str:
    """Return grid path for single configuration."""
    return f"subgrids/{LABELS[nf]}-{j}.pineappl.lz4"


def grid_path(nf: int) -> str:
    """Return combined grid path."""
    return f"{LABELS[nf]}.pineappl.lz4"


def compute(nl: int, ndata: int, processes: int = -1) -> None:
    """Compute grids."""
    if processes <= 0:
        processes = max(os.cpu_count() + processes, 1)
    start = time.perf_counter()
    with Pool(processes) as p:
        p.starmap(
            compute_subgrid,
            zip(
                range(ndata),
                [ndata] * ndata,
                [nl] * ndata,
                np.geomspace(SH_MIN, SH_MAX, ndata),
            ),
        )
    delta = time.perf_counter() - start
    print("---")
    print(f"computed {ndata} grids in {delta/60:.2f} min")


def compute_subgrid(j: int, ndata: int, nl: int, sh: float) -> None:
    """Compute a subgrid."""
    m2: float = MASSES[nl]
    print(f"j = {j:d}/{ndata:d}, S = {sh:e} GeV^2")
    lhapdf.setVerbosity(0)
    start = time.perf_counter()
    # init object
    mk = MaunaKea.MaunaKea(m2, nl, MaunaKea.ORDER_ALL, MaunaKea.LUMI_ALL)
    mk.intCfg.calls = 5000
    mk.setHadronicS(sh)
    mk.setPDF(PDFS[nl], 0)
    mk.setCentralScaleRatio(2.0)
    # fill the grid
    mk.run()
    int_out = mk.getIntegrationOutput()
    delta = time.perf_counter() - start
    print(f"sigma_tot = {int_out.result:e} +- {int_out.error:e} [pb]")
    print(f"took {delta:.2f} s and chi2iter={int_out.chi2iter}")
    # save
    mk.write(sub_grid_path(nl, j))


def merge(nl: int, ndata: int) -> None:
    """Merge grids."""
    # merge grids according to their c.o.m. energy
    base = None
    for j in range(ndata):
        sub_grid_path_ = pathlib.Path(sub_grid_path(nl, j))
        grid = pineappl.grid.Grid.read(sub_grid_path_)
        # rebin by S_h
        sqrt_sh = np.sqrt(float(grid.key_values()["MaunaKea/hadronicS"]))
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
    merged_path = pathlib.Path(grid_path(nl))
    print(f"Write combined grid to {merged_path}")
    base.write_lz4(merged_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("nl", help="Number of light flavors")
    parser.add_argument("ndata", help="Number of points")
    parser.add_argument("-c", "--compute", help="Compute grids", action="store_true")
    parser.add_argument("-m", "--merge", help="Merge grids", action="store_true")
    parser.add_argument("--processes", default=-1, help="Number of parallel threads")
    args = parser.parse_args()
    if args.compute:
        compute(int(args.nl), int(args.ndata), int(args.processes))
    if args.merge:
        merge(int(args.nl), int(args.ndata))
