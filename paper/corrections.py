"""Compute correction factors."""

import argparse
import pathlib

import lhapdf
import pandas as pd
import pineappl

from run import LABELS, grid_path

CORRECTIONS = {3: [None, "pbar", None, None, None]}


def run(m2: float, nl: int) -> None:
    grid_path_ = pathlib.Path(grid_path(m2, nl))
    pdf = "NNPDF40_nnlo_pch_as_01180_nf_3"
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF(pdf)
    grid = pineappl.grid.Grid.read(grid_path_)
    pp = grid.convolute_with_one(
        2212,
        central_pdf.xfxQ2,
        central_pdf.alphasQ2,
    )
    # We can reuse the same grid also for p-pbar collision if we adjust the PDG ID.
    # We need to average over the two possibilities, as e.g. the LO (q-qbar) channel under
    # charge conjungation of one hadron will become (ubar-ubar), (u,u), ...,
    # but we only save a single combination (u,ubar), ... in the grid.
    ppbar = (
        grid.convolute_with_two(
            2212,
            central_pdf.xfxQ2,
            -2212,
            central_pdf.xfxQ2,
            central_pdf.alphasQ2,
        )
        + grid.convolute_with_two(
            -2212,
            central_pdf.xfxQ2,
            2212,
            central_pdf.xfxQ2,
            central_pdf.alphasQ2,
        )
    ) / 2.0
    df = pd.DataFrame()
    df["sqrt_s"] = grid.bin_left(0)
    df["pp"] = pp
    df["ppbar"] = 0.0
    df["R"] = 1.0
    for j, mod in enumerate(CORRECTIONS[abs(nl)]):
        if mod is None:
            continue
        df.loc[j, "ppbar"] = ppbar[j]
        df.loc[j, "R"] = pp[j] / ppbar[j]
    df.to_csv(f"data/{LABELS[nl]}-{m2:.2f}-{pdf}-corrections.csv")


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("m2", type=float, help="Mass of heavy quark")
    parser.add_argument("nl", type=int, help="Number of light flavors")
    args = parser.parse_args()
    m2_: float = float(args.m2)
    nl_: int = int(args.nl)
    run(m2_, nl_)


if __name__ == "__main__":
    main()
