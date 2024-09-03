"""Compute correction factors."""

import argparse
import pathlib

import lhapdf
import numpy as np
import pandas as pd
import pineappl

from run import DATA, LABELS, PDFS, grid_path


def run(m2: float, nl: int, pdf: str) -> None:
    grid_path_ = pathlib.Path(grid_path(m2, nl))
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
    sqrt_s = grid.bin_left(0)
    data = []
    for cfg in DATA[abs(nl)]:
        idx = np.argsort(np.abs(sqrt_s - cfg.sqrts))[0]
        el = {"no": cfg.no, "sqrt_s": cfg.sqrts, "pp": pp[idx], "pX": pp[idx]}
        # apply correction
        if cfg.correction == -1:
            el["pX"] = ppbar[idx]
        data.append(el)
    df = pd.DataFrame.from_records(data)
    df["R"] = df["pp"] / df["pX"]
    df.to_csv(f"data/{LABELS[nl]}-{m2:.2f}-{pdf}-corrections.csv")


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("m2", type=float, help="Mass of heavy quark")
    parser.add_argument("nl", type=int, help="Number of light flavors")
    parser.add_argument("--pdf", help="PDF set used for computing pp")
    args = parser.parse_args()
    m2_: float = float(args.m2)
    nl_: int = int(args.nl)
    # determine PDF
    pdf_ = None
    if args.pdf:
        pdf_ = args.pdf
    else:
        pdf_ = PDFS[abs(nl_)].get(f"{m2_:.2f}")
    if pdf_ is None or len(pdf_.strip()) <= 0:
        raise ValueError("No PDF set given!")
    run(m2_, nl_, pdf_)


if __name__ == "__main__":
    main()
