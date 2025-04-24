"""Compute correction factors."""

import argparse
import pathlib
from collections.abc import Sequence

import lhapdf
import numpy as np
import pandas as pd
import pineappl

from run import DATA, LABELS, PDFS, A, grid_path


def change_collision(
    grid: pineappl.grid.Grid, id_a: int, pdf_a: lhapdf.PDF, id_b: int, pdf_b: lhapdf.PDF
) -> Sequence[float]:
    """Change collision system.

    We can reuse the same grid also for p-X collision if we adjust the PDG ID.
    We need to average over the two possibilities, as e.g. the LO (q-qbar) channel under
    charge conjungation of one hadron will become (ubar-ubar), (u,u), ...,
    but we only save a single combination (u,ubar), ... in the grid.
    """
    return (
        grid.convolute_with_two(
            id_a,
            pdf_a.xfxQ2,
            id_b,
            pdf_b.xfxQ2,
            pdf_a.alphasQ2,
        )
        + grid.convolute_with_two(
            id_b,
            pdf_b.xfxQ2,
            id_a,
            pdf_a.xfxQ2,
            pdf_a.alphasQ2,
        )
    ) / 2.0


def run(mass: float, nl: int, pdf: str) -> None:
    """Compute corrections."""
    grid_path_ = pathlib.Path(grid_path(mass, nl))
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF(pdf)

    grid = pineappl.grid.Grid.read(grid_path_)
    pp = grid.convolute_with_one(
        2212,
        central_pdf.xfxQ2,
        central_pdf.alphasQ2,
    )
    ppbar = change_collision(grid, 2212, central_pdf, -2212, central_pdf)
    nuclear = {
        A["Au"]: change_collision(
            grid, 2212, central_pdf, 2212, lhapdf.mkPDF("EPPS21nlo_CT18Anlo_Au197")
        ),
        A["Pb"]: change_collision(
            grid, 2212, central_pdf, 2212, lhapdf.mkPDF("EPPS21nlo_CT18Anlo_Pb208")
        ),
        # no EPPS PDFs for S, we use Calcium40 instead
        A["S"]: change_collision(
            grid, 2212, central_pdf, 2212, lhapdf.mkPDF("EPPS21nlo_CT18Anlo_Ca40")
        ),
    }
    sqrt_s = grid.bin_left(0)
    # iterate experiments
    data = []
    for cfg in DATA[abs(nl)]:
        # find current energy
        idx = np.argsort(np.abs(sqrt_s - cfg.sqrts))[0]
        el = {
            "no": cfg.no,
            "sqrt_s": cfg.sqrts,
            "pp": pp[idx],
            "pX": pp[idx],
            "X": cfg.correction,
        }
        # apply correction
        if cfg.correction is None:
            pass
        elif cfg.correction == -1:
            el["pX"] = ppbar[idx]
        else:
            el["pX"] = nuclear[cfg.correction][idx]
        data.append(el)
    df = pd.DataFrame.from_records(data)
    df["R"] = df["pp"] / df["pX"]
    df.to_csv(f"data/{LABELS[nl]}-{mass:.2f}-{pdf}-corrections.csv")


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("mass", type=float, help="Mass of heavy quark")
    parser.add_argument("nl", type=int, help="Number of light flavors")
    parser.add_argument("--pdf", help="PDF set used for computing pp")
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
    run(mass_, nl_, pdf_)


if __name__ == "__main__":
    main()
