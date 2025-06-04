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
    grid: pineappl.grid.Grid,
    id_a: int,
    pdf_a: lhapdf.PDF,
    id_b: int,
    pdf_b: lhapdf.PDF,
    order_mask=np.array([], dtype=bool),
) -> Sequence[float]:
    """Change collision system.

    We can reuse the same grid also for p-X collision if we adjust the PDG ID.
    We need to average over the two possibilities, as e.g. the LO (q-qbar) channel under
    charge conjungation of one hadron will become (ubar-ubar), (u,u), ...,
    but we only save a single combination (u,ubar), ... in the grid.
    """
    return (
        grid.convolute_with_two(
            id_a, pdf_a.xfxQ2, id_b, pdf_b.xfxQ2, pdf_a.alphasQ2, order_mask=order_mask
        )
        + grid.convolute_with_two(
            id_b, pdf_b.xfxQ2, id_a, pdf_a.xfxQ2, pdf_a.alphasQ2, order_mask=order_mask
        )
    ) / 2.0


def run(mass: float, nl: int, pdf: str, use_nnlo_grids: bool = False) -> None:
    """Compute corrections."""
    grid_path_ = pathlib.Path(grid_path(mass, nl))
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF(pdf)
    grid = pineappl.grid.Grid.read(grid_path_)
    order_mask = pineappl.grid.Order.create_mask(
        grid.orders(), 2 - 1 + (2 if use_nnlo_grids else 1), 0, True
    )
    pp = grid.convolute_with_one(
        2212, central_pdf.xfxQ2, central_pdf.alphasQ2, order_mask=order_mask
    )
    ppbar = change_collision(
        grid, 2212, central_pdf, -2212, central_pdf, order_mask=order_mask
    )
    nuclear = {
        A["Au"]: change_collision(
            grid,
            2212,
            central_pdf,
            2212,
            lhapdf.mkPDF("EPPS21nlo_CT18Anlo_Au197"),
            order_mask=order_mask,
        ),
        A["Pb"]: change_collision(
            grid,
            2212,
            central_pdf,
            2212,
            lhapdf.mkPDF("EPPS21nlo_CT18Anlo_Pb208"),
            order_mask=order_mask,
        ),
        # no EPPS PDFs for S, we use Calcium40 instead
        A["S"]: change_collision(
            grid,
            2212,
            central_pdf,
            2212,
            lhapdf.mkPDF("EPPS21nlo_CT18Anlo_Ca40"),
            order_mask=order_mask,
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
    nnlo_tag = "-nnlo" if use_nnlo_grids else ""
    df.to_csv(f"data/{LABELS[nl]}-{mass:.2f}-{pdf}-corrections{nnlo_tag}.csv")


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("mass", type=float, help="Mass of heavy quark")
    parser.add_argument("nl", type=int, help="Number of light flavors")
    parser.add_argument("--pdf", help="PDF set used for computing pp")
    parser.add_argument("--nnlo-grids", help="Use grids at NNLO", action="store_true")
    args = parser.parse_args()
    mass_: float = float(args.mass)
    nl_: int = int(args.nl)
    use_nnlo_grids: bool = bool(args.nnlo_grids)
    # determine PDF
    pdf_ = None
    if args.pdf:
        pdf_ = args.pdf
    else:
        pdf_ = PDFS[abs(nl_)].get(f"{mass_:.2f}")
    if pdf_ is None or len(pdf_.strip()) <= 0:
        raise ValueError("No PDF set given!")
    run(mass_, nl_, pdf_, use_nnlo_grids)


if __name__ == "__main__":
    main()
