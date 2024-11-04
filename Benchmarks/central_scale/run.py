"""Check `setCentralScaleRatio`."""

import lhapdf
import numpy as np
import pandas as pd
import pineappl

import MaunaKea

XI: float = 2.0
MOD_GRID = "up.pineappl.lz4"
CENTRAL_GRID = "central.pineappl.lz4"


def run(ratio: float, grid: str) -> None:
    """Compute a given grid."""
    nl: int = 3  # number of light flavors
    m2: float = 1.51**2  # [GeV^2] mass of heavy quark
    S_h: float = 8e3**2  # [GeV^2] collider energy
    # init object
    mk = MaunaKea.MaunaKea(m2, nl, MaunaKea.ORDER_LO, MaunaKea.LUMI_ALL)
    # configure
    # mk.intCfg.calls = 5000
    mk.intCfg.verbosity = 3
    mk.setHadronicS(S_h)
    mk.setPDF("NNPDF40_nnlo_as_01180", 0)
    mk.setCentralScaleRatio(ratio)
    # fill the grid
    mk.run()
    output = mk.getIntegrationOutput()
    print(f"sigma_tot = {output.result:e} +- {output.error:e} [pb]")
    # save
    mk.write(grid)


def build_grids() -> None:
    """Generate the two candidates."""
    run(XI, MOD_GRID)
    run(1.0, CENTRAL_GRID)


def compare() -> None:
    """Compare the two grids."""
    mod = pineappl.grid.Grid.read(MOD_GRID)
    central = pineappl.grid.Grid.read(CENTRAL_GRID)
    pdf = lhapdf.mkPDF("NNPDF40_nnlo_as_01180", 0)
    xis = np.geomspace(1.0, XI * XI, 5)
    res_central = central.convolute_with_one(
        2212, pdf.xfxQ2, pdf.alphasQ2, xi=[(xi, xi) for xi in xis]
    )
    res_mod = mod.convolute_with_one(
        2212, pdf.xfxQ2, pdf.alphasQ2, xi=[(xi / XI, xi / XI) for xi in xis]
    )
    df = pd.DataFrame({"xi": xis, "central": res_central, "mod": res_mod})
    df["rel."] = (df["mod"] - df["central"]) / df["central"]
    print(df)


if __name__ == "__main__":
    build_grids()
    compare()
