"""Top cross sections measured at LHC experiments."""

import MaunaKea


def run() -> None:
    """Compute all grids."""
    for m in [170.0, 172.5, 175.0]:  # [GeV^2] mass of heavy quark
        tagM = f"{m}".replace(".", "p")
        for sqrtS_h in [5e3, 7e3, 8e3, 13e3]:  # [GeV] collider energy
            tagS = f"{sqrtS_h/1e3}".replace(".", "p")
            grid(m * m, sqrtS_h * sqrtS_h, f"sigma_tot-{tagS}-{tagM}.pineappl.lz4")


def grid(m2: float, S_h: float, path: str) -> None:
    """Compute a single grid."""
    nl: int = 5  # number of light flavors
    # init object
    mk = MaunaKea.MaunaKea(m2, nl, MaunaKea.ORDER_ALL, MaunaKea.LUMI_ALL)
    # configure
    # mk.intCfg.calls = 5000
    mk.intCfg.verbosity = 3
    mk.setHadronicS(S_h)
    mk.setPDF("NNPDF40_nnlo_as_01180", 0)
    # mk.setCentralScaleRatio(2.0)
    # fill the grid
    mk.run()
    output = mk.getIntegrationOutput()
    print(f"sigma_tot = {output.result:e} +- {output.error:e} [pb]")
    # save
    print(f"Save grid to {path}")
    mk.write(path)


if __name__ == "__main__":
    run()
