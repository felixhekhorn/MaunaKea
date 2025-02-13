import MaunaKea


def run() -> None:
    nl: int = 5  # number of light flavors
    m2: float = 173.2**2  # [GeV^2] mass of heavy quark
    S_h: float = 8e3**2  # [GeV^2] collider energy
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
    mk.write("top-8-173p2.pineappl.lz4")


if __name__ == "__main__":
    run()
