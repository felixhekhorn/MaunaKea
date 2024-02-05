import MaunaKea


def toppp():
    nl: int = 5
    m2: float = 172.5**2
    S_h: float = 7e3**2
    # init object
    mk = MaunaKea.MaunaKea(m2, nl, 1, 1)
    # mk.intCfg.calls = 5000;
    mk.intCfg.verbosity = 3
    mk.setHadronicS(S_h)
    mk.setPDF("NNPDF40_nnlo_as_01180", 0)
    # mk.setGridCentralScaleRatio(2.);
    # fill the grid
    mk.run()
    intOut = mk.getIntegrationOutput()
    print("sigma_tot = {:e} +- {:e} [pb]\n".format(intOut.result, intOut.error))
    # save
    mk.write("MaunaKea.pineappl.lz4")


if __name__ == "__main__":
    toppp()
