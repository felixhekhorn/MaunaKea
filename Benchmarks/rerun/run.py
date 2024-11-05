"""Rerun a existing grid."""

import sys

import pineappl

import MaunaKea


def rerun(grid_path: str, new_path: str) -> None:
    """Recompute a given grid."""
    grid = pineappl.grid.Grid.read(grid_path)
    kv = grid.key_values()
    sep = "/" if "MaunaKea/m2" in kv else "::"
    # init object
    mk = MaunaKea.MaunaKea(
        float(kv[f"MaunaKea{sep}m2"]),
        int(kv[f"MaunaKea{sep}nl"]),
        int(kv[f"MaunaKea{sep}order_mask"]),
        int(kv[f"MaunaKea{sep}lumi_mask"]),
    )
    # configure
    intCfg = kv[f"MaunaKea{sep}IntegrationConfig"]
    for line in intCfg.splitlines():
        var, val = line.rsplit(":", 1)
        mk.intCfg.__setattr__(var, int(val))
    mk.intCfg.verbosity = 3
    mk.setHadronicS(float(kv[f"MaunaKea{sep}hadronicS"]))
    pdf = kv[f"MaunaKea{sep}PDF"]
    if "#" in pdf:
        pdf_set = pdf.rsplit("#", 1)
    elif "/" in pdf:
        pdf_set = pdf.rsplit("/", 1)
    else:
        raise ValueError("Could not read PDF")
    mk.setPDF(pdf_set[0], int(pdf_set[1]))
    mk.setCentralScaleRatio(float(kv[f"MaunaKea{sep}xi"]))
    # fill the grid
    mk.run()
    output = mk.getIntegrationOutput()
    print(f"sigma_tot = {output.result:e} +- {output.error:e} [pb]")
    # save
    mk.write(new_path)


if __name__ == "__main__":
    rerun(sys.argv[1], sys.argv[2])
