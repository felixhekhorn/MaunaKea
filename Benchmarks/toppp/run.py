import argparse
import pathlib
import subprocess

import lhapdf
import pandas as pd
import pineappl
from jinja2 import Environment, FileSystemLoader, select_autoescape

import MaunaKea

GRID_PATH = pathlib.Path("toppp.pineappl.lz4")
TOPPP_EXE = pathlib.Path(__file__).parents[3] / "top++2.0" / "top++"


def compute_mauna_kea(m: float, sqrt_s: int) -> None:
    """Compute MaunaKea numbers."""
    # init object
    mk = MaunaKea.MaunaKea(m * m, 5, MaunaKea.ORDER_ALL, MaunaKea.LUMI_ALL)
    mk.intCfg.calls = 50000
    mk.setHadronicS(sqrt_s * sqrt_s)
    mk.setPDF("NNPDF40_nnlo_as_01180", 0)
    # fill the grid
    mk.run()
    int_out = mk.getIntegrationOutput()
    print(f"sigma_tot = {int_out.result:e} +- {int_out.error:e} [pb]\n")
    # save
    mk.write(str(GRID_PATH))


def compute_toppp(m: float, sqrt_s: int) -> None:
    """Compute top++ numbers."""
    # prepare runcard
    env = Environment(loader=FileSystemLoader("."), autoescape=select_autoescape())
    template = env.get_template("top++.cfg.template")
    p_in = pathlib.Path("top++.cfg")
    my_vars = {
        "m": m,
        "sqrt_s": sqrt_s,
        "use_lo": False,
        "use_nlo": False,
        "use_nnlo": False,
        "pc": "ALL",
    }
    labs = ["gg", "qqbar", "qg", "qq", "qqbarprime", "qqprime"]
    data = {}
    for ch in labs:
        res = []
        for pto in range(2 + 1):
            v = my_vars.copy()
            v["use_" + ("n" * pto) + "lo"] = True
            v["pc"] = ch
            p_in.write_text(template.render(**v), encoding="utf8")
            # run
            subprocess.run([TOPPP_EXE, p_in], check=True)
            # read
            p_out = pathlib.Path("top++.res")
            out = p_out.read_text(encoding="utf8").splitlines()
            res.append(float(out[15].split()[-1]))
        data[ch] = res
    df = pd.DataFrame.from_records(data)
    df.to_csv("top++.csv")


def compare() -> None:
    """Compare numbers."""
    # MaunaKea
    grid_path = pathlib.Path(GRID_PATH)
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF("NNPDF40_nnlo_as_01180", 0)
    grid = pineappl.grid.Grid.read(grid_path)
    data = {}
    labs = ["gg", "qqbar", "qg", "qq", "qqbarprime", "qqprime"]
    for idx, ch in enumerate(labs):
        lm = [False] * len(labs)
        lm[idx] = True
        me = []
        for pto in range(2 + 1):
            om_true = pineappl.grid.Order.create_mask(grid.orders(), pto + 1, 0, True)
            om_lower = pineappl.grid.Order.create_mask(grid.orders(), pto, 0, True)
            om = [a ^ b for (a, b) in zip(om_true, om_lower)]
            me.append(
                grid.convolute_with_one(
                    2212,
                    central_pdf.xfxQ2,
                    central_pdf.alphasQ2,
                    order_mask=om,
                    lumi_mask=lm,
                )[0]
            )
        data[ch] = me
    # top++
    toppp = pd.read_csv("top++.csv", index_col=0)
    print("top++\n", "-" * 5)
    print(toppp)
    mk = pd.DataFrame.from_records(data)
    print()
    print("MaunaKea\n", "-" * 8)
    print(mk)
    diff = (mk - toppp) / toppp
    print()
    print("rel. diff\n", "-" * 8)
    print(diff)


def main() -> None:
    """CLI entry point"""
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", help="MaunaKea|top++|compare")
    args = parser.parse_args()

    m: float = 172.5
    sqrt_s: float = 7e3
    mode = args.mode.strip().lower()
    if mode == "maunakea":
        compute_mauna_kea(m, sqrt_s)
    elif mode == "top++":
        compute_toppp(m, sqrt_s)
    elif mode == "compare":
        compare()
    else:
        raise ValueError("unkown mode")


if __name__ == "__main__":
    main()
