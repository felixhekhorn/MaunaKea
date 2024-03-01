import argparse
import pathlib
import subprocess

import lhapdf
import pandas as pd
import pineappl
from jinja2 import Environment, FileSystemLoader, select_autoescape

import MaunaKea

GRID_PATH = pathlib.Path("Hathor.pineappl.lz4")


def compute_mauna_kea(nl: int, m: float, sqrt_s: int, pdf: str) -> None:
    """Compute MaunaKea numbers."""
    # init object
    mk = MaunaKea.MaunaKea(m * m, nl, MaunaKea.ORDER_ALL, MaunaKea.LUMI_ALL)
    mk.intCfg.calls = 50000
    mk.setHadronicS(sqrt_s * sqrt_s)
    mk.setPDF(pdf, 0)
    # fill the grid
    mk.run()
    int_out = mk.getIntegrationOutput()
    print(f"sigma_tot = {int_out.result:e} +- {int_out.error:e} [pb]\n")
    # save
    mk.write(str(GRID_PATH))


def compute_hathor(nl: int, m: float, sqrt_s: int, pdf: str) -> None:
    """Compute Hathor numbers."""
    # prepare runcard
    env = Environment(loader=FileSystemLoader("."), autoescape=select_autoescape())
    template = env.get_template("run.cpp.template")
    p_in = pathlib.Path("run.cpp")
    my_vars = {
        "nl": nl,
        "m": m,
        "sqrt_s": sqrt_s,
        "pdf": pdf,
        "scheme": 0,
    }
    labs = {
        "gg": "GG",
        "qqbar": "QQB",
        "qg": "QG",
        "qq": "QQ",
        "qqbarprime": "QQPB",
        "qqprime": "QQP",
    }
    data = {}
    for ch in labs.keys():
        res = []
        for pto in range(2 + 1):
            v = my_vars.copy()
            # combine the "scheme" together
            my_scheme = labs.copy()
            del my_scheme[ch]
            my_scheme = ["NO" + e for e in my_scheme.values()]
            my_scheme.append("N" * pto + "LO")
            v["scheme"] = " | ".join(["Hathor::" + e for e in my_scheme])
            p_in.write_text(template.render(**v), encoding="utf8")
            # compile
            subprocess.run(["make"], check=True)
            # run
            result = subprocess.run("./run.exe", check=True, capture_output=True)
            # read
            res.append(float(result.stdout.splitlines()[-3].split()[-2]))
        data[ch] = res
    df = pd.DataFrame.from_records(data)
    df.to_csv("Hathor.csv")


def compare(pdf: str) -> None:
    """Compare numbers."""
    # MaunaKea
    grid_path = pathlib.Path(GRID_PATH)
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF(pdf, 0)
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
    # Hathor
    hathor = pd.read_csv("Hathor.csv", index_col=0)
    print("Hathor\n", "-" * 5)
    print(hathor)
    mk = pd.DataFrame.from_records(data)
    print()
    print("MaunaKea\n", "-" * 8)
    print(mk)
    diff = (mk - hathor) / hathor
    print()
    print("rel. diff\n", "-" * 8)
    print(diff)


def main() -> None:
    """CLI entry point"""
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", help="MaunaKea|Hathor|compare")
    args = parser.parse_args()

    nl: float = 3
    m: float = 1.51
    sqrt_s: float = 7e3
    pdf = "NNPDF40_nlo_pch_as_01180_nf_3"
    mode = args.mode.strip().lower()
    if mode == "maunakea":
        compute_mauna_kea(nl, m, sqrt_s, pdf)
    elif mode == "hathor":
        compute_hathor(nl, m, sqrt_s, pdf)
    elif mode == "compare":
        compare(pdf)
    else:
        raise ValueError("unkown mode")


if __name__ == "__main__":
    main()
