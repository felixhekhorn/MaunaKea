import argparse
import pathlib
import subprocess

import lhapdf
import pandas as pd
import pineappl
from jinja2 import Environment, FileSystemLoader, select_autoescape

import MaunaKea

TOPPP_NF5_EXE = pathlib.Path(__file__).parents[3] / "top++2.0" / "top++"
TOPPP_NF3_EXE = pathlib.Path(__file__).parents[3] / "top++2.0-c" / "top++"


def grid_path(nl: int):
    """MaunaKea output file name."""
    return pathlib.Path(f"toppp-{nl}.pineappl.lz4")


def compute_mauna_kea(nl: int, m: float, sqrt_s: int, pdf: str) -> None:
    """Compute MaunaKea numbers."""
    # init object
    mk = MaunaKea.MaunaKea(m * m, nl, MaunaKea.ORDER_ALL, MaunaKea.LUMI_ALL)
    mk.intCfg.calls = 50000
    mk.intCfg.verbosity = 3
    mk.setHadronicS(sqrt_s * sqrt_s)
    mk.setPDF(pdf, 0)
    # fill the grid
    mk.run()
    int_out = mk.getIntegrationOutput()
    print(f"sigma_tot = {int_out.result:e} +- {int_out.error:e} [pb]\n")
    # save
    mk.write(str(grid_path(nl)))


TOPPP_CFG_TEMPLATE = Environment(
    loader=FileSystemLoader("."), autoescape=select_autoescape()
).get_template("top++.cfg.template")
TOPPP_DEFAULT_VARS = {
    "m": 175.0,
    "sqrt_s": 7e3,
    "svrange": 1.0,
    "murs": (1.0,),
    "mufs": (1.0,),
    "use_lo": False,
    "use_nlo": False,
    "use_nnlo": False,
    "pc": "ALL",
    "pdf": "NNPDF40_nnlo_as_01180",
}


def run_toppp(
    n_output: int,
    nl: int,
    m: float,
    sqrt_s: int,
    pto: int,
    pdf: str,
    ch: str = "ALL",
    murs=(1.0,),
    mufs=(1.0,),
) -> None:
    """Run top++."""
    # prepare runcard
    p_in = pathlib.Path("top++.cfg")
    v = TOPPP_DEFAULT_VARS.copy()
    v["m"] = m
    v["sqrt_s"] = sqrt_s
    v["svrange"] = max(mufs) * max(murs)
    v["murs"] = murs
    v["mufs"] = mufs
    v["use_" + ("n" * pto) + "lo"] = True
    v["pc"] = ch
    v["pdf"] = pdf
    p_in.write_text(TOPPP_CFG_TEMPLATE.render(**v), encoding="utf8")
    # run
    if nl == 3:
        exe = TOPPP_NF3_EXE
    elif nl == 5:
        exe = TOPPP_NF5_EXE
    else:
        raise ValueError(f"No executable for nl={nl}")
    subprocess.run([exe, p_in], check=True)
    # read file
    p_out = pathlib.Path("top++.res")
    out = p_out.read_text(encoding="utf8").splitlines()
    # parse numbers
    res = []
    for j in range(n_output):
        res.append(float(out[15 + j].split()[-1]))
    return res


def compute_toppp(nl: int, m: float, sqrt_s: int, pdf: str) -> None:
    """Compute top++ numbers."""
    labs = ["gg", "qqbar", "qg", "qq", "qqbarprime", "qqprime"]
    data = {}
    for ch in labs:
        res = []
        for pto in range(2 + 1):
            res.append(run_toppp(1, nl, m, sqrt_s, pto, pdf, ch)[0])
        data[ch] = res
    df = pd.DataFrame.from_records(data)
    df.to_csv("top++.csv")


def compute_toppp_sv(
    nl: int, m: float, sqrt_s: int, pto: int, pdf: str, ch: str, abs_xi: float
) -> None:
    """Compute SV top++ numbers."""
    xis = [1.0 / abs_xi, 1.0, abs_xi]
    data = []
    res = list(reversed(run_toppp(9, nl, m, sqrt_s, pto, pdf, ch, murs=xis, mufs=xis)))
    for muf in xis:
        for mur in xis:
            data.append(
                {"pto": pto, "ch": ch, "muf": muf, "mur": mur, "top++": res.pop()}
            )
    df = pd.DataFrame.from_records(data)
    df.to_csv("top++_sv.csv")


def compare(nl: int, pdf: str) -> None:
    """Compare numbers."""
    # MaunaKea
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF(pdf, 0)
    grid = pineappl.grid.Grid.read(grid_path(nl))
    data = {}
    labs = ["gg", "qqbar", "qg", "qq", "qqbarprime", "qqprime"]
    for idx, ch in enumerate(labs):
        lm = [False] * len(labs)
        lm[idx] = True
        me = []
        for pto in range(2 + 1):
            # extract correction
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
    print("top++\n", "-" * 5, sep="")
    print(toppp)
    mk = pd.DataFrame.from_records(data)
    print()
    print("MaunaKea\n", "-" * 8, sep="")
    print(mk)
    diff = (mk - toppp) / toppp
    print()
    print("rel. diff\n", "-" * 8, sep="")
    print(diff)


def compare_sv(nl: int, pto: int, pdf: str, ch: str, abs_xi: float) -> None:
    """Compare SV numbers."""
    # MaunaKea
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF(pdf, 0)
    grid = pineappl.grid.Grid.read(grid_path(nl))
    labs = ["gg", "qqbar", "qg", "qq", "qqbarprime", "qqprime"]
    lm = [False] * len(labs)
    lm[labs.index(ch)] = True
    me = []
    # extract correction
    om_true = pineappl.grid.Order.create_mask(grid.orders(), pto + 1, 0, True)
    om_lower = pineappl.grid.Order.create_mask(grid.orders(), pto, 0, True)
    om = [a ^ b for (a, b) in zip(om_true, om_lower)]
    lin_xis = [1.0 / abs_xi, 1.0, abs_xi]
    xis = []
    for muf in lin_xis:
        for mur in lin_xis:
            xis.append((mur, muf))
    me = grid.convolute_with_one(
        2212,
        central_pdf.xfxQ2,
        central_pdf.alphasQ2,
        order_mask=om,
        lumi_mask=lm,
        xi=xis,
    )
    # load top++
    df = pd.read_csv("top++_sv.csv", index_col=0)
    df["MaunaKea"] = me
    df["rel. error"] = (df["MaunaKea"] - df["top++"]) / df["top++"]
    print(df)


def main() -> None:
    """CLI entry point"""
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", help="MaunaKea|top++|compare|top++-sv|compare-sv")
    parser.add_argument("nl", help="number of light flavors")
    args = parser.parse_args()

    nl: int = int(args.nl)
    m: float = 172.5
    pdf = "NNPDF40_nnlo_as_01180"
    if nl == 3:
        m = 1.51
        pdf = "NNPDF40_nlo_pch_as_01180_nf_3"
    sqrt_s: float = 7e3
    mode: str = args.mode.strip().lower()
    pto: int = 2
    ch: str = "qg"
    abs_xi: float = 2.0

    if mode == "maunakea":
        compute_mauna_kea(nl, m, sqrt_s, pdf)
    elif mode == "top++":
        compute_toppp(nl, m, sqrt_s, pdf)
    elif mode == "top++-sv":
        compute_toppp_sv(nl, m, sqrt_s, pto, pdf, ch, abs_xi)
    elif mode == "compare":
        compare(nl, pdf)
    elif mode == "compare-sv":
        compare_sv(nl, pto, pdf, ch, abs_xi)
    else:
        raise ValueError("unkown mode")


if __name__ == "__main__":
    main()
