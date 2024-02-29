import pathlib
import subprocess

import lhapdf
import pineappl
from jinja2 import Environment, FileSystemLoader, select_autoescape

import MaunaKea

GRID_PATH = pathlib.Path("toppp.pineappl.lz4")
TOPPP_EXE = pathlib.Path(__file__).parents[3] / "top++2.0" / "top++"


def compute_mauna_kea(m: float, sqrt_s: int) -> None:
    """Compute MaunaKea numbers."""
    # init object
    mk = MaunaKea.MaunaKea(m * m, 5, MaunaKea.ORDER_LO, MaunaKea.LUMI_ALL)
    # mk.intCfg.calls = 50000
    mk.setHadronicS(sqrt_s * sqrt_s)
    mk.setPDF(f"NNPDF40_nnlo_as_01180", 0)
    # fill the grid
    mk.run()
    int_out = mk.getIntegrationOutput()
    print(f"sigma_tot = {int_out.result:e} +- {int_out.error:e} [pb]\n")
    # save
    mk.write(GRID_PATH)


def compute_toppp(m: float, sqrt_s: int) -> None:
    """Compute top++ numbers."""
    # prepare runcard
    env = Environment(loader=FileSystemLoader("."), autoescape=select_autoescape())
    template = env.get_template("top++.cfg.template")
    p = pathlib.Path("top++.cfg")
    p.write_text(template.render(m=m, sqrt_s=sqrt_s))
    # run
    subprocess.run([TOPPP_EXE, p])


def compare() -> None:
    # MaunaKea
    grid_path = pathlib.Path(GRID_PATH)
    lhapdf.setVerbosity(0)
    central_pdf = lhapdf.mkPDF("NNPDF40_nnlo_as_01180", 0)
    grid = pineappl.grid.Grid.read(grid_path)
    central = grid.convolute_with_one(2212, central_pdf.xfxQ2, central_pdf.alphasQ2)[0]
    # top++
    p = pathlib.Path("top++.res")
    out = p.read_text().splitlines()
    toppp_res = float(out[14].split()[-1])
    print(f"MaunaKea = {central} pb")
    print(f"top+++ = {toppp_res} pb")
    diff = (central - toppp_res) / toppp_res
    print(f"diff = {diff}")


def main() -> None:
    m: float = 172.5
    sqrt_s: float = 7e3
    # compute_mauna_kea(m,sqrt_s)
    # compute_toppp(m,sqrt_s)
    compare()


if __name__ == "__main__":
    main()
