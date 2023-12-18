import dataclasses
import pathlib

import lhapdf
import numpy as np
import pineappl


@dataclasses.dataclass()
class Ref:
    m_top: float = 172.5
    xiF: float = 1.0
    xiR: float = 1.0
    toppp: float = np.NaN
    MaunaKea: float = np.NaN


here = pathlib.Path(__file__).parent

# read top++
topout_path = here.parent / "top++2.0" / "out.txt"


def parse() -> list[Ref]:
    topout = topout_path.read_text(encoding="utf8")
    start_marker = "Computing the cross-section:"
    start = topout.find(start_marker)
    end = topout.find("......... Cross-section computed!")
    raw_lines = topout[start + len(start_marker) : end].strip().split("\n")
    lines = []
    for l in raw_lines:
        raw_els = l.split(",")
        els = []
        for raw_el in raw_els:
            el = raw_el.split("=")
            els.append(el[-1])
        els[-1] = els[-1][: -len(" [pb].")]
        lines.append(Ref(*[float(e) for e in els]))
    return lines


# put me on top
g_path = here / "MaunaKea.pineappl.lz4"
g = pineappl.grid.Grid.read(g_path)
pdf = lhapdf.mkPDF("NNPDF40_nnlo_as_01180", 0)
els = parse()
for el in els:
    el.MaunaKea = g.convolute_with_one(
        2212, pdf.xfxQ2, pdf.alphasQ2, xi=[(el.xiR, el.xiF)]
    )[0]

# compute diff
diff = [(el.toppp - el.MaunaKea) / el.toppp for el in els]


# Thanks https://stackoverflow.com/questions/33964913/equivalent-of-polyfit-for-a-2d-polynomial-in-python
def polyfit2d(x, y, z, kx=3, ky=3, order=None):
    """
    Two dimensional polynomial fitting by least squares.
    Fits the functional form f(x,y) = z.

    Notes
    -----
    Resultant fit can be plotted with:
    np.polynomial.polynomial.polygrid2d(x, y, soln.reshape((kx+1, ky+1)))

    Parameters
    ----------
    x, y: array-like, 1d
        x and y coordinates.
    z: np.ndarray, 2d
        Surface to fit.
    kx, ky: int, default is 3
        Polynomial order in x and y, respectively.
    order: int or None, default is None
        If None, all coefficients up to maxiumum kx, ky, ie. up to and including x^kx*y^ky, are considered.
        If int, coefficients up to a maximum of kx+ky <= order are considered.

    Returns
    -------
    Return paramters from np.linalg.lstsq.

    soln: np.ndarray
        Array of polynomial coefficients.
    residuals: np.ndarray
    rank: int
    s: np.ndarray

    """

    # grid coords
    x, y = np.meshgrid(x, y)
    # coefficient array, up to x^kx, y^ky
    coeffs = np.ones((kx + 1, ky + 1))

    # solve array
    a = np.zeros((coeffs.size, x.size))

    # for each coefficient produce array x^i, y^j
    for index, (j, i) in enumerate(np.ndindex(coeffs.shape)):
        # do not include powers greater than order
        if order is not None and i + j > order:
            arr = np.zeros_like(x)
        else:
            arr = coeffs[i, j] * x**i * y**j
        a[index] = arr.ravel()

    # do leastsq fitting and return leastsq result
    return np.linalg.lstsq(a.T, np.ravel(z), rcond=None)


xifs = np.unique(np.array([el.xiF for el in els]))
xirs = np.unique(np.array([el.xiR for el in els]))
np.testing.assert_allclose(xifs, xirs)
xis = 2.0 * np.log(xifs)
lxis = len(xis)

print("rel. diff top++ vs. MaunaKea")
print(np.array(diff).reshape(lxis, lxis))

coeff_tp, _r, _rank, _s = polyfit2d(
    xis, xis, np.array([(el.toppp) for el in els]).reshape(lxis, lxis), 2, 2, 2
)
coeff_mk, _r, _rank, _s = polyfit2d(
    xis, xis, np.array([(el.MaunaKea) for el in els]).reshape(lxis, lxis), 2, 2, 2
)
coeff_tp = coeff_tp.reshape(3, 3)
coeff_mk = coeff_mk.reshape(3, 3)
print("top++")
print(coeff_tp)
print("MaunaKea")
print(coeff_mk)
