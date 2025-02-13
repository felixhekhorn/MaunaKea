# To run this script:
# 1. you may want to create a virtual environment, e.g. `python3 -m venv env`, and activate it, e.g. `. env/bin activate`
# 2. install pineappl python interface: `pip install pineappl`
import lhapdf
import numpy as np
import pineappl

# setup LHAPDF and load PDFs
lhapdf.setVerbosity(0)
p_pdf = lhapdf.mkPDF("CT18ANLO", 0)
nuclear_pdf = lhapdf.mkPDF("EPPS21nlo_CT18Anlo_Pb208", 0)
# load grid
grid = pineappl.grid.Grid.read("ccbar-1.30.pineappl.lz4")
# bins are in sqrt(s)
sqrts = grid.bin_left(0)
# restrict to NLO
order_mask = pineappl.grid.Order.create_mask(grid.orders(), 2, 0, True)
# choose scale
xi = 1.0
# do it! we just need to undo one optimization, since the grid was generated for p-p:
# We can reuse the same grid also for p-X collision if we average over the two input state.
# E.g. the LO (q-qbar) channel is just that and we need the q in coming from both sides.
pPb = (
    grid.convolve_with_two(
        2212,
        p_pdf.xfxQ2,
        2212,
        nuclear_pdf.xfxQ2,
        p_pdf.alphasQ2,
        order_mask=order_mask,
        xi=[(xi, xi)],
    )
    + grid.convolve_with_two(
        2212,
        nuclear_pdf.xfxQ2,
        2212,
        p_pdf.xfxQ2,
        p_pdf.alphasQ2,
        order_mask=order_mask,
        xi=[(xi, xi)],
    )
) / 2.0
# output
data = np.array([sqrts, pPb])
np.savetxt(f"Hannu-xi_{xi}.txt", data.T)
print(data.T)
