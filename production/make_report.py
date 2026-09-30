#!/usr/bin/env python3
# Build a summary PDF (one page per plot: image + parameter list) and a ZIP (PNGs + PDF).
# Uses only Pillow.  Run from the project root:  python3 production/make_report.py
import os, glob, zipfile
from PIL import Image, ImageDraw, ImageFont

def find_font(size):
    for p in sorted(glob.glob("/usr/share/fonts/**/*.ttf", recursive=True)):
        try:
            return ImageFont.truetype(p, size)
        except Exception:
            pass
    return ImageFont.load_default()

TITLE_F = find_font(28); BODY_F = find_font(20)
outdir = "data/summary"; os.makedirs(outdir, exist_ok=True)

# QQ heatmaps were made up to these (live-checkpoint) times
QQ_TMAX = {"0p6": 68, "0p8": 56, "1p0": 60, "1p2": 72, "1p4": 72}
QQ_TH   = {"0p6": 0.6, "0p8": 0.8, "1p0": 1.0, "1p2": 1.2, "1p4": 1.4}
QQ_T2P  = {"0p6": 0.3, "0p8": 0.4, "1p0": 0.5, "1p2": -0.4, "1p4": -0.3}

plots = []
# 1) cumulative radiated energy vs T at L=20
plots.append(("data/erad/plots/cumulative_L20.png",
  "Cumulative radiated energy vs T  (L=20)",
  ["Study: energy radiation from +-1/2 static probe charges",
   "Observable: E_rad(T) = integral_0^T <J(x_D,t)> dt  (J = energy current T^01)",
   "N = 512... (erad) N=384 sites, ag = 0.2 (physical box 76.8), m/g = 0 (solid) and 1 (dashed)",
   "L = 20 sites (g*x = 4), theta=0 vacuum, probe charge q = +-1/2 (dtheta2pi = 0.5)",
   "T = 50, dt = 0.2, maxbond = 128 (E_rad converged: 96~256 within 0.007%)",
   "Detectors x_D = physical 15 and 20 (outside the string, boundary-clean)"]))
# 2) erad E-field heatmaps
for mg, mgv in (("0p0","0"),("1p0","1")):
    for L in (10,20,30,40):
        plots.append((f"data/erad/plots/heat_E_mg{mg}_L{L}.png",
          f"erad electric-field heatmap  (m/g={mgv}, L={L})",
          ["Study: energy radiation from +-1/2 static probe charges",
           "Observable: electric field E(x,t)  [1-2-1 spatial smoothing]",
           f"N = 384 sites, ag = 0.2, m/g = {mgv}",
           f"L = {L} sites (g*x = {L*0.2:g}), theta=0 vacuum, q = +-1/2 (dtheta2pi = 0.5)",
           "T = 50, dt = 0.2, maxbond = 128",
           ("Screening: imposed field disperses/radiates" if mgv=="0"
            else "Confinement: persistent flux tube, little radiated")]))
# 3) QQ E-field heatmaps
for tag in ("0p6","0p8","1p0","1p2","1p4"):
    plots.append((f"data/qq_sweep/plots/qqheat_E_theta{tag}.png",
      f"QQ electric-field heatmap  (theta/pi = {QQ_TH[tag]})",
      ["Study: charge-transfer QQ, Wilson-line quench |psi> = W|vac>",
       "Observable: electric field E(x,t)  [1-2-1 spatial smoothing]",
       "N = 512 sites, ag = 0.2 (physical box 102.4), m/g = 1",
       f"theta/pi = {QQ_TH[tag]}  (theta2pi folded = {QQ_T2P[tag]})",
       "Wilson line g*x = 2 (10 sites, centered)",
       f"T target = 90; shown up to t ~ {QQ_TMAX[tag]} (live checkpoint, run ongoing)",
       "dt = 0.2, maxbond = 256; detectors R = 15, 20 sites (for the QQ correlator)"]))

pages = []
PAGE_W = 1120
for path, title, params in plots:
    if not os.path.exists(path):
        print("MISSING:", path); continue
    img = Image.open(path).convert("RGB")
    s = 1000 / img.width
    img = img.resize((int(img.width*s), int(img.height*s)))
    text_h = 44 + 26*len(params)
    page = Image.new("RGB", (PAGE_W, img.height + text_h + 40), "white")
    page.paste(img, ((PAGE_W - img.width)//2, 20))
    d = ImageDraw.Draw(page); y = img.height + 30
    d.text((30, y), title, fill="black", font=TITLE_F); y += 40
    for ln in params:
        d.text((44, y), "- " + ln, fill=(20,20,20), font=BODY_F); y += 26
    pages.append(page)
    print("added:", path)

pdf = os.path.join(outdir, "summary.pdf")
pages[0].save(pdf, save_all=True, append_images=pages[1:], resolution=100.0)
zf = os.path.join(outdir, "summary.zip")
with zipfile.ZipFile(zf, "w", zipfile.ZIP_DEFLATED) as z:
    z.write(pdf, "summary.pdf")
    for path,_,_ in plots:
        if os.path.exists(path): z.write(path, "plots/"+os.path.basename(path))
print(f"wrote {pdf} ({len(pages)} pages) and {zf}")
