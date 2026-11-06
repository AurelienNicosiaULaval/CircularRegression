# Run from the package root to render the README and pkgdown logo.
rsvg::rsvg_png(
  "figure/hex-logo.svg",
  file = "man/figures/logo.png",
  width = 1600,
  height = 1600
)
