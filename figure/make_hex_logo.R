draw_hex_logo <- function() {
  op <- par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", bg = "transparent")
  on.exit(par(op), add = TRUE)

  plot.new()
  plot.window(xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.15), asp = 1)

  th <- seq(0, 2 * pi, length.out = 7) + pi / 6
  hx <- cos(th)
  hy <- sin(th)

  polygon(1.06 * hx, 1.06 * hy, col = "#0B2239", border = NA)
  polygon(0.95 * hx, 0.95 * hy, col = "#F7FAFC", border = NA)

  r_circle <- 0.42

  symbols(
    0,
    0,
    circles = r_circle,
    inches = FALSE,
    add = TRUE,
    fg = "#0B2239",
    lwd = 5,
    bg = NA
  )

  axis_angles <- c(0, pi / 2, pi, 3 * pi / 2)
  segments(
    x0 = r_circle * cos(axis_angles),
    y0 = r_circle * sin(axis_angles),
    x1 = -r_circle * cos(axis_angles),
    y1 = -r_circle * sin(axis_angles),
    col = "#C3D3E0",
    lwd = 1.5
  )

  r_arc <- 0.57
  draw_arc <- function(start_deg, end_deg, col, lwd = 9) {
    tt <- seq(start_deg, end_deg, length.out = 120) * pi / 180
    lines(r_arc * cos(tt), r_arc * sin(tt), col = col, lwd = lwd, lend = "round")
  }

  # Keep the lower center free so the package title does not collide
  draw_arc(210, 245, "#2A9D8F", lwd = 9)
  draw_arc(295, 330, "#2A9D8F", lwd = 9)
  draw_arc(20, 145, "#E76F51", lwd = 9)

  target <- 55 * pi / 180
  arrows(
    0,
    0,
    0.42 * cos(target),
    0.42 * sin(target),
    length = 0.12,
    lwd = 4,
    col = "#0B2239"
  )

  pts <- c(15, 80, 145, 230, 300) * pi / 180
  points(r_circle * cos(pts), r_circle * sin(pts), pch = 21, cex = 1.25, lwd = 1.2,
         bg = "#F7FAFC", col = "#0B2239")

  pkg_name <- "CircularRegression"
  avail_width <- 1.16
  pkg_cex <- min(
    1.12,
    0.92 * avail_width / strwidth(pkg_name, cex = 1, font = 2, units = "user")
  )
  text(0, -0.62, pkg_name, cex = pkg_cex, col = "#0B2239", font = 2)
}

png("figure/hex-logo.png", width = 1600, height = 1600, res = 320, bg = "transparent")
draw_hex_logo()
dev.off()

ok_svg <- FALSE
if (file.exists("figure/hex-logo.svg")) {
  file.remove("figure/hex-logo.svg")
}

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite("figure/hex-logo.svg", width = 6, height = 6, bg = "transparent")
  draw_hex_logo()
  dev.off()
  ok_svg <- file.exists("figure/hex-logo.svg")
} else {
  suppressWarnings(
    tryCatch(
      {
        svg("figure/hex-logo.svg", width = 6, height = 6, bg = "transparent")
        draw_hex_logo()
        dev.off()
        ok_svg <- file.exists("figure/hex-logo.svg")
      },
      error = function(e) {
        ok_svg <<- FALSE
      }
    )
  )
}

if (!ok_svg) {
  message("SVG export skipped on this machine (no svglite/cairo backend).")
}
