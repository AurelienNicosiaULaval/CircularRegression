# Circular diagnostic helpers shared by fitted-model plot methods.

utils::globalVariables(".data")

.require_circular_diagnostic_packages <- function(arrange = FALSE) {
  required <- c("ggplot2", "ggcircular")
  if (isTRUE(arrange)) required <- c(required, "gridExtra")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Circular diagnostic plots require ",
      paste(sprintf("'%s'", missing), collapse = ", "),
      ". Install the missing package(s) and try again.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.radian_breaks <- function() c(-pi, -pi / 2, 0, pi / 2, pi)

.radian_labels <- function() c("-pi", "-pi/2", "0", "pi/2", "pi")

.diagnostic_theme <- function() {
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 11, hjust = 0),
    plot.subtitle = ggplot2::element_text(size = 8.5, hjust = 0),
    axis.title = ggplot2::element_text(size = 9),
    axis.text = ggplot2::element_text(size = 8)
  )
}

.shared_circular_bandwidth <- function(residuals, bandwidth = NULL) {
  if (!is.null(bandwidth)) {
    if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
        !is.finite(bandwidth) || bandwidth <= 0) {
      stop("'bandwidth' must be a single finite positive number.", call. = FALSE)
    }
    return(as.numeric(bandwidth))
  }

  residuals <- as.numeric(residuals)
  residuals <- residuals[is.finite(residuals)]
  if (length(residuals) == 0L) {
    stop("At least one finite residual is required to select a bandwidth.", call. = FALSE)
  }

  residuals <- ggcircular::normalize_angle(.wrap_angle(residuals))
  rbar <- ggcircular::mean_resultant_length(residuals, na.rm = TRUE)
  kappa <- suppressWarnings(circular::A1inv(rbar))
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    kappa <- 1
  }
  kappa <- max(0.25, kappa)
  1 / sqrt(kappa)
}

.shared_rose_radial_limit <- function(plots) {
  values <- unlist(
    lapply(
      plots,
      function(plot) {
        built <- ggplot2::ggplot_build(plot)
        unlist(
          lapply(
            built$data,
            function(layer) {
              c(layer$y, layer$ymax)
            }
          ),
          use.names = FALSE
        )
      }
    ),
    use.names = FALSE
  )
  values <- values[is.finite(values) & values >= 0]
  if (length(values) == 0L) return(1)
  1.02 * max(1, values)
}

.apply_rose_radial_limit <- function(plot, radial_limit) {
  if (!is.numeric(radial_limit) || length(radial_limit) != 1L ||
      !is.finite(radial_limit) || radial_limit <= 0) {
    stop("'radial_limit' must be a single finite positive number.", call. = FALSE)
  }
  plot + ggplot2::scale_y_continuous(
    limits = c(0, radial_limit),
    expand = ggplot2::expansion(mult = c(0, 0))
  )
}

.vonmises_cdf <- function(angle, kappa) {
  angle <- as.numeric(angle)
  if (any(!is.finite(angle)) || any(angle < 0 | angle > 2 * pi)) {
    stop("CDF angles must lie in [0, 2*pi].", call. = FALSE)
  }
  kappa <- as.numeric(kappa)
  if (length(kappa) == 1L) kappa <- rep(kappa, length(angle))
  if (length(kappa) != length(angle) || any(!is.finite(kappa)) || any(kappa < 0)) {
    stop("A finite non-negative von Mises concentration is required.", call. = FALSE)
  }

  cdf <- numeric(length(angle))
  high_concentration <- kappa > 100
  if (any(high_concentration)) {
    # Scaled quadrature retains the von Mises density rather than replacing it
    # with a normal approximation. A 12-standard-deviation tail is negligible.
    cdf[high_concentration] <- vapply(which(high_concentration), function(i) {
      k <- kappa[i]
      d <- angle[i] - pi
      area <- stats::integrate(function(t) {
        exp(-2 * k * sin(t / (2 * sqrt(k)))^2) /
          (2 * pi * besselI(k, 0, expon.scaled = TRUE) * sqrt(k))
      }, 0, min(abs(d) * sqrt(k), 12), rel.tol = 1e-10)$value
      0.5 + sign(d) * area
    }, numeric(1))
  }
  if (any(!high_concentration)) {
    low_index <- which(!high_concentration)
    if (length(unique(kappa[low_index])) == 1L) {
      cdf[low_index] <- suppressWarnings(circular::pvonmises(
        angle[low_index],
        mu = pi,
        kappa = kappa[low_index[[1L]]],
        from = 0
      ))
    } else {
      cdf[low_index] <- vapply(
        low_index,
        function(i) suppressWarnings(circular::pvonmises(
          angle[[i]],
          mu = pi,
          kappa = kappa[[i]],
          from = 0
        )),
        numeric(1)
      )
    }
  }
  cdf[angle == 0] <- 0
  cdf[angle == 2 * pi] <- 1
  pmin(pmax(cdf, 0), 1)
}

.vonmises_pit <- function(residuals, kappa) {
  residuals <- .wrap_angle(as.numeric(residuals))
  kappa <- as.numeric(kappa)
  if (length(kappa) == 1L) kappa <- rep(kappa, length(residuals))
  if (length(kappa) != length(residuals) || any(!is.finite(kappa)) || any(kappa < 0)) {
    stop("A finite non-negative von Mises concentration is required for the PIT diagnostic.", call. = FALSE)
  }

  # circular::pvonmises() uses a CDF cut. Shifting the residual by pi fixes the
  # cut at the antipode of the expected residual direction (zero).
  shifted <- residuals + pi
  pit <- .vonmises_cdf(shifted, kappa)
  eps <- sqrt(.Machine$double.eps)
  pmin(pmax(pit, eps), 1 - eps)
}

.vonmises_bin_profile <- function(kappa, bins) {
  if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) || bins < 2 || bins != floor(bins)) {
    stop("'bins' must be a finite integer greater than one.", call. = FALSE)
  }
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) || kappa < 0) {
    stop("A finite non-negative von Mises concentration is required.", call. = FALSE)
  }

  edges <- seq(0, 2 * pi, length.out = as.integer(bins) + 1L)
  grid_density <- max(20 * sqrt(kappa), 64 / (2 * pi))
  mass <- vapply(
    seq_len(length(edges) - 1L),
    function(i) {
      n_grid <- max(33L, ceiling((edges[[i + 1L]] - edges[[i]]) * grid_density) + 1L)
      grid <- seq(edges[[i]], edges[[i + 1L]], length.out = n_grid)
      density <- suppressWarnings(circular::dvonmises(grid, mu = 0, kappa = kappa))
      sum(diff(grid) * (utils::head(density, -1L) + utils::tail(density, -1L)) / 2)
    },
    numeric(1)
  )
  mass <- mass / sum(mass)

  data.frame(
    residual = edges,
    density = c(mass / diff(edges), utils::tail(mass / diff(edges), 1L))
  )
}

.residual_rose_plot <- function(
  residuals,
  kappa = NULL,
  bins = 24,
  bandwidth = NULL,
  radial_limit = NULL,
  rug = TRUE,
  title = "Circular residuals",
  subtitle = NULL
) {
  .require_circular_diagnostic_packages()
  residual_data <- data.frame(residual = ggcircular::normalize_angle(.wrap_angle(residuals)))

  plot <- ggplot2::ggplot(residual_data, ggplot2::aes(x = .data$residual)) +
    ggcircular::geom_rose(
      bins = bins,
      normalize = "density",
      fill = "grey80",
      colour = "grey35"
    ) +
    ggcircular::geom_circular_density(
      bw = bandwidth,
      colour = "#2C7FB8",
      linewidth = 0.8
    ) +
    ggcircular::geom_mean_direction(colour = "#D95F0E", linewidth = 0.7) +
    ggcircular::scale_x_circular_radians(
      limits = c(0, 2 * pi),
      breaks = c(0, pi / 2, pi, 3 * pi / 2)
    ) +
    ggcircular::coord_circular() +
    ggcircular::theme_circular() +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = "Density") +
    .diagnostic_theme()

  if (isTRUE(rug)) {
    plot <- plot + ggcircular::geom_circular_rug(colour = "grey20", alpha = 0.35)
  }
  if (!is.null(kappa)) {
    kappa <- as.numeric(kappa)
    if (length(kappa) == 1L && is.finite(kappa) && kappa >= 0) {
      model_profile <- .vonmises_bin_profile(kappa, bins)
      plot <- plot + ggplot2::geom_step(
        data = model_profile,
        mapping = ggplot2::aes(x = .data$residual, y = .data$density),
        inherit.aes = FALSE,
        colour = "#D95F0E",
        linewidth = 0.9,
        direction = "hv"
      )
    }
  }
  if (!is.null(radial_limit)) {
    plot <- .apply_rose_radial_limit(plot, radial_limit)
  }
  plot
}

.pit_qq_plot <- function(residuals, kappa, title, subtitle = NULL) {
  keep <- is.finite(residuals)
  if (length(kappa) > 1L) {
    if (length(kappa) != length(residuals)) {
      stop("The PIT diagnostic requires one concentration per residual.", call. = FALSE)
    }
    kappa <- kappa[keep]
  }
  residuals <- residuals[keep]
  if (!length(residuals)) {
    stop("At least one finite residual is required for the PIT diagnostic.", call. = FALSE)
  }
  pit <- sort(.vonmises_pit(residuals, kappa))
  diagnostic_data <- data.frame(expected = stats::ppoints(length(pit)), observed = pit)

  ggplot2::ggplot(diagnostic_data, ggplot2::aes(x = .data$expected, y = .data$observed)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey45", linetype = "dashed") +
    ggplot2::geom_point(colour = "#2C7FB8", alpha = 0.75) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Expected uniform quantile",
      y = "von Mises PIT quantile"
    ) +
    ggplot2::theme_minimal() +
    .diagnostic_theme()
}

.residual_vs_fitted_plot <- function(residual_data, title) {
  ggplot2::ggplot(residual_data, ggplot2::aes(x = .data$.fitted, y = .data$.resid)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey45", linetype = "dashed") +
    ggplot2::geom_point(colour = "#2C7FB8", alpha = 0.7) +
    ggcircular::scale_x_circular_radians(
      limits = c(-pi, pi),
      breaks = .radian_breaks(),
      labels = .radian_labels()
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-pi, pi),
      breaks = .radian_breaks(),
      labels = .radian_labels()
    ) +
    ggplot2::labs(title = title, x = "Fitted direction", y = "Wrapped residual") +
    ggplot2::theme_minimal() +
    .diagnostic_theme()
}

.residual_by_index_plot <- function(residual_data, title) {
  line <- if (sum(is.finite(residual_data$.resid)) > 1L) {
    ggplot2::geom_line(colour = "grey70", linewidth = 0.35, na.rm = TRUE)
  } else NULL
  ggplot2::ggplot(residual_data, ggplot2::aes(x = .data$.index, y = .data$.resid)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey45", linetype = "dashed") +
    line +
    ggplot2::geom_point(colour = "#2C7FB8", alpha = 0.7, size = 1.3, na.rm = TRUE) +
    ggplot2::scale_y_continuous(
      limits = c(-pi, pi),
      breaks = .radian_breaks(),
      labels = .radian_labels()
    ) +
    ggplot2::labs(title = title, x = "Observation order", y = "Wrapped residual") +
    ggplot2::theme_minimal() +
    .diagnostic_theme()
}

.plot_fixed_circular_diagnostics <- function(x, residuals, fitted, kappa, type, bins, title) {
  type <- match.arg(type, c("all", "rose", "pit", "fitted", "index"))
  .require_circular_diagnostic_packages(arrange = identical(type, "all"))
  residual_data <- data.frame(
    .fitted = .wrap_angle(as.numeric(fitted)),
    .resid = .wrap_angle(as.numeric(residuals)),
    .index = seq_along(residuals)
  )
  finite_residual <- is.finite(residual_data$.resid)
  if (!any(finite_residual)) {
    stop("At least one finite residual is required for circular diagnostics.", call. = FALSE)
  }
  note <- if (all(finite_residual)) NULL else sprintf(
    "Residuals are undefined for %d of %d observations.",
    sum(!finite_residual), nrow(residual_data)
  )
  pit_subtitle <- if (length(kappa) == 1L) {
    "Fitted common concentration"
  } else {
    "Fitted observation-specific concentration"
  }
  panel_builders <- list(
    rose = function() .residual_rose_plot(
      residuals = residual_data$.resid[finite_residual],
      kappa = if (length(kappa) == 1L) kappa else NULL,
      bins = bins,
      title = "Circular residuals",
      subtitle = if (length(kappa) == 1L) {
        "Rose, density, and fitted von Mises bin expectation"
      } else {
        "Rose and kernel density"
      }
    ),
    pit = function() .pit_qq_plot(residual_data$.resid, kappa, "PIT Q-Q plot", pit_subtitle),
    fitted = function() .residual_vs_fitted_plot(
      residual_data[finite_residual, , drop = FALSE], "Residuals vs fitted direction"
    ),
    index = function() .residual_by_index_plot(residual_data, "Residuals by observation order")
  )
  selected <- if (type == "all") names(panel_builders) else type
  panels <- lapply(panel_builders[selected], function(build) build())

  if (type == "all") {
    gridExtra::grid.arrange(
      grobs = panels,
      ncol = 2,
      top = grid::textGrob(title, gp = grid::gpar(fontsize = 12)),
      bottom = if (is.null(note)) NULL else grid::textGrob(note, gp = grid::gpar(fontsize = 9))
    )
  } else {
    print(panels[[type]] + ggplot2::labs(caption = note))
  }
  invisible(x)
}
