#' Plot MODGEE result
#'
#' Visualize cumulative extinction, species extinction probabilities, or
#' conditional expected extinction years from an existing MODGEE fit.
#'
#' @param result A result list returned by \code{fit_modgee()} using the current
#'   output format. Older objects must be refitted or updated to that format.
#' @param type A single string: \code{"curve"} (default),
#'   \code{"probability"}, or \code{"time"}. Full names are required.
#'
#' @return A \code{ggplot} object. The function neither saves a file nor
#'   explicitly prints the plot. Use \code{print()} inside loops and
#'   \code{ggplot2::ggsave()} to save it.
#'
#' @details
#' All plots use a white background, a closed black panel border, no grid,
#' larger axis text, and a two-column legend below the panel. Groups use
#' green, orange, blue and purple first; additional groups extend the palette.
#'
#' \strong{Cumulative extinction curves}
#'
#' \code{type = "curve"} plots \code{result$curve$time} against
#' \code{result$curve$cumulative_extinction_rate}. The axes are labelled
#' \emph{Year} and \emph{Cumulative extinction rate}. Rates are shown as
#' percentages. The overall curve is a thick black line. For grouped fits,
#' coloured group curves from \code{result$group_results$curve} are also
#' shown, and the overall curve is labelled \emph{Overall} in the legend.
#' The fitted overall curve is used directly, not an unweighted average
#' of group curves. No smoothing is added. When bootstrap intervals are
#' available, a grey band shows the overall pointwise confidence interval
#' from \code{result$bootstrap$curve}, at the confidence level used in fitting.
#' Group curves retain their original colours without additional bands.
#' Missing interval bounds leave gaps in the band.
#'
#' \strong{Extinction probability histograms}
#'
#' \code{type = "probability"} counts recorded species using
#' \code{result$species$p_extinct}. Bins have width 0.1, with boundaries
#' \code{[0, 0.1]}, \code{(0.1, 0.2]}, ..., \code{(0.9, 1]}.
#' The axes are \emph{Extinction probability} and \emph{Number of species}.
#' Grouped fits use stacked colours, so each bar's total is the total
#' number of species in that interval. Counts are not normalized by group size.
#'
#' \strong{Expected extinction years}
#'
#' \code{type = "time"} plots \code{estimated_extinction_year} against
#' \code{p_extinct}, with one point per eligible species. The axes are
#' \emph{Expected extinction year} and \emph{Extinction probability}.
#' Estimated years are conditional means given extinction, not confirmed
#' extinction dates. They are read directly from the species table and are
#' already mapped to the numeric time labels supplied during fitting.
#' If interval indices were supplied instead of years, the plot retains those
#' indices despite the default year labels; labels can be changed with
#' \code{ggplot2::labs()}. No jitter or regression line is added.
#'
#' \strong{Groups and missing estimates}
#'
#' Group membership comes from \code{result$species$belonged_group}.
#' Group order and the colour mapping use \code{result$group_results$summary}.
#' The same group has the same colour across the three plot types.
#' All-zero species are excluded from species-level plots. Groups containing
#' only all-zero species have no fitted curve and no recorded species to plot.
#' Dark extinction estimates are not converted into artificial species records.
#'
#' For probability histograms, missing or non-finite estimates are omitted
#' with a warning giving their count. Time plots silently omit species with
#' missing or non-finite estimates or zero extinction probability, because
#' their conditional extinction time is undefined.
#' Missing curve estimates, unknown groups, duplicate species identifiers,
#' invalid probabilities, or no eligible observations result in an error.
#' Only the fields needed by the requested plot are required.
#'
#' @importFrom rlang .data
#' @export
#'
#' @seealso
#' \strong{Functions:}
#' \code{\link{read_matrix}}, \code{\link{fit_modgee}}, \code{\link{save_result}}
#'
#' \strong{Examples:}
#' \code{\link{example_file}}, \code{\link{example_M}}, \code{\link{example_groups}}, \code{\link{example_result}}, \code{\link{example_grouped_result}}, \code{\link{example_bootstrap_result}}
#' @section Usage Examples:
#' \preformatted{
#' plot_modgee(example_result)
#' plot_modgee(example_result, type = "probability")
#' plot_modgee(example_result, type = "time")
#'
#' plot_modgee(example_grouped_result)
#' plot_modgee(example_grouped_result, type = "probability")
#' plot_modgee(example_grouped_result, type = "time")
#'
#' # Precomputed bootstrap fit with an overall confidence band.
#' plot_modgee(example_bootstrap_result)
#'
#' # Save a plot to a PNG file.
#' p <- plot_modgee(example_bootstrap_result)
#' ggplot2::ggsave("modgee_curve.png", plot = p,
#'                 width = 8, height = 6, units = "in", dpi = 300)
#' }
plot_modgee <- function(result, type = "curve") {
  choices <- c("curve", "probability", "time")
  if (!is.character(type) || length(type) != 1L ||
      is.na(type) || !type %in% choices)
    stop('type must be "curve", "probability", or "time".', call. = FALSE)
  if (!is.list(result) || is.data.frame(result))
    stop("result must be a fit_modgee result list.", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Install ggplot2 before using plot_modgee().", call. = FALSE)
  if (utils::packageVersion("ggplot2") < "3.4.0")
    stop("plot_modgee requires ggplot2 >= 3.4.0.", call. = FALSE)

  require_columns <- function(x, cols, label) {
    if (!is.data.frame(x) || !all(cols %in% names(x)))
      stop(label, " must be a data frame containing: ",
           paste(cols, collapse = ", "), ".", call. = FALSE)
  }
  probability <- function(x, label) {
    if (!is.numeric(x)) stop(label, " must be numeric.", call. = FALSE)
    if (any(is.finite(x) & (x < -1e-8 | x > 1 + 1e-8)))
      stop(label, " contains values outside [0, 1].", call. = FALSE)
    x[is.finite(x)] <- pmin(1, pmax(0, x[is.finite(x)]))
    x
  }
  valid_labels <- function(x) {
    is.character(x) && !anyNA(x) && all(nzchar(trimws(x)))
  }
  percent <- function(x) paste0(format(round(100 * x, 6), trim = TRUE), "%")
  grouped <- !is.null(result$group_results)
  group_names <- character()
  palette <- character()
  if (grouped) {
    if (!is.list(result$group_results))
      stop("result$group_results must be a list.", call. = FALSE)
    gs <- result$group_results$summary
    require_columns(gs, "group", "result$group_results$summary")
    group_names <- as.character(gs$group)
    if (!length(group_names) || !valid_labels(group_names) || anyDuplicated(group_names))
      stop("Group names must be nonempty and unique.", call. = FALSE)
    base_colours <- c("#238443", "#E08214", "#386CB0", "#984EA3",
                      "#D55E00", "#56B4E9", "#882255", "#999933")
    n <- length(group_names)
    palette <- stats::setNames(if (n <= length(base_colours)) base_colours[seq_len(n)]
                              else grDevices::hcl.colors(n, "Dark 3"), group_names)
  }
  theme <- ggplot2::theme_bw(base_size = 14, base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.7),
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black", size = 12),
      axis.title = ggplot2::element_text(size = 14),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 9)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 9)),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      legend.key = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 12, 10, 10)
    )

  if (type == "curve") {
    cols <- c("time", "cumulative_extinction_rate")
    overall <- result$curve
    require_columns(overall, cols, "result$curve")
    check_curve <- function(x, label) {
      if (!nrow(x) || !is.numeric(x$time) || any(!is.finite(x$time)))
        stop(label, " must contain finite numeric times.", call. = FALSE)
      x$cumulative_extinction_rate <- probability(x$cumulative_extinction_rate, label)
      if (any(!is.finite(x$cumulative_extinction_rate)) || anyDuplicated(x$time))
        stop(label, " contains missing rates or duplicate times.", call. = FALSE)
      x <- x[order(x$time), , drop = FALSE]
      if (any(diff(x$cumulative_extinction_rate) < -1e-8))
        stop(label, " is not cumulative (rates decrease).", call. = FALSE)
      x
    }
    overall <- check_curve(overall, "Overall curve")
    overall <- overall[, cols, drop = FALSE]
    # Internal IDs cannot collide with user-supplied group names.
    overall$series <- "overall"
    curves <- overall
    colours <- c(overall = "black")
    widths <- c(overall = 1.3)
    labels <- c(overall = "Overall")
    if (grouped) {
      gc <- result$group_results$curve
      require_columns(gc, c("group", cols), "result$group_results$curve")
      gc$group <- as.character(gc$group)
      if (!nrow(gc) || !valid_labels(gc$group) || any(!gc$group %in% group_names))
        stop("Group curves contain missing or unknown groups.", call. = FALSE)
      active <- group_names[group_names %in% gc$group]
      parts <- lapply(active, function(g) {
        x <- check_curve(gc[gc$group == g, , drop = FALSE], paste("Curve for", g))
        if (!identical(as.numeric(x$time), as.numeric(overall$time)))
          stop("Group and overall curve times must match.", call. = FALSE)
        x <- x[, cols, drop = FALSE]
        x$series <- paste0("group_", match(g, group_names))
        x
      })
      # Overall is drawn last so its thicker line stays visible.
      curves <- rbind(do.call(rbind, parts), overall)
      ids <- paste0("group_", match(active, group_names))
      colours <- c(stats::setNames(unname(palette[active]), ids), colours)
      widths <- c(stats::setNames(rep(0.8, length(ids)), ids), widths)
      labels <- c(stats::setNames(ifelse(active == "Overall", "Overall (group)", active), ids), labels)
    }
    curves$series <- factor(curves$series, levels = names(colours))
    band <- NULL
    if (!is.null(result$bootstrap)) {
      ci <- result$bootstrap$curve
      require_columns(ci, c("time", "lower", "upper"), "result$bootstrap$curve")
      if (!is.numeric(ci$time) || any(!is.finite(ci$time)) ||
          anyDuplicated(ci$time) || !setequal(ci$time, overall$time))
        stop("Bootstrap and overall curve times must match.", call. = FALSE)
      ci <- ci[match(overall$time, ci$time), , drop = FALSE]
      ci$lower <- probability(ci$lower, "Bootstrap lower bounds")
      ci$upper <- probability(ci$upper, "Bootstrap upper bounds")
      valid <- is.finite(ci$lower) & is.finite(ci$upper)
      if (any(ci$lower[valid] > ci$upper[valid]))
        stop("Bootstrap lower bounds exceed upper bounds.", call. = FALSE)
      if (any(!valid))
        warning(sum(!valid), " bootstrap intervals omitted because bounds are missing or non-finite.", call. = FALSE)
      # Separate runs prevent ribbons from spanning missing intervals.
      ci$segment <- cumsum(!valid)
      band <- ci[valid, , drop = FALSE]
    }
    upper <- min(1, max(0.05, max(c(curves$cumulative_extinction_rate,
                                  band$upper)) * 1.05))
    p <- ggplot2::ggplot(curves, ggplot2::aes(
      x = .data$time, y = .data$cumulative_extinction_rate,
      colour = .data$series, linewidth = .data$series, group = .data$series))
    if (!is.null(band) && nrow(band)) {
      p <- p + ggplot2::geom_ribbon(
        data = band,
        mapping = ggplot2::aes(x = .data$time, ymin = .data$lower,
                               ymax = .data$upper, group = .data$segment),
        inherit.aes = FALSE, fill = "#222222", alpha = 0.18,
        colour = NA, show.legend = FALSE)
    }
    p <- p +
      ggplot2::geom_line() +
      ggplot2::scale_colour_manual(values = colours, breaks = names(colours), labels = labels) +
      ggplot2::scale_linewidth_manual(values = widths, breaks = names(widths), labels = labels) +
      ggplot2::scale_y_continuous(labels = percent, limits = c(0, upper)) +
      ggplot2::labs(x = "Year", y = "Cumulative extinction rate",
                    colour = "Group", linewidth = "Group") + theme
    p <- p + ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = 2, byrow = TRUE),
      linewidth = ggplot2::guide_legend(ncol = 2, byrow = TRUE)
    )
    if (!grouped) p <- p + ggplot2::theme(legend.position = "none")
    return(p)
  }

  sp <- result$species
  required <- c("species", "detected_timesteps", "p_extinct")
  if (type == "time") required <- c(required, "estimated_extinction_year")
  if (grouped) required <- c(required, "belonged_group")
  require_columns(sp, required, "result$species")
  if (!valid_labels(as.character(sp$species)) || anyDuplicated(sp$species))
    stop("Species identifiers must be nonempty and unique.", call. = FALSE)
  counts <- sp$detected_timesteps
  if (!is.numeric(counts) || any(!is.finite(counts)) ||
      any(counts < 0 | counts != floor(counts)))
    stop("detected_timesteps must contain nonnegative integers.", call. = FALSE)
  sp <- sp[counts > 0, , drop = FALSE]
  sp$p_extinct <- probability(sp$p_extinct, "Species probabilities")
  if (grouped) {
    sp$belonged_group <- as.character(sp$belonged_group)
    if (!valid_labels(sp$belonged_group) || any(!sp$belonged_group %in% group_names))
      stop("Recorded species contain missing or unknown groups.", call. = FALSE)
    sp$group <- factor(sp$belonged_group, levels = group_names)
  } else sp$group <- factor(rep("All species", nrow(sp)))
  keep <- is.finite(sp$p_extinct)
  if (type == "time") {
    if (!is.numeric(sp$estimated_extinction_year))
      stop("estimated_extinction_year must be numeric.", call. = FALSE)
    keep <- keep & is.finite(sp$estimated_extinction_year) & sp$p_extinct > 0
  }
  if (type != "time" && any(!keep)) warning(sum(!keep), " recorded species omitted because required estimates are missing, non-finite, or undefined.", call. = FALSE)
  sp <- sp[keep, , drop = FALSE]
  if (!nrow(sp)) stop("No species have valid estimates for this plot.", call. = FALSE)

  if (type == "probability") {
    sp$bin <- cut(sp$p_extinct, breaks = seq(0, 1, by = 0.1),
                  include.lowest = TRUE, right = TRUE, labels = FALSE)
    bins <- as.data.frame(table(bin = factor(sp$bin, levels = 1:10),
                                group = droplevels(sp$group)))
    bins$midpoint <- (as.integer(bins$bin) - 0.5) / 10
    p <- ggplot2::ggplot(bins, ggplot2::aes(x = .data$midpoint, y = .data$Freq,
                                          fill = .data$group)) +
      ggplot2::geom_col(width = 0.1, position = "stack") +
      ggplot2::scale_x_continuous(labels = percent, breaks = seq(0, 1, by = 0.1),
                                 limits = c(0, 1), expand = ggplot2::expansion(mult = 0)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.06)), breaks = function(limits) {
        b <- pretty(limits); b[b >= 0 & abs(b - round(b)) < 1e-8]
      }) + ggplot2::labs(x = "Extinction probability", y = "Number of species", fill = "Group")
    if (grouped) p <- p + ggplot2::scale_fill_manual(values = palette)
    else p <- p + ggplot2::scale_fill_manual(values = "#238443", guide = "none")
  } else {
    p <- ggplot2::ggplot(sp, ggplot2::aes(x = .data$estimated_extinction_year,
                                        y = .data$p_extinct, colour = .data$group)) +
      ggplot2::geom_point(size = 1.8, alpha = 0.6) +
      ggplot2::scale_y_continuous(labels = percent, limits = c(0, 1)) +
      ggplot2::labs(x = "Expected extinction year", y = "Extinction probability", colour = "Group")
    if (grouped) p <- p + ggplot2::scale_colour_manual(values = palette)
    else p <- p + ggplot2::scale_colour_manual(values = "#238443", guide = "none")
  }
  p <- p + theme + ggplot2::guides(
    fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE),
    colour = ggplot2::guide_legend(ncol = 2, byrow = TRUE)
  )
  if (!grouped) p <- p + ggplot2::theme(legend.position = "none")
  p
}
