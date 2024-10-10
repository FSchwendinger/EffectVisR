
# Authors: Fabian Schwendinger & Eric Lichtenstein

#' Effect Size Plot - Version 1
#'
#' Creates a plot showing an arrow representing an effect size, with a colored range
#' representing the 95% confidence interval (CI). The color transition is smooth across the CI,
#' and a dashed white line separates positive and negative effects. The plot provides a clear
#' visual summary of the effect size and its confidence range.
#'
#' @param effect_size Numeric, the estimated effect size.
#' @param ci_lower Numeric, the lower bound of the 95% confidence interval.
#' @param ci_upper Numeric, the upper bound of the 95% confidence interval.
#' @param palette Character, color palette to use for the confidence range. Options are
#'   "viridis", "magma", "plasma", "cividis", or "grey".
#' @param arrow_color Character, color for the arrow indicating the effect size.
#' @param save_path Character, path to save the plot. If NULL, the plot will not be saved.
#' @param axis_title_size Numeric, font size for the axis title.
#' @param axis_label_size Numeric, font size for the axis labels.
#' @param reverse Logical, whether to reverse the color scale (default is FALSE).
#'
#' @import ggplot2
#' @importFrom viridisLite viridis
#' @importFrom grDevices gray.colors
#' @return A ggplot2 object.
#' @name es_plot_v1
#' @export
#' @examples
#' # Example usage of es_plot_v1
#' library(EffectVisR)
#'
#' # Define the effect size and confidence intervals
#' effect_size <- 0.4
#' ci_lower <- 0.1
#' ci_upper <- 0.7
#'
#' # Create the plot using the "viridis" color palette
#' EffectVisR::es_plot_v1(
#'   effect_size = effect_size,
#'   ci_lower = ci_lower,
#'   ci_upper = ci_upper,
#'   palette = "viridis",
#'   arrow_color = "black"
#' )
#'
#' # Create the plot using the "magma" color palette and save the plot
#' EffectVisR::es_plot_v1(
#'   effect_size = 0.4,
#'   ci_lower = 0.1,
#'   ci_upper = 0.7,
#'   palette = "magma",
#'   arrow_color = "black",
#'   save_path = "effect_size_plot.png"
#' )

utils::globalVariables(c("x", "y", "z"))

es_plot_v1 <- function(effect_size, ci_lower, ci_upper, palette = "viridis",
                             arrow_color = "black", save_path = NULL,
                             axis_title_size = 12, axis_label_size = 10) {

  # Helper function to cap effect sizes and CIs within -1 to 1 range
  cap_values <- function(x) {
    pmin(pmax(x, -1), 1)
  }

  # Cap the effect size and CIs between -1 and 1
  capped_effect_size <- cap_values(effect_size)
  capped_ci_lower <- cap_values(ci_lower)
  capped_ci_upper <- cap_values(ci_upper)

  # Map the capped values to angles between -90° and 90° (in radians)
  angle <- capped_effect_size * (pi / 2)
  ci_angle_lower <- capped_ci_lower * (pi / 2)
  ci_angle_upper <- capped_ci_upper * (pi / 2)

  # Create a grid of points in the region defined by the CI
  r <- seq(0, 1, length.out = 250)  # Radial values
  th <- seq(-pi / 2, pi / 2, length.out = 250)  # Full angle values for color gradient across entire range
  grid_data <- expand.grid(r = r, th = th)

  # Calculate x and y positions for points in polar coordinates
  grid_data$x <- grid_data$r * cos(grid_data$th)
  grid_data$y <- grid_data$r * sin(grid_data$th)

  # Function to get the color palette based on user input
  get_palette <- function(palette, n) {
    switch(palette,
           "viridis" = viridisLite::viridis(n, option = "viridis", direction = ifelse(reverse, -1, 1)),
           "magma" = viridisLite::viridis(n, option = "magma", direction = ifelse(reverse, -1, 1)),
           "plasma" = viridisLite::viridis(n, option = "plasma", direction = ifelse(reverse, -1, 1)),
           "cividis" = viridisLite::viridis(n, option = "cividis", direction = ifelse(reverse, -1, 1)),
           "grey" = grDevices::gray.colors(n, start = ifelse(reverse, 0.9, 0.1), end = ifelse(reverse, 0.1, 0.9)),
           viridisLite::viridis(n)  # Default palette if none is specified
    )
  }

  # Apply a gradient across the entire range but only show within the CI
  gradient_colors <- get_palette(palette, length(unique(grid_data$th)))

  # Color only within the CI
  grid_data$z <- ifelse(grid_data$th >= ci_angle_lower & grid_data$th <= ci_angle_upper,
                        gradient_colors[as.numeric(cut(grid_data$th, breaks = 250))],
                        NA)

  # Create a half-circle for the scale
  scale_circle_data <- data.frame(
    x = cos(seq(-pi / 2, pi / 2, length.out = 100)),
    y = sin(seq(-pi / 2, pi / 2, length.out = 100))
  )

  # Create the plot with points shaded using gradient colors within the 95% CI
  p <- ggplot2::ggplot(grid_data) +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = z), width = 0.01, height = 0.01, na.rm = TRUE) +  # Show colors only in CI
    ggplot2::scale_fill_identity() +

    # Add a half-circle to highlight the scale
    ggplot2::geom_path(data = scale_circle_data, ggplot2::aes(x = x, y = y), color = "black", linewidth = 1) +

    # Add an arrow to represent the effect size with user-defined color
    ggplot2::annotate("segment", x = 0, y = 0, xend = cos(angle), yend = sin(angle),
                      size = 2.5, color = arrow_color) +  # Thicker segment for needle
    ggplot2::annotate("point", x = cos(angle), y = sin(angle),
                      size = 5, color = arrow_color, shape = 21, fill = "white") +  # Circle at the tip

    # Axis lines along the radius using annotate
    ggplot2::annotate("segment", x = 0, y = 0, xend = cos(pi/2), yend = sin(pi/2), linetype = "dashed", color = "gray") +  # 90°
    ggplot2::annotate("segment", x = 0, y = 0, xend = cos(pi/4), yend = sin(pi/4), linetype = "dashed", color = "gray") +  # 45°
    ggplot2::annotate("segment", x = 0, y = 0, xend = cos(0), yend = sin(0), linetype = "dashed", color = "gray") +  # 0°
    ggplot2::annotate("segment", x = 0, y = 0, xend = cos(-pi/4), yend = sin(-pi/4), linetype = "dashed", color = "gray") +  # -45°
    ggplot2::annotate("segment", x = 0, y = 0, xend = cos(-pi/2), yend = sin(-pi/2), linetype = "dashed", color = "gray") +  # -90°

    # Labels at the end of each line using annotate
    ggplot2::annotate("text", x = cos(pi/2), y = sin(pi/2), label = "\u2265 1", vjust = -1, size = axis_label_size) +  # 90° -> 1
    ggplot2::annotate("text", x = cos(pi/4), y = sin(pi/4), label = "0.5", vjust = -1, size = axis_label_size) +  # 45° -> 0.5
    ggplot2::annotate("text", x = cos(0), y = sin(0), label = "0", hjust = -1, size = axis_label_size) +  # 0° -> 0
    ggplot2::annotate("text", x = cos(-pi/4), y = sin(-pi/4), label = "-0.5", vjust = 2, size = axis_label_size) +  # -45° -> -0.5
    ggplot2::annotate("text", x = cos(-pi/2), y = sin(-pi/2), label = "\u2264 -1", vjust = 2, size = axis_label_size) +  # -90° -> -1

    # Line at x = 0
    ggplot2::annotate("segment", x = 0, xend = 0, y = -1, yend = 1, color = "black", size = 1) +

    # Add a color scale to the left of the plot
    ggplot2::geom_tile(data = data.frame(y = seq(-1, 1, length.out = 100), z = get_palette(palette, 100)),
                       ggplot2::aes(x = -0.05, y = y, fill = z), width = 0.05, height = 0.02) +  # Color scale as vertical line

    ggplot2::expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) +  # Adjust plot limits to fit the color scale

    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::annotate("text", x = -.2, y = 0, label = "Effect size (95% CI)", angle = 90, hjust = 0.5, size = axis_title_size) +
    ggplot2::theme(legend.position = "none")

  # Save the plot to the specified path with transparent background, if provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(filename = save_path, plot = p, dpi = 400, bg = "transparent", width = 8, height = 8)
  }

  return(p)
}



