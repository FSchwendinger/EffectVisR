
# Authors: Fabian Schwendinger & Eric Lichtenstein
# Contributors: Denis Infanger

#' Visualize Effect Sizes with Confidence Intervals in Polar Plot Halves
#' Effect Size Plot - Version 2
#'
#' This function creates a polar plot half showing an arrow representing an effect size
#' and a colored fan representing the 95% confidence interval (CI). The plot supports both linear and logarithmic effect
#' sizes, allowing visualization of metrics such as odds ratios, hazard ratios, or standardized effect sizes.
#' A smooth color transition across the CI enhances interpretability, and the plot provides a clear visual
#' summary of the effect size and its uncertainty.
#'
#'Shading in Version 2
#'How it works: In v2, the shading varies based on the distance from the effect size (estimate):
#'The shading is darker near the effect size, representing higher certainty around the estimate.
#'The shading becomes lighter as you move further away from the effect size within the CI range.
#'Why this is useful: This method emphasizes the most likely region of the CI (near the estimate), while visually de-emphasizing areas near the CI boundaries.
#'When to choose v2:
#'If you want to focus attention on the precision of the estimate.
#'For a more informative and intuitive visualization of the CI's relationship to the effect size.
#'
#' @param effect_size Numeric, the estimated effect size.
#' @param ci_lower Numeric. The lower bound of the 95% confidence interval. Must be less than `ci_upper`.
#' @param ci_upper Numeric. The upper bound of the 95% confidence interval. Must be greater than `ci_lower`.
#' @param palette Character. The color palette to use for the confidence range. Options are:
#'   - `"viridis"`: A perceptually uniform colormap.
#'   - `"magma"`: A colormap suitable for perceptual uniformity.
#'   - `"plasma"`: A high-contrast colormap.
#'   - `"cividis"`: A colormap optimized for colorblind users.
#'   - `"grey"`: A greyscale palette.
#'   - `"RdYlGn"`: A diverging palette from red to green.
#'   Default is `"grey"`.
#' @param arrow_color Character. The color of the arrow representing the effect size. Default is "black".
#' @param save_path Character. Path to save the plot as an image. If `NULL` (default), the plot is not saved.
#' @param axis_title_size Numeric. Font size for the axis title. Default is 12.
#' @param axis_label_size Numeric. Font size for the axis labels. Default is 10.
#' @param reverse Logical. If `TRUE`, reverses the color scale. Default is `FALSE`.
#' @param limit Numeric. Specifies the range for capping the effect size plot on each end.
#'   For logarithmic scales, this represents the maximum value after log-transformation. Default is 1.
#' @param eff_type Character. Specifies the scale of the effect size:
#'   - `"linear"`: Effect sizes are interpreted directly.
#'   - `"log"`: Effect sizes are log-transformed (e.g., odds ratios, hazard ratios).
#'   Default is `"linear"`.
#' @param num_labels Numeric (1, 2, or 3). Specifies the number of axis labels to display. Default is 1.
#' @param area_arc Numeric. Specifies the proportion of the CI fan to display.
#'   - `0`: Displays a full circular fan.
#'   - Values closer to 1 produce narrower arcs.
#'   Must be within the range `[0, 1)`. Default is `0`.
#'
#' @seealso
#'   - [ggplot2](https://ggplot2.tidyverse.org/) for additional plotting options.
#'   - [viridisLite](https://cran.r-project.org/package=viridisLite) for color palette details.
#'   - [EffectVisR](https://github.com/FSchwendinger/EffectVisR) for related visualization tools.
#'
#' @import ggplot2
#' @importFrom viridisLite viridis
#' @importFrom grDevices gray.colors
#' @return A `ggplot2` object representing the effect size plot with confidence intervals. The plot
#'   can be further customized or saved using `ggplot2` functions.
#' @name es_plot_v2
#' @export
#' @examples
#' # Minimalistic Example: Create a basic plot with default settings
#' library(EffectVisR)
#' EffectVisR::es_plot_v2(
#'   effect_size = 0.5,
#'   ci_lower = 0.2,
#'   ci_upper = 0.8
#' )

#' # Example 1: Create a plot with a linear effect size and default settings
#' effect_size <- 0.5
#' ci_lower <- 0.2
#' ci_upper <- 0.8
#' EffectVisR::es_plot_v2(
#'   effect_size = effect_size,
#'   ci_lower = ci_lower,
#'   ci_upper = ci_upper,
#'   palette = "viridis",
#'   arrow_color = "blue",
#'   save_path = NULL,          # Do not save the plot
#'   axis_title_size = 12,      # Default axis title size
#'   axis_label_size = 10,      # Default axis label size
#'   reverse = FALSE,           # Do not reverse the palette
#'   limit = 1,                 # Default effect size limit
#'   eff_type = "linear",       # Use a linear scale
#'   num_labels = 1,            # Single label on the axis
#'   area_arc = 0               # Full fan
#' )

#' # Example 2: Create a plot with a logarithmic effect size and custom settings
#' EffectVisR::es_plot_v2(
#'   effect_size = 1.5,         # Effect size on a log scale
#'   ci_lower = 1.2,
#'   ci_upper = 2.0,
#'   palette = "RdYlGn",        # Use the "RdYlGn" palette
#'   arrow_color = "darkred",   # Use a dark red arrow
#'   save_path = "log_effect_size_plot.png",  # Save the plot
#'   axis_title_size = 14,      # Larger axis title font size
#'   axis_label_size = 12,      # Larger axis label font size
#'   reverse = TRUE,            # Reverse the color palette
#'   limit = 10,                # Custom effect size limit
#'   eff_type = "log",          # Use a logarithmic scale
#'   num_labels = 3,            # Display three labels
#'   area_arc = 0.95            # Narrow arc for CI area
#' )

#' # Example 3: Customizing the color scale with a grey palette
#' EffectVisR::es_plot_v2(
#'   effect_size = -0.3,        # Negative effect size
#'   ci_lower = -0.6,
#'   ci_upper = -0.1,
#'   palette = "grey",          # Use a greyscale palette
#'   arrow_color = "black",     # Black arrow
#'   save_path = NULL,          # Do not save the plot
#'   axis_title_size = 10,      # Smaller axis title font
#'   axis_label_size = 8,       # Smaller axis label font
#'   reverse = FALSE,           # Do not reverse the palette
#'   limit = 1,                 # Default limit
#'   eff_type = "linear",       # Use a linear scale
#'   num_labels = 2,            # Two labels on the axis
#'   area_arc = 0.5             # Semi-circular arc for CI area
#' )



utils::globalVariables(c("x", "y", "z"))

# --- Function definition ---


es_plot_v2 <- function(effect_size,
                       ci_lower,
                       ci_upper,
                       palette = c("grey", "viridis", "magma", "plasma", "cividis", "RdYlGn"),
                       arrow_color = "black",
                       reverse = FALSE,
                       save_path = NULL,
                       axis_title_size = 12,
                       axis_label_size = 10,
                       limit = 1,
                       eff_type = c("linear", "log"),
                       num_labels = 1,
                       area_arc = 0) {


  # Function that linearly scales values of x from range [min, max] to range [a, b]
  scale_lin <- function(x, a, b, min, max) {
    (((b - a) * (x - min)) / (max - min)) + a
  }

  # Helper function to cap effect sizes and CIs within provided limits
  cap_values <- function(x, limit) {
    pmin(pmax(x, -limit), limit)
  }

  # Get effect size type
  eff_type <- match.arg(eff_type)


  # --- Input validation ---


  # CI bounds check
  if (ci_lower >= ci_upper) {
    stop("Error: `ci_lower` must be less than `ci_upper`.")
  }

  # Effect size bounds check
  if (effect_size < ci_lower || effect_size > ci_upper) {
    stop("Error: `effect_size` must lie within the range [`ci_lower`, `ci_upper`].")
  }

  if (eff_type %in% "log" && limit <= 1) {
    stop("Plotting limit must exceed 1 for effect sizes on logarithmic scales.")
  }

  if (area_arc < 0 || area_arc >= 1) {
    stop("`area_arc` must be a numeric value between 0 and less than 1.")
  }

  if (num_labels < 1 || num_labels > 3) {
    stop("`num_labels` must be an integer between 1 and 3.")
  }

  # Limit must be positive
  limit <- abs(limit)

  # Because axis is on the log scale, transform limit to log scale too
  if (eff_type %in% "log") {
    limit <- log(limit)
  }

  # Cap the effect size and CIs between provided limits
  capped_effect_size <- cap_values(effect_size, limit)
  capped_ci_lower <- cap_values(ci_lower, limit)
  capped_ci_upper <- cap_values(ci_upper, limit)

  # Map the capped values to angles between -90° and 90° (in radians)
  angle <- capped_effect_size * (pi / (2 * limit))
  ci_angle_lower <- capped_ci_lower * (pi / (2 * limit))
  ci_angle_upper <- capped_ci_upper * (pi / (2 * limit))

  # Create a grid of points in the region defined by the CI
  r <- seq(area_arc, 1, length.out = 250)  # Radial values
  th <- seq(ci_angle_lower, ci_angle_upper, length.out = 250)  # Angle values between CI
  grid_data <- expand.grid(r = r, th = th)

  # Calculate x and y positions for points in polar coordinates
  grid_data$x <- grid_data$r * cos(grid_data$th)
  grid_data$y <- grid_data$r * sin(grid_data$th)

  # Get palette, allowing for partially matched inputs
  palette <- match.arg(palette)

  # Function to get the color palette based on user input
  get_palette <- function(palette, n) {
    switch(
      palette,
      "viridis" = viridisLite::viridis(n, option = "viridis", direction = ifelse(reverse, -1, 1)),
      "magma" = viridisLite::viridis(n, option = "magma", direction = ifelse(reverse, -1, 1)),
      "plasma" = viridisLite::viridis(n, option = "plasma", direction = ifelse(reverse, -1, 1)),
      "cividis" = viridisLite::viridis(n, option = "cividis", direction = ifelse(reverse, -1, 1)),
      "grey" = grDevices::gray.colors(
        n,
        start = ifelse(reverse, 0.9, 0.1),
        end = ifelse(reverse, 0.1, 0.9)
      ),
      "RdYlGn" = grDevices::colorRampPalette(c("#d73027", "#ffffbf", "#1a9850"))(n),
      viridisLite::viridis(n)  # Default palette if none is specified
    )
  }

  # Create a gradient where the color is darkest at the effect size and lightest at CI limits
  distance_from_estimate <- abs(grid_data$th - angle)
  grid_data$z <- get_palette(palette, length(unique(grid_data$th)))[as.numeric(cut(distance_from_estimate, breaks = 250))]

  # Create a half-circle for the scale
  scale_circle_data <- data.frame(x = cos(seq(-pi / 2, pi / 2, length.out = 100)), y = sin(seq(-pi / 2, pi / 2, length.out = 100)))

  # Add axis labels and ticks based on the effect type
  if (eff_type == "log") {
    axis_vec <- seq(-limit, limit, length.out = num_labels + 1)
    axis_labels <- exp(axis_vec)
  } else {
    axis_vec <- seq(-limit, limit, length.out = num_labels + 1)
    axis_labels <- axis_vec
  }

  axis_angles <- axis_vec * (pi / (2 * limit))

  ticks_data <- data.frame(
    x_start = cos(axis_angles),
    y_start = sin(axis_angles),
    x_end = cos(axis_angles) * 1.05,
    y_end = sin(axis_angles) * 1.05
  )


  # --- Plot creation ---


  # Create the plot with shading only within the 95% CI
  p <- ggplot2::ggplot(grid_data) +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = z),
                       width = 0.01,
                       height = 0.01) +  # Faster tile-based shading
    ggplot2::scale_fill_identity() +  # Use the exact fill specified in 'z'

    # Add a half-circle to highlight the scale
    ggplot2::geom_path(
      data = scale_circle_data,
      ggplot2::aes(x = x, y = y),
      color = "black",
      linewidth = 1
    ) +

    # Add an arrow to represent the effect size with user-defined color
    ggplot2::annotate(
      "segment",
      x = 0,
      y = 0,
      xend = cos(angle),
      yend = sin(angle),
      linewidth = 2.5,
      color = arrow_color
    ) +  # Thicker segment for needle
    ggplot2::annotate(
      "point",
      x = cos(angle),
      y = sin(angle),
      size = 5,
      color = arrow_color,
      shape = 21,
      fill = "white"
    )  # Circle at the tip

  # Axis lines along the radius using annotate
  if (eff_type %in% "log") {
    # Define a small positive lower limit to avoid log(0)
    min_log <- -limit
    max_log <- limit

    # Determine axis_vec_log based on num_labels
    if (num_labels == 3) {
      axis_vec_log <- c(
        min_log * 0.75,
        min_log * 0.5,
        min_log * 0.25,
        0,
        max_log * 0.25,
        max_log * 0.5,
        max_log * 0.75
      )
    } else if (num_labels == 2) {
      axis_vec_log <- c(min_log * 0.666,
                        min_log * 0.333,
                        0,
                        max_log * 0.333,
                        max_log * 0.666)
    } else if (num_labels == 1) {
      axis_vec_log <- c(0, min_log)
    }

    # Scale these log values to angles
    axis_angle_vec <- scale_lin(
      axis_vec_log,
      a = -pi / 2,
      b = pi / 2,
      min = min_log,
      max = max_log
    )


    # Convert log values back to original scale for labels
    axis_labels <- round(exp(axis_vec_log), 1)

  } else if (eff_type %in% "linear") {
    axis_vec <- seq(-limit + 0.5, limit - 0.5, by = 0.5)
    axis_angle_vec <- scale_lin(
      axis_vec,
      a = -pi / 2,
      b = pi / 2,
      min = -limit,
      max = limit
    )
    axis_labels <- axis_vec
  }

  # Calculate hjust and vjust for label alignment
  vjust_vec <- 1 / 2 - 1 / 1.15 * sin(axis_angle_vec)
  hjust_vec <- 1 / 2 - 1 / 1.15 * cos(axis_angle_vec)

  # Prepare data for ticks at the label positions based on axis_angle_vec
  tick_length <- 0.025 # Length of ticks (you can adjust this)
  ticks_data <- data.frame(
    x_start = cos(axis_angle_vec),
    y_start = sin(axis_angle_vec),
    x_end = cos(axis_angle_vec) * (1 + tick_length),
    # Extend outward by tick_length
    y_end = sin(axis_angle_vec) * (1 + tick_length)   # Extend outward by tick_length
  )

  # Add dashed axis lines and labels
  for (i in seq_along(axis_angle_vec)) {
    p <- p +
      ggplot2::annotate(
        "segment",
        x = 0,
        y = 0,
        xend = cos(axis_angle_vec[i]),
        yend = sin(axis_angle_vec[i]),
        linetype = "dashed",
        color = "gray"
      ) +
      ggplot2::annotate(
        "text",
        x = cos(axis_angle_vec[i]),
        y = sin(axis_angle_vec[i]),
        label = axis_labels[i],
        hjust = hjust_vec[i],
        vjust = vjust_vec[i],
        size = axis_label_size
      )
  }

  # Labels for min and max
  p <- p + ggplot2::annotate(
    "text",
    x = cos(pi / 2),
    y = sin(pi / 2),
    label = paste0("\u2265 ", ifelse(eff_type %in% "linear", limit, exp(limit))),
    vjust = -1,
    size = axis_label_size
  ) +  # 90° -> limit
    ggplot2::annotate(
      "text",
      x = cos(-pi / 2),
      y = sin(-pi / 2),
      label = ifelse(eff_type %in% "linear", paste0("\u2264 ", -limit), ""),
      vjust = 2,
      size = axis_label_size
    ) +  # -90° -> -limit

    # Add ticks at the labelled positions on the circle
    ggplot2::geom_segment(
      data = ticks_data,
      ggplot2::aes(
        x = x_start,
        y = y_start,
        xend = x_end,
        yend = y_end
      ),
      color = "black"
    ) +

    # Line at x = 0
    ggplot2::annotate(
      "segment",
      x = 0,
      xend = 0,
      y = -1,
      yend = 1,
      color = "black",
      linewidth = 1
    ) +

    ggplot2::expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) +  # Adjust plot limits to fit the top label

    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::annotate(
      "text",
      x = -.2,
      y = 0,
      label = "Effect size (95% CI)",
      angle = 90,
      hjust = 0.5,
      size = axis_title_size
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(
        t = 5.5,
        r = 5.5,
        b = 5.5,
        l = -100
      )  # Adjust the left margin
    )


  # ---Save the plot (optional) ---


    if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      dpi = 400,
      bg = "transparent",
      width = 8,
      height = 8
    )
  }

  return(p)
}

# es_plot_v2(effect_size = 0.9, 0.2, 2, eff_type = "log", area_arc = 0.95, limit = 20, num_labels = 3, palette = "RdYlGn")

