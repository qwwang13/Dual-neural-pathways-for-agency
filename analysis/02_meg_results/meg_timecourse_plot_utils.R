suppressPackageStartupMessages({
  library(tidyverse)
  library(grid)
})

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

roi12 <- c(
  "DLPFCR", "PreSMAL", "PreSMAR", "SMAR",
  "InsulaR", "PrecentralR", "PostcentralL", "PostcentralR",
  "SPGL", "PrecuneusR", "AGR", "MTGR"
)

pal_action_line <- c(
  "Voluntary"   = "#1F86D1",
  "Involuntary" = "#FF9A2F"
)

pal_action_fill <- c(
  "Voluntary"   = "#9FD0F2",
  "Involuntary" = "#FFC980"
)

pal_outcome_line <- c(
  "Voluntary"   = "#74C95C",
  "Involuntary" = "#F783A0"
)

pal_outcome_fill <- c(
  "Voluntary"   = "#B4E39A",
  "Involuntary" = "#F9D1DC"
)

y_breaks_1dp <- function(x) {
  rng <- range(x, finite = TRUE)
  if (!all(is.finite(rng))) return(NULL)
  pretty(rng, n = 4)
}

read_tc_summary <- function(path, condition_label, roi_order = roi12) {
  dat <- read_csv(path, show_col_types = FALSE)

  if (!("Time" %in% names(dat)))   stop("Missing column: Time in ", path)
  if (!("region" %in% names(dat))) stop("Missing column: region in ", path)
  if (!("Value" %in% names(dat)))  stop("Missing column: Value in ", path)
  if (!("SE" %in% names(dat)))     stop("Missing column: SE in ", path)

  dat %>%
    transmute(
      Time = as.numeric(Time),
      region = as.character(region),
      Value = as.numeric(Value),
      SE = as.numeric(SE),
      condition = condition_label
    ) %>%
    filter(region %in% roi_order) %>%
    mutate(
      region = factor(region, levels = roi_order),
      condition = factor(condition, levels = c("Voluntary", "Involuntary")),
      ymin = Value - SE,
      ymax = Value + SE
    )
}

plot_mean_SE_tc <- function(dat,
                            pal_line,
                            pal_fill,
                            title_text,
                            ncol = 4,
                            ribbon_order = c("Involuntary", "Voluntary"),
                            ribbon_alpha = c(
                              "Voluntary" = 0.54,
                              "Involuntary" = 0.50
                            ),
                            border_lw = 0.20, base_size = 6,
                            line_width = .35, font_family = "Arial",
                            highlight_regions = character(),
                            highlight_border_color = "#9C5357", highlight_border_lw = .50) {

  dat <- dat %>%
    mutate(
      condition = factor(
        condition,
        levels = c("Voluntary", "Involuntary")
      )
    )

  dat_rib1 <- dat %>% filter(condition == ribbon_order[1])
  dat_rib2 <- dat %>% filter(condition == ribbon_order[2])

  dat_line_v <- dat %>% filter(condition == "Voluntary")
  dat_line_i <- dat %>% filter(condition == "Involuntary")
  border_data <- distinct(dat, region)
  stopifnot(!anyDuplicated(highlight_regions),
    all(highlight_regions %in% as.character(border_data$region)))

  p <- ggplot() +
    geom_ribbon(
      data = dat_rib1,
      aes(
        x = Time,
        ymin = ymin,
        ymax = ymax,
        fill = condition,
        group = condition
      ),
      alpha = unname(ribbon_alpha[ribbon_order[1]]),
      color = NA
    ) +
    geom_ribbon(
      data = dat_rib2,
      aes(
        x = Time,
        ymin = ymin,
        ymax = ymax,
        fill = condition,
        group = condition
      ),
      alpha = unname(ribbon_alpha[ribbon_order[2]]),
      color = NA
    ) +
    geom_line(
      data = dat_line_v,
      aes(
        x = Time,
        y = Value,
        color = condition,
        group = condition
      ),
      linewidth = line_width,
      lineend = "round"
    ) +
    geom_line(
      data = dat_line_i,
      aes(
        x = Time,
        y = Value,
        color = condition,
        group = condition
      ),
      linewidth = line_width,
      lineend = "round"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.20,
      color = "grey45"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.15,
      color = "grey80"
    ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      fill = NA,
      colour = "black",
      linewidth = border_lw
    ) +
    facet_wrap(
      ~ region,
      labeller=labeller(region=c(DLPFCR="DLPFC-R",PreSMAL="PreSMA-L",PreSMAR="PreSMA-R",SMAR="SMA-R",
        InsulaR="Insula-R",PrecentralR="PreCG-R",PostcentralL="PoCG-L",PostcentralR="PoCG-R",
        SPGL="SPG-L",PrecuneusR="Precuneus-R",AGR="AG-R",MTGR="MTG-R")),
      ncol = ncol,
      drop = FALSE,
      scales = "free_y"
    ) +
    scale_color_manual(
      values = pal_line,
      name = NULL, breaks=c("Voluntary","Involuntary")
    ) +
    scale_fill_manual(
      values = pal_fill,
      name = NULL, breaks=c("Voluntary","Involuntary")
    ) +
    scale_x_continuous(
      breaks = seq(-200, 1000, 200),
      limits = c(-200, 1000),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = y_breaks_1dp,
      labels = function(x) sprintf("%.1f", x),
      expand = expansion(mult = c(0.06, 0.06))
    ) +
    labs(
      title = title_text,
      x = "Time (ms)",
      y = "Amplitude (z-score)"
    ) +
    guides(fill="none") +
    theme_bw(base_size = base_size, base_family = font_family) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = base_size+1,
        face = "bold"
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        size = base_size+1,
        face = "plain",
        color = "black",
        margin = margin(b = 3, t = 1)
      ),
      axis.title = element_text(
        size = base_size+1,
        face = "bold"
      ),
      axis.text = element_text(
        size = base_size,
        color = "black"
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.15
      ),
      axis.ticks.length = unit(.7, "mm"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = base_size),
      legend.key.height = unit(0.25, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.key = element_rect(
        fill = "white",
        colour = NA
      ),
      panel.grid = element_blank(),
      panel.background = element_rect(
        fill = "white",
        color = NA
      ),
      panel.border = element_blank(),
      axis.line = element_blank(),
      panel.spacing = unit(0.9, "lines"),
      plot.margin = margin(3, 14, 3, 4)
    )

  if(length(highlight_regions)) {
    p <- p + geom_rect(
      data = filter(border_data, as.character(region) %in% highlight_regions),
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = NA, colour = highlight_border_color,
      linewidth = highlight_border_lw, show.legend = FALSE
    )
  }
  return(p)
}
