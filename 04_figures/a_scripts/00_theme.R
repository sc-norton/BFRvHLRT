# =============================================================================
# 00_theme.R
# Shared color palette and ggplot theme for CRC_pro figures.
#
# Source this near the top of every figure script instead of redefining
# `colors` / `theme_cr()` locally:
#
#   source("04_figures/a_scripts/00_theme.R")
#
# Canonical defaults are base_size = 12, title_hjust = 0.5 (left over from
# the majority of existing figure scripts). Override either argument at a
# given call site (e.g. theme_cr(base_size = 9, title_hjust = 0)) for dense
# panels or left-aligned titles without touching this file.
#
# Figure-specific colors (e.g. a `tissue` or `trajectory` palette used by
# only one script) should stay local to that script, not live here.
# =============================================================================

library(ggplot2)

colors <- list(
  red          = "#C0392B", red_light    = "#E6A19A", red_dark     = "#7B241C",
  blue         = "#2471A3", blue_light   = "#A9CCE3", blue_dark    = "#154360",
  green        = "#4E8B6F", green_light  = "#A9C9B8", green_dark   = "#2E5E4E",
  pink         = "#D98C9C", pink_light   = "#F2C6CF", pink_dark    = "#A65C6B",
  purple       = "#7B5EA7", purple_light = "#C5B4DD", purple_dark  = "#4A3A73",
  yellow       = "#D4A017", yellow_light = "#F1D77A", yellow_dark  = "#9A6F0A",
  orange       = "#C4693A", orange_light = "#E5A07B", orange_dark  = "#8A4321",
  aqua         = "#5DAFA3", aqua_light   = "#A9D6CF", aqua_dark    = "#3A7F76",
  brown        = "#8C5A3C", brown_light  = "#C7A189", brown_dark   = "#5A3A28",
  red_pale     = "#FADBD8", blue_pale    = "#D6EAF8",
  gray_light   = "#F5F5F5", gray         = "#9E9E9E", gray_dark    = "#4A4A4A",
  black        = "#000000", white        = "#FFFFFF",
  auburn       = c(navy = "#0D234A", auburn = "#E96A2C")
)

theme_cr <- function(base_size = 12, title_hjust = 0.5) {
  theme_classic(base_size = base_size) +
    theme(
      text            = element_text(face = "bold"),
      axis.title      = element_text(face = "bold", size = base_size),
      axis.text       = element_text(face = "bold", size = base_size),
      strip.text      = element_text(face = "bold", size = base_size),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold", size = base_size),
      legend.text     = element_text(size = base_size),
      panel.spacing   = unit(1.5, "lines"),
      panel.border    = element_blank(),
      plot.title      = element_text(face = "bold", size = base_size + 1, hjust = title_hjust)
    )
}
