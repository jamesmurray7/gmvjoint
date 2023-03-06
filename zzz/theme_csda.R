theme_csda <- function(base_family = "Arial", base_size = 12){
  theme_light(base_family = base_family, base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 8.5, colour = "black"),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8),
      legend.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

theme_csda_bw <- function(base_family = "Arial", base_size = 12){
  theme_bw(base_family = base_family, base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = base_size * 0.9, colour = "black"),
      axis.text = element_text(size = base_size * 0.66),
      legend.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

.nice.orange <- "#FEB24C"