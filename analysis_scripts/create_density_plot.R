library(ggplot2)
library(readr)

angle_data <- read_tsv("outputs/measurement_angles.tsv")

angle_data[["size_class"]] <- factor(
  angle_data[["size_class"]],
  levels = c("Tiny", "Small", "Intermediate", "Large", "Bulky")
)
angle_data[["angle"]] <- ((angle_data[["angle"]] + 180) %% 360) - 180

size_class_colors <- c(
  "Tiny"         = "#4299e1",
  "Small"        = "#48bb78",
  "Intermediate" = "#f6ad55",
  "Large"        = "#f56565",
  "Bulky"        = "#c53030"
)

n <- nrow(angle_data)

orientation_plot <- angle_data |>
  ggplot(aes(x = angle, color = size_class)) +
  geom_density(size = 1.2, adjust = 1.2) +
  scale_color_manual(values = size_class_colors) +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, by = 50)) +
  labs(
    title = sprintf("Tripeptide (XRX) in Helix (n = %d)", n),
    x = "Angle between adjacent C-alpha → Centroid vectors [°]",
    y = "Norm. Freq. [A.U.]",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "#bdbdbd", color = NA),
    plot.background  = element_rect(fill = "#bdbdbd", color = NA),
    panel.grid.major.x = element_line(color = "#dddddd", linetype = "dotted"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position   = c(0.2, 0.8),
    legend.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("outputs/orientation_density.png", plot = orientation_plot, width = 10, height = 6, dpi = 300)
print(orientation_plot)
