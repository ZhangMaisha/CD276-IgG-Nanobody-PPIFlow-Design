library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

# =========================
# Read data
# =========================

df <- read_excel("PRODIGY_manual_summary.xlsx")

# =========================
# Clean labels
# =========================

plot_df <- df %>%
  mutate(
    candidate = as.character(candidate),
    
    system = case_when(
      str_detect(candidate, regex("IgG|antibody", ignore_case = TRUE)) ~ "IgG antibody",
      str_detect(candidate, regex("NB|nanobody", ignore_case = TRUE)) ~ "Nanobody",
      TRUE ~ candidate
    ),
    
    stage = case_when(
      str_detect(candidate, regex("recheck|AF3", ignore_case = TRUE)) ~ "After AF3",
      TRUE ~ "Before AF3"
    ),
    
    stage = factor(stage, levels = c("Before AF3", "After AF3")),
    system = factor(system, levels = c("IgG antibody", "Nanobody")),
    dg_label = sprintf("%.1f", dg_kcal_mol)
  )

# =========================
# Plot
# =========================

LINE_WIDTH <- 0.45

p <- ggplot(plot_df, aes(x = stage, y = dg_kcal_mol, group = system, colour = system)) +
  
  geom_hline(
    yintercept = -9,
    linetype = "dashed",
    linewidth = LINE_WIDTH,
    colour = "grey35"
  ) +
  
  geom_line(linewidth = LINE_WIDTH) +
  geom_point(size = 2.8) +
  
  geom_text(
    aes(label = dg_label),
    vjust = -1,
    size = 3.4,
    show.legend = FALSE
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = -8.85,
    label = expression(Delta * "G cutoff = -9 kcal/mol"),
    size = 3.1,
    colour = "grey25"
  ) +
  
  scale_colour_manual(
    values = c(
      "IgG antibody" = "#D55E00",
      "Nanobody" = "#0072B2"
    )
  ) +
  
  scale_y_continuous(
    name = expression("Prodigy-predicted " * Delta * "G (kcal/mol)"),
    expand = expansion(mult = c(0.08, 0.16))
  ) +
  
  labs(
    x = NULL,
    colour = NULL
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    axis.text = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 12, 10, 12)
  )

print(p)

ggsave(
  filename = "Figure2A_Prodigy_DeltaG.png",
  plot = p,
  width = 5.2,
  height = 3.6,
  dpi = 600
)