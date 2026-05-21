# ============================================================
# RMSF-only analysis for GROMACS 0–50 ns MD results
# Working directory should contain:
#   antibody_rank5/ or antibody_rank5_10ns/
#   nanobody_rank7/ or nanobody_rank7_10ns/
# Output folders:
#   nanobody/
#   antibody/
#   combined/
#   summary/
# ============================================================

library(ggplot2)
library(dplyr)
library(readr)

# ------------------------------------------------------------
# Detect input folders
# ------------------------------------------------------------

find_dir <- function(candidates) {
  for (d in candidates) {
    if (dir.exists(d)) return(d)
  }
  stop("Cannot find directory: ", paste(candidates, collapse = " or "))
}

systems <- list(
  antibody = find_dir(c("antibody_rank5", "antibody_rank5_10ns")),
  nanobody = find_dir(c("nanobody_rank7", "nanobody_rank7_10ns"))
)

dir.create("antibody", showWarnings = FALSE)
dir.create("nanobody", showWarnings = FALSE)
dir.create("combined", showWarnings = FALSE)
dir.create("summary", showWarnings = FALSE)

# ------------------------------------------------------------
# Colours and theme
# ------------------------------------------------------------

system_cols <- c(
  antibody = "#D55E00",
  nanobody = "#0072B2"
)

chain_cols <- c(
  antigen = "#009E73",
  binder = "#CC79A7",
  heavy = "#0072B2",
  light = "#E69F00",
  antigen_group = "#009E73"
)

clean_system_label <- function(x) {
  ifelse(x == "antibody", "IgG antibody", "Nanobody")
}

clean_chain_label <- function(chain) {
  dplyr::recode(
    chain,
    antigen = "antigen chain",
    binder = "binder chain",
    heavy = "heavy chain",
    light = "light chain",
    antigen_group = "antigen group",
    .default = chain
  )
}

theme_report <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.3),
      axis.line = element_line(colour = "black", linewidth = 0.45),
      axis.ticks = element_line(colour = "black", linewidth = 0.45),
      axis.text = element_text(size = 14, colour = "black"),
      axis.title = element_text(size = 16, colour = "black"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      legend.key.size = unit(0.55, "cm"),
      legend.background = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(8, 10, 8, 10)
    )
}

# ------------------------------------------------------------
# Read XVG
# ------------------------------------------------------------

read_xvg <- function(file) {
  if (!file.exists(file)) {
    message("Missing: ", file)
    return(NULL)
  }
  
  lines <- readLines(file, warn = FALSE)
  lines <- lines[!grepl("^[@#]", lines)]
  lines <- lines[nchar(trimws(lines)) > 0]
  
  if (length(lines) == 0) {
    message("No data in: ", file)
    return(NULL)
  }
  
  df <- read.table(text = lines, header = FALSE)
  
  if (ncol(df) < 2) {
    message("Not enough columns in: ", file)
    return(NULL)
  }
  
  df <- df[, 1:2]
  colnames(df) <- c("residue", "rmsf")
  return(df)
}

choose_file <- function(system_dir, candidates) {
  for (f in candidates) {
    full <- file.path(system_dir, f)
    if (file.exists(full)) return(f)
  }
  return(NULL)
}

summarise_rmsf <- function(df, system_name, chain_name, file_name) {
  if (is.null(df)) return(NULL)
  
  data.frame(
    system = system_name,
    chain = chain_name,
    file = file_name,
    n_points = nrow(df),
    mean_rmsf = mean(df$rmsf, na.rm = TRUE),
    sd_rmsf = sd(df$rmsf, na.rm = TRUE),
    min_rmsf = min(df$rmsf, na.rm = TRUE),
    max_rmsf = max(df$rmsf, na.rm = TRUE)
  )
}

# ------------------------------------------------------------
# RMSF file definitions
# ------------------------------------------------------------

all_rmsf_defs <- list(
  rmsf = c("rmsf_0_50ns.xvg", "rmsf_20ns.xvg", "rmsf_10ns.xvg")
)

chain_rmsf_defs <- list(
  nanobody = list(
    antigen = c("rmsf_0_50ns_Chain_A.xvg", "rmsf_20ns_antigen_chainA.xvg", "rmsf_10ns_antigen_chainA.xvg"),
    binder  = c("rmsf_0_50ns_Chain_B.xvg", "rmsf_20ns_nanobody_chainB.xvg", "rmsf_10ns_nanobody_chainB.xvg")
  ),
  antibody = list(
    antigen = c("rmsf_0_50ns_Chain_A.xvg", "rmsf_20ns_antigen_chainA.xvg", "rmsf_10ns_antigen_chainA.xvg"),
    heavy   = c("rmsf_0_50ns_Chain_B.xvg", "rmsf_20ns_heavy_chainB.xvg", "rmsf_10ns_heavy_chainB.xvg"),
    light   = c("rmsf_0_50ns_Chain_C.xvg", "rmsf_20ns_light_chainC.xvg", "rmsf_10ns_light_chainC.xvg"),
    binder  = c("rmsf_0_50ns_Binder.xvg"),
    antigen_group = c("rmsf_0_50ns_Antigen.xvg")
  )
)

# ------------------------------------------------------------
# Individual RMSF plot
# ------------------------------------------------------------

plot_rmsf_single <- function(df, system_name, chain_name = "all chains", shade_epitope = FALSE) {
  sys_label <- clean_system_label(system_name)
  
  if (chain_name == "all chains") {
    plot_col <- system_cols[[system_name]]
    y_lab <- paste0("RMSF of ", sys_label, " full complex (nm)")
  } else {
    plot_col <- chain_cols[[chain_name]]
    y_lab <- paste0("RMSF of ", sys_label, " ", clean_chain_label(chain_name), " (nm)")
  }
  
  p <- ggplot(df, aes(x = residue, y = rmsf)) +
    geom_line(linewidth = 0.75, colour = plot_col) +
    labs(
      x = "Residue index",
      y = y_lab
    ) +
    theme_report()
  
  if (shade_epitope) {
    p <- p +
      annotate(
        "rect",
        xmin = 242, xmax = 248,
        ymin = -Inf, ymax = Inf,
        fill = "#D55E00",
        alpha = 0.28
      ) +
      annotate(
        "rect",
        xmin = 282, xmax = 285,
        ymin = -Inf, ymax = Inf,
        fill = "#F0C000",
        alpha = 0.38
      ) +
      geom_line(linewidth = 0.75, colour = plot_col)
  }
  
  return(p)
}

summary_list <- list()

# ------------------------------------------------------------
# Individual all-chain RMSF and chain-specific RMSF
# ------------------------------------------------------------

for (sys in names(systems)) {
  sys_dir <- systems[[sys]]
  out_dir <- sys
  
  message("\nProcessing RMSF for ", sys, " from ", sys_dir)
  
  chosen <- choose_file(sys_dir, all_rmsf_defs$rmsf)
  if (!is.null(chosen)) {
    df <- read_xvg(file.path(sys_dir, chosen))
    if (!is.null(df)) {
      p <- plot_rmsf_single(df, sys, "all chains", shade_epitope = FALSE)
      out_file <- file.path(out_dir, "rmsf.png")
      ggsave(out_file, p, width = 7.2, height = 5.0, dpi = 400)
      message("Saved: ", out_file)
      
      summary_list[[paste(sys, "rmsf_all", sep = "_")]] <- summarise_rmsf(
        df, sys, "all chains", chosen
      )
    }
  }
  
  for (chain_name in names(chain_rmsf_defs[[sys]])) {
    chosen <- choose_file(sys_dir, chain_rmsf_defs[[sys]][[chain_name]])
    
    if (is.null(chosen)) {
      message("No chain RMSF found for ", sys, " / ", chain_name)
      next
    }
    
    df <- read_xvg(file.path(sys_dir, chosen))
    if (is.null(df)) next
    
    shade <- chain_name %in% c("antigen", "antigen_group")
    
    p <- plot_rmsf_single(df, sys, chain_name, shade_epitope = shade)
    
    out_file <- file.path(out_dir, paste0("rmsf_chain_", chain_name, ".png"))
    ggsave(out_file, p, width = 7.2, height = 5.0, dpi = 400)
    message("Saved: ", out_file)
    
    summary_list[[paste(sys, "rmsf_chain", chain_name, sep = "_")]] <- summarise_rmsf(
      df, sys, chain_name, chosen
    )
  }
}

# ------------------------------------------------------------
# Combined all-chain RMSF
# ------------------------------------------------------------

dfs_all <- list()

for (sys in names(systems)) {
  sys_dir <- systems[[sys]]
  chosen <- choose_file(sys_dir, all_rmsf_defs$rmsf)
  
  if (!is.null(chosen)) {
    df <- read_xvg(file.path(sys_dir, chosen))
    if (!is.null(df)) {
      df$system <- sys
      dfs_all[[sys]] <- df
    }
  }
}

if (length(dfs_all) > 0) {
  dat_all <- bind_rows(dfs_all)
  
  p_all <- ggplot(dat_all, aes(x = residue, y = rmsf, colour = system)) +
    geom_line(linewidth = 0.75) +
    scale_colour_manual(
      values = system_cols,
      labels = c(antibody = "IgG antibody full complex", nanobody = "Nanobody full complex")
    ) +
    labs(
      x = "Residue index",
      y = "RMSF of full complex (nm)"
    ) +
    theme_report()
  
  ggsave(
    file.path("combined", "combined_rmsf.png"),
    p_all,
    width = 7.4,
    height = 5.2,
    dpi = 400
  )
  
  message("Saved: combined/combined_rmsf.png")
}

# ------------------------------------------------------------
# Combined antigen-chain RMSF with epitope shading
# Overlay antibody and nanobody antigen chains in one plot
# ------------------------------------------------------------

dfs_antigen <- list()

for (sys in names(systems)) {
  sys_dir <- systems[[sys]]
  chosen <- choose_file(sys_dir, chain_rmsf_defs[[sys]]$antigen)
  
  if (!is.null(chosen)) {
    df <- read_xvg(file.path(sys_dir, chosen))
    if (!is.null(df)) {
      df$system <- sys
      dfs_antigen[[sys]] <- df
    }
  }
}

if (length(dfs_antigen) > 0) {
  dat_antigen <- bind_rows(dfs_antigen)
  
  p_antigen <- ggplot(dat_antigen, aes(x = residue, y = rmsf, colour = system)) +
    annotate(
      "rect",
      xmin = 242, xmax = 248,
      ymin = -Inf, ymax = Inf,
      fill = "#D55E00",
      alpha = 0.28
    ) +
    annotate(
      "rect",
      xmin = 282, xmax = 285,
      ymin = -Inf, ymax = Inf,
      fill = "#F0C000",
      alpha = 0.38
    ) +
    geom_vline(xintercept = c(242, 248), colour = "#D55E00", linewidth = 0.35, alpha = 0.75) +
    geom_vline(xintercept = c(282, 285), colour = "#B8860B", linewidth = 0.35, alpha = 0.75) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(
      values = system_cols,
      labels = c(
        antibody = "IgG antibody CD276 antigen chain",
        nanobody = "Nanobody CD276 antigen chain"
      )
    ) +
    labs(
      x = "CD276 residue index",
      y = "RMSF of CD276 antigen chain (nm)"
    ) +
    theme_report() +
    theme(
      legend.position = "bottom"
    )
  
  ggsave(
    file.path("combined", "combined_antigen_rmsf_epitope.png"),
    p_antigen,
    width = 7.4,
    height = 4.8,
    dpi = 400
  )
  
  message("Saved: combined/combined_antigen_rmsf_epitope.png")
}

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary_df <- bind_rows(summary_list)
write_csv(summary_df, file.path("summary", "RMSF_metric_summary.csv"))

message("\nRMSF-only plots and summaries finished.")
message("Key report-ready figure:")
message(" - combined/combined_antigen_rmsf_epitope.png")