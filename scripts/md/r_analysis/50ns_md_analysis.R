# ============================================================
# Plot GROMACS 0–50 ns MD analysis results, excluding RMSF
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
# Theme and colours
# ------------------------------------------------------------

system_cols <- c(
  antibody = "#D55E00",
  nanobody = "#0072B2"
)

clean_system_label <- function(x) {
  ifelse(x == "antibody", "IgG antibody", "Nanobody")
}

theme_report <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.3),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.background = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA)
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
  colnames(df) <- c("x", "y")
  return(df)
}

choose_file <- function(system_dir, candidates) {
  for (f in candidates) {
    full <- file.path(system_dir, f)
    if (file.exists(full)) return(f)
  }
  return(NULL)
}

prepare_axis <- function(df, metric) {
  if (is.null(df)) return(NULL)
  
  if (metric == "residue") {
    df$residue <- df$x
    return(df)
  }
  
  if (metric == "potential") {
    df$step <- df$x
    return(df)
  }
  
  if (metric == "rmsd") {
    df$time_ns <- df$x
    return(df)
  }
  
  if (max(df$x, na.rm = TRUE) > 100) {
    df$time_ns <- df$x / 1000
  } else {
    df$time_ns <- df$x
  }
  
  return(df)
}

# ------------------------------------------------------------
# Metric definitions, excluding RMSF
# ------------------------------------------------------------

metric_defs <- list(
  rmsd = list(
    candidates = c("rmsd_0_50ns.xvg", "rmsd_20ns.xvg", "rmsd_10ns.xvg", "rmsd.xvg"),
    ylab = "Backbone RMSD (nm)",
    metric = "rmsd"
  ),
  rmsd_xtal = list(
    candidates = c("rmsd_0_50ns_xtal.xvg", "rmsd_20ns_xtal.xvg", "rmsd_10ns_xtal.xvg", "rmsd_xtal.xvg"),
    ylab = "Backbone RMSD vs minimised structure (nm)",
    metric = "rmsd"
  ),
  rg = list(
    candidates = c("gyrate_0_50ns.xvg", "gyrate_20ns.xvg", "gyrate_10ns.xvg", "gyrate.xvg"),
    ylab = "Radius of gyration (nm)",
    metric = "gyrate"
  ),
  sasa = list(
    candidates = c("sasa_0_50ns.xvg", "sasa_20ns.xvg", "sasa_10ns.xvg"),
    ylab = expression(SASA~(nm^2)),
    metric = "sasa"
  ),
  sasa_residue = list(
    candidates = c("sasa_residue_0_50ns.xvg", "sasa_residue_20ns.xvg", "sasa_residue_10ns.xvg"),
    ylab = expression("Residue-level SASA"~(nm^2)),
    metric = "residue"
  ),
  mindist = list(
    candidates = c("mindist_0_50ns.xvg", "mindist_20ns.xvg", "mindist_10ns.xvg"),
    ylab = "Minimum binder-antigen distance (nm)",
    metric = "mindist"
  ),
  contacts = list(
    candidates = c("contacts_0_50ns.xvg", "contacts_20ns.xvg", "contacts_10ns.xvg"),
    ylab = "Number of interface contacts",
    metric = "contacts"
  ),
  hbonds = list(
    candidates = c("hbonds_0_50ns.xvg", "hbonds_20ns.xvg", "hbonds_10ns.xvg"),
    ylab = "Number of interface hydrogen bonds",
    metric = "hbonds"
  ),
  potential = list(
    candidates = c("potential.xvg"),
    ylab = "Potential energy (kJ/mol)",
    metric = "potential"
  ),
  temperature = list(
    candidates = c("temperature.xvg"),
    ylab = "Temperature (K)",
    metric = "temperature"
  ),
  pressure = list(
    candidates = c("pressure.xvg"),
    ylab = "Pressure (bar)",
    metric = "pressure"
  ),
  density = list(
    candidates = c("density.xvg"),
    ylab = expression(Density~(kg~m^{-3})),
    metric = "density"
  )
)

# ------------------------------------------------------------
# Plot functions
# ------------------------------------------------------------

plot_single <- function(df, system_name, metric_name, ylab) {
  df <- prepare_axis(df, metric_name)
  plot_col <- system_cols[[system_name]]
  system_label <- clean_system_label(system_name)
  
  if (metric_name == "residue") {
    p <- ggplot(df, aes(x = residue, y = y)) +
      geom_line(linewidth = 0.7, colour = plot_col) +
      labs(
        x = "Residue index",
        y = paste0(system_label, ": ", ylab)
      ) +
      theme_report()
  } else if (metric_name == "potential") {
    p <- ggplot(df, aes(x = step, y = y)) +
      geom_line(linewidth = 0.7, colour = plot_col) +
      labs(
        x = "Energy minimisation step",
        y = paste0(system_label, ": ", ylab)
      ) +
      theme_report()
  } else {
    p <- ggplot(df, aes(x = time_ns, y = y)) +
      geom_line(linewidth = 0.8, colour = plot_col) +
      labs(
        x = "Time (ns)",
        y = paste0(system_label, ": ", ylab)
      ) +
      theme_report()
  }
  
  return(p)
}

summarise_metric <- function(df, system_name, metric_name, file_name) {
  if (is.null(df)) return(NULL)
  
  y <- df$y
  n <- length(y)
  last20_start <- max(1, floor(n * 0.8))
  last20 <- y[last20_start:n]
  
  data.frame(
    system = system_name,
    metric = metric_name,
    file = file_name,
    n_points = n,
    mean_all = mean(y, na.rm = TRUE),
    sd_all = sd(y, na.rm = TRUE),
    min_all = min(y, na.rm = TRUE),
    max_all = max(y, na.rm = TRUE),
    last_value = tail(y, 1),
    mean_last20 = mean(last20, na.rm = TRUE),
    sd_last20 = sd(last20, na.rm = TRUE)
  )
}

# ------------------------------------------------------------
# Individual plots
# ------------------------------------------------------------

summary_list <- list()

for (sys in names(systems)) {
  sys_dir <- systems[[sys]]
  out_dir <- sys
  
  message("\nProcessing ", sys, " from ", sys_dir)
  
  for (metric_name in names(metric_defs)) {
    def <- metric_defs[[metric_name]]
    chosen <- choose_file(sys_dir, def$candidates)
    
    if (is.null(chosen)) {
      message("No file found for ", sys, " / ", metric_name)
      next
    }
    
    df <- read_xvg(file.path(sys_dir, chosen))
    if (is.null(df)) next
    
    p <- plot_single(df, sys, def$metric, def$ylab)
    
    out_file <- file.path(out_dir, paste0(metric_name, ".png"))
    ggsave(out_file, p, width = 7.2, height = 5.0, dpi = 400)
    message("Saved: ", out_file)
    
    summary_list[[paste(sys, metric_name, sep = "_")]] <- summarise_metric(
      df, sys, metric_name, chosen
    )
  }
}

# ------------------------------------------------------------
# Combined plots
# ------------------------------------------------------------

for (metric_name in names(metric_defs)) {
  def <- metric_defs[[metric_name]]
  dfs <- list()
  
  for (sys in names(systems)) {
    sys_dir <- systems[[sys]]
    chosen <- choose_file(sys_dir, def$candidates)
    
    if (!is.null(chosen)) {
      df <- read_xvg(file.path(sys_dir, chosen))
      if (!is.null(df)) {
        df <- prepare_axis(df, def$metric)
        df$system <- sys
        dfs[[sys]] <- df
      }
    }
  }
  
  if (length(dfs) == 0) next
  
  dat <- bind_rows(dfs)
  
  if (def$metric == "residue") {
    p <- ggplot(dat, aes(x = residue, y = y, colour = system)) +
      geom_line(linewidth = 0.75) +
      scale_colour_manual(
        values = system_cols,
        labels = c(antibody = "IgG antibody", nanobody = "Nanobody")
      ) +
      labs(
        x = "Residue index",
        y = def$ylab
      ) +
      theme_report()
  } else if (def$metric == "potential") {
    p <- ggplot(dat, aes(x = step, y = y, colour = system)) +
      geom_line(linewidth = 0.75) +
      scale_colour_manual(
        values = system_cols,
        labels = c(antibody = "IgG antibody", nanobody = "Nanobody")
      ) +
      labs(
        x = "Energy minimisation step",
        y = def$ylab
      ) +
      theme_report()
  } else {
    p <- ggplot(dat, aes(x = time_ns, y = y, colour = system)) +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(
        values = system_cols,
        labels = c(antibody = "IgG antibody", nanobody = "Nanobody")
      ) +
      labs(
        x = "Time (ns)",
        y = def$ylab
      ) +
      theme_report()
  }
  
  out_file <- file.path("combined", paste0("combined_", metric_name, ".png"))
  ggsave(out_file, p, width = 7.4, height = 5.2, dpi = 400)
  message("Saved: ", out_file)
}

# ------------------------------------------------------------
# Summary table
# ------------------------------------------------------------

summary_df <- bind_rows(summary_list)
write_csv(summary_df, file.path("summary", "MD_metric_summary_no_rmsf.csv"))

message("\nNon-RMSF plots and summaries finished.")