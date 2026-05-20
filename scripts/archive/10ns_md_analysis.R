# ============================================================
# Plot GROMACS MD results for CD276 antibody/nanobody project
# Output folders:
#   nanobody/
#   antibody/
#   combined/
# ============================================================

# install.packages(c("ggplot2", "dplyr", "readr", "patchwork"))

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# ------------------------------------------------------------
# 1. Working directory
# ------------------------------------------------------------

setwd("~/Desktop/CMML3/ICA2/MD_analysis_updated")

systems <- list(
  nanobody = "nanobody_rank7_10ns",
  antibody = "antibody_rank5_10ns"
)

dir.create("nanobody", showWarnings = FALSE)
dir.create("antibody", showWarnings = FALSE)
dir.create("combined", showWarnings = FALSE)
dir.create("summary", showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Read XVG file
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
    message("No data: ", file)
    return(NULL)
  }
  
  df <- read.table(text = lines, header = FALSE)
  
  if (ncol(df) < 2) {
    message("Not enough columns: ", file)
    return(NULL)
  }
  
  df <- df[, 1:2]
  colnames(df) <- c("x", "y")
  return(df)
}

# ------------------------------------------------------------
# 3. Time conversion
# ------------------------------------------------------------

prepare_x_axis <- function(df, metric) {
  if (is.null(df)) return(NULL)
  
  if (grepl("rmsf", metric, ignore.case = TRUE)) {
    df$residue <- df$x
    return(df)
  }
  
  if (grepl("potential", metric, ignore.case = TRUE)) {
    df$step <- df$x
    return(df)
  }
  
  # RMSD generated with -tu ns usually already uses ns
  if (grepl("rmsd", metric, ignore.case = TRUE)) {
    df$time_ns <- df$x
    return(df)
  }
  
  # Other time series may be in ps; convert if values are large
  if (max(df$x, na.rm = TRUE) > 100) {
    df$time_ns <- df$x / 1000
  } else {
    df$time_ns <- df$x
  }
  
  return(df)
}

# ------------------------------------------------------------
# 4. Metric definitions
# ------------------------------------------------------------

metrics <- list(
  list(
    file = "potential.xvg",
    metric = "potential",
    title = "Potential energy",
    xlab = "Energy minimization step",
    ylab = "Potential energy (kJ/mol)"
  ),
  list(
    file = "temperature.xvg",
    metric = "temperature",
    title = "Temperature during NVT",
    xlab = "Time (ns)",
    ylab = "Temperature (K)"
  ),
  list(
    file = "pressure.xvg",
    metric = "pressure",
    title = "Pressure during NPT",
    xlab = "Time (ns)",
    ylab = "Pressure (bar)"
  ),
  list(
    file = "density.xvg",
    metric = "density",
    title = "Density during NPT",
    xlab = "Time (ns)",
    ylab = expression(Density~(kg~m^{-3}))
  ),
  list(
    file = "rmsd_10ns.xvg",
    metric = "rmsd_10ns",
    title = "Backbone RMSD, 0-10 ns",
    xlab = "Time (ns)",
    ylab = "Backbone RMSD (nm)"
  ),
  list(
    file = "rmsd_10ns_xtal.xvg",
    metric = "rmsd_10ns_xtal",
    title = "Backbone RMSD vs minimized structure, 0-10 ns",
    xlab = "Time (ns)",
    ylab = "Backbone RMSD (nm)"
  ),
  list(
    file = "gyrate_10ns.xvg",
    metric = "gyrate_10ns",
    title = "Radius of gyration, 0-10 ns",
    xlab = "Time (ns)",
    ylab = "Rg (nm)"
  ),
  list(
    file = "rmsf_10ns.xvg",
    metric = "rmsf_10ns",
    title = "RMSF, all chains, 0-10 ns",
    xlab = "Residue index",
    ylab = "RMSF (nm)"
  ),
  list(
    file = "sasa_10ns.xvg",
    metric = "sasa_10ns",
    title = "SASA, 0-10 ns",
    xlab = "Time (ns)",
    ylab = expression(SASA~(nm^2))
  ),
  list(
    file = "mindist_10ns.xvg",
    metric = "mindist_10ns",
    title = "Minimum binder-antigen distance, 0-10 ns",
    xlab = "Time (ns)",
    ylab = "Minimum distance (nm)"
  ),
  list(
    file = "contacts_10ns.xvg",
    metric = "contacts_10ns",
    title = "Interface contacts, 0-10 ns",
    xlab = "Time (ns)",
    ylab = "Number of contacts"
  ),
  list(
    file = "hbonds_10ns.xvg",
    metric = "hbonds_10ns",
    title = "Interface hydrogen bonds, 0-10 ns",
    xlab = "Time (ns)",
    ylab = "Number of hydrogen bonds"
  )
)

# chain-specific RMSF files, if available
chain_rmsf_files <- list(
  nanobody = c(
    "rmsf_10ns_antigen_chainA.xvg",
    "rmsf_10ns_nanobody_chainB.xvg"
  ),
  antibody = c(
    "rmsf_10ns_antigen_chainA.xvg",
    "rmsf_10ns_heavy_chainB.xvg",
    "rmsf_10ns_light_chainC.xvg"
  )
)

# ------------------------------------------------------------
# 5. Plot one file
# ------------------------------------------------------------

plot_one <- function(system_name, system_dir, file_name, metric, title, xlab, ylab, output_dir) {
  file <- file.path(system_dir, file_name)
  df <- read_xvg(file)
  if (is.null(df)) return(NULL)
  
  df <- prepare_x_axis(df, metric)
  
  if (grepl("rmsf", metric, ignore.case = TRUE)) {
    p <- ggplot(df, aes(x = residue, y = y)) +
      geom_line(linewidth = 0.6, colour = "black") +
      theme_bw(base_size = 14) +
      labs(
        x = "Residue index",
        y = ylab,
        title = paste0(system_name, ": ", title)
      )
  } else if (grepl("potential", metric, ignore.case = TRUE)) {
    p <- ggplot(df, aes(x = step, y = y)) +
      geom_line(linewidth = 0.6, colour = "black") +
      theme_bw(base_size = 14) +
      labs(
        x = xlab,
        y = ylab,
        title = paste0(system_name, ": ", title)
      )
  } else {
    p <- ggplot(df, aes(x = time_ns, y = y)) +
      geom_line(linewidth = 0.6, colour = "black") +
      theme_bw(base_size = 14) +
      labs(
        x = "Time (ns)",
        y = ylab,
        title = paste0(system_name, ": ", title)
      )
  }
  
  out_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(file_name), ".png"))
  ggsave(out_file, p, width = 7, height = 5, dpi = 300)
  message("Saved: ", out_file)
  
  return(df)
}

# ------------------------------------------------------------
# 6. Generate individual plots and summary
# ------------------------------------------------------------

summary_rows <- list()

summarise_metric <- function(df, system_name, file_name, metric) {
  if (is.null(df)) return(NULL)
  
  y <- df$y
  n <- length(y)
  last20_start <- max(1, floor(n * 0.8))
  last20 <- y[last20_start:n]
  
  data.frame(
    system = system_name,
    file = file_name,
    metric = metric,
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

for (system_name in names(systems)) {
  system_dir <- systems[[system_name]]
  output_dir <- system_name
  
  message("\nProcessing ", system_name)
  
  for (m in metrics) {
    df <- plot_one(
      system_name = system_name,
      system_dir = system_dir,
      file_name = m$file,
      metric = m$metric,
      title = m$title,
      xlab = m$xlab,
      ylab = m$ylab,
      output_dir = output_dir
    )
    
    s <- summarise_metric(df, system_name, m$file, m$metric)
    if (!is.null(s)) {
      summary_rows[[paste(system_name, m$file, sep = "_")]] <- s
    }
  }
  
  # chain-specific RMSF
  for (f in chain_rmsf_files[[system_name]]) {
    if (file.exists(file.path(system_dir, f))) {
      df <- plot_one(
        system_name = system_name,
        system_dir = system_dir,
        file_name = f,
        metric = "rmsf_chain_specific",
        title = paste0("Chain-specific RMSF: ", tools::file_path_sans_ext(f)),
        xlab = "Residue index",
        ylab = "RMSF (nm)",
        output_dir = output_dir
      )
      
      s <- summarise_metric(df, system_name, f, "rmsf_chain_specific")
      if (!is.null(s)) {
        summary_rows[[paste(system_name, f, sep = "_")]] <- s
      }
    }
  }
}

summary_df <- bind_rows(summary_rows)
write_csv(summary_df, file.path("summary", "MD_metric_summary.csv"))

# ------------------------------------------------------------
# 7. Combined plots
# ------------------------------------------------------------

make_combined_plot <- function(file_name, metric, title, ylab, x_type = "time") {
  data_list <- list()
  
  for (system_name in names(systems)) {
    file <- file.path(systems[[system_name]], file_name)
    df <- read_xvg(file)
    
    if (!is.null(df)) {
      df <- prepare_x_axis(df, metric)
      df$system <- system_name
      data_list[[system_name]] <- df
    }
  }
  
  if (length(data_list) == 0) {
    message("No data for combined plot: ", file_name)
    return(NULL)
  }
  
  dat <- bind_rows(data_list)
  
  if (x_type == "residue") {
    p <- ggplot(dat, aes(x = residue, y = y, colour = system)) +
      geom_line(linewidth = 0.6) +
      theme_bw(base_size = 14) +
      labs(
        x = "Residue index",
        y = ylab,
        title = title
      )
  } else if (x_type == "step") {
    p <- ggplot(dat, aes(x = step, y = y, colour = system)) +
      geom_line(linewidth = 0.6) +
      theme_bw(base_size = 14) +
      labs(
        x = "Energy minimization step",
        y = ylab,
        title = title
      )
  } else {
    p <- ggplot(dat, aes(x = time_ns, y = y, colour = system)) +
      geom_line(linewidth = 0.6) +
      theme_bw(base_size = 14) +
      labs(
        x = "Time (ns)",
        y = ylab,
        title = title
      )
  }
  
  out_file <- file.path("combined", paste0("combined_", tools::file_path_sans_ext(file_name), ".png"))
  ggsave(out_file, p, width = 7, height = 5, dpi = 300)
  message("Saved: ", out_file)
  
  return(p)
}

make_combined_plot("rmsd_10ns.xvg", "rmsd_10ns", "Backbone RMSD, 0-10 ns", "Backbone RMSD (nm)")
make_combined_plot("gyrate_10ns.xvg", "gyrate_10ns", "Radius of gyration, 0-10 ns", "Rg (nm)")
make_combined_plot("sasa_10ns.xvg", "sasa_10ns", "SASA, 0-10 ns", expression(SASA~(nm^2)))
make_combined_plot("mindist_10ns.xvg", "mindist_10ns", "Minimum binder-antigen distance, 0-10 ns", "Minimum distance (nm)")
make_combined_plot("contacts_10ns.xvg", "contacts_10ns", "Interface contacts, 0-10 ns", "Number of contacts")
make_combined_plot("hbonds_10ns.xvg", "hbonds_10ns", "Interface hydrogen bonds, 0-10 ns", "Number of hydrogen bonds")
make_combined_plot("rmsf_10ns.xvg", "rmsf_10ns", "RMSF, all chains, 0-10 ns", "RMSF (nm)", x_type = "residue")
make_combined_plot("potential.xvg", "potential", "Potential energy minimization", "Potential energy (kJ/mol)", x_type = "step")
make_combined_plot("temperature.xvg", "temperature", "Temperature during NVT", "Temperature (K)")
make_combined_plot("pressure.xvg", "pressure", "Pressure during NPT", "Pressure (bar)")
make_combined_plot("density.xvg", "density", "Density during NPT", expression(Density~(kg~m^{-3})))

# ------------------------------------------------------------
# 8. Create a short interpretation template
# ------------------------------------------------------------

interpretation_text <- "
Recommended report interpretation:

1. RMSD:
A stable complex should show an initial increase followed by a plateau. Continuous increase may indicate structural drift.

2. Rg:
A stable complex should maintain a relatively constant radius of gyration. A continuous increase suggests loosening or unfolding.

3. Interface contacts:
Stable antigen-binder association should maintain interface contacts over time. A drop toward zero suggests dissociation.

4. Hydrogen bonds:
Hydrogen bonds can fluctuate, but persistent or recurrent interfacial H-bonds support interface stability.

5. RMSF:
Framework regions are expected to be less flexible. CDR loops may show higher RMSF, but excessive peaks combined with loss of contacts/hydrogen bonds suggest unstable binding loops.

6. SASA:
Stable SASA suggests no major exposure or unfolding. A large continuous increase can indicate structural loosening.

7. Temperature / pressure / density:
These are primarily simulation quality-control metrics. Temperature should remain near 300 K and density should stabilise near water-like density.
"

writeLines(interpretation_text, con = file.path("summary", "MD_interpretation_notes.txt"))

message("\nAll done.")
message("Outputs:")
message(" - nanobody/")
message(" - antibody/")
message(" - combined/")
message(" - summary/MD_metric_summary.csv")
message(" - summary/MD_interpretation_notes.txt")