################################################################################
#  make_tables_figures.R
#
#  Reads the summary files produced by sim_01 ... sim_08 and writes the
#  tables (LaTeX) and figures (PDF) used in Section 6 and in the Supplementary Material of the
#  paper into ../output/tables/ and ../output/figures/.
#
#    tab_size.tex         compact size table (main text)
#    tab_size_full.tex    full size table, all settings (appendix)
#    tab_misclass.tex     power under symmetric typing errors
#    tab_confounding.tex  rejection rates in the confounding scenario
#    fig_power.pdf        power curves, pair vs centroid statistic
#    (v0.5) tab_confounding_scenarios.tex (sim_04b), tab_segregation_cv.tex
#           (sim_05b), tab_robustness.tex, tab_robustness_full.tex,
#           tab_robustness_size.tex (sim_08)
#
#  Base R only. Run after the simulation scripts, from the R/ folder:
#     source("make_tables_figures.R")
################################################################################

source("plot_style.R")
dir_tab <- file.path("..", "output", "tables")
dir_fig <- file.path("..", "output", "figures")
dir.create(dir_tab, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

fmt <- function(x, digits = 3) {
  out <- formatC(x, format = "f", digits = digits)
  out[is.na(x)] <- "--"
  return(out)
}

write_lines <- function(lines, file) {
  writeLines(lines, file)
  cat("written:", file, "\n")
}

process_name <- c(hpp = "Homogeneous Poisson", ipp = "Inhomogeneous Poisson",
                  thomas = "Thomas cluster")
label_name <- c(oracle = "true types", est_easy = "estimated, SNR 1",
                est_hard = "estimated, SNR 0.05",
                est_very_hard = "estimated, SNR 0.03")

## =============================================================================
## Size tables (sim_01)
## =============================================================================
f1 <- file.path("results", "sim_01_size_summary.csv")
if (file.exists(f1)) {
  s1 <- read.csv(f1)
  m <- s1$n_rep[1]
  se <- sqrt(0.05 * 0.95 / m)

  get_value <- function(proc, n, lab, w, column) {
    sel <- s1$process == proc & s1$n == n & s1$labels == lab & s1$weights == w
    return(s1[sel, column])
  }

  ## ---- compact table (main text) ---------------------------------------------
  ## rows: process x n; columns: global pair test with unit and inverse
  ## weights for the three typings, centroid test with true types.
  lines <- c("\\begin{tabular}{llrrrrrrr}",
             "\\toprule",
             " & & \\multicolumn{3}{c}{Pair, unit weights} & \\multicolumn{3}{c}{Pair, inverse weights} & Centroid \\\\",
             "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-9}",
             "Process & $n$ & true & SNR 1 & SNR 0.05 & true & SNR 1 & SNR 0.05 & true \\\\",
             "\\midrule")
  for (proc in c("hpp", "ipp", "thomas")) {
    for (n in sort(unique(s1$n))) {
      row <- paste0(process_name[[proc]], " & ", n)
      for (w in c("none", "inverse_intensity")) {
        for (lab in c("oracle", "est_easy", "est_hard")) {
          row <- paste0(row, " & ", fmt(get_value(proc, n, lab, w, "size_global")))
        }
      }
      row <- paste0(row, " & ",
                    fmt(get_value(proc, n, "oracle", "none", "size_centroid_bonf")),
                    " \\\\")
      lines <- c(lines, row)
    }
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}",
             paste0("% Monte Carlo s.e. at 0.05 with ", m, " replicates: ",
                    fmt(se, 4)))
  write_lines(lines, file.path(dir_tab, "tab_size.tex"))

  ## ---- full table (appendix) ------------------------------------------------
  s1 <- s1[order(s1$process, s1$n, s1$labels, s1$weights), ]
  lines <- c("\\begin{tabular}{llllrrrr}",
             "\\toprule",
             "Process & $n$ & Labels & Weights & ARI & Global & Type 1 & Centroid \\\\",
             "\\midrule")
  for (r in 1:nrow(s1)) {
    lines <- c(lines, paste0(
      process_name[[s1$process[r]]], " & ", s1$n[r], " & ",
      label_name[[s1$labels[r]]], " & ",
      ifelse(s1$weights[r] == "none", "unit", "inverse"), " & ",
      fmt(s1$mean_ari[r], 2), " & ",
      fmt(s1$size_global[r]), " & ",
      fmt(s1$size_type1[r]), " & ",
      fmt(s1$size_centroid_bonf[r]), " \\\\"))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(dir_tab, "tab_size_full.tex"))
}

## =============================================================================
## Figure: power curves (sim_02)
## =============================================================================
f2 <- file.path("results", "sim_02_power_summary.csv")
if (file.exists(f2)) {
  s2 <- read.csv(f2)
  patterns <- c("one_bump", "two_bumps", "band")
  titles <- c("One hotspot", "Two hotspots", "Elongated band")
  modes <- c("oracle", "est_easy", "est_hard", "est_very_hard")
  modes <- modes[modes %in% s2$labels]
  cols <- col_accuracy[seq_along(modes)]
  pchs <- pch_accuracy[seq_along(modes)]
  ari <- tapply(s2$mean_ari, s2$labels, mean)[modes]
  acc_label <- c(oracle = "true types", est_easy = "SNR 1",
                 est_hard = "SNR 0.05", est_very_hard = "SNR 0.03")

  open_pdf(file.path(dir_fig, "fig_power.pdf"), height = 2.9)
  set_panels(c(1, 3), mar = c(2.9, 3.1, 1.5, 0.4), oma = c(2.6, 0, 0, 0))
  for (p in seq_along(patterns)) {
    plot(NA, xlim = range(s2$beta), ylim = c(0, 1), bty = "n",
         xlab = expression(beta), ylab = "Rejection rate")
    panel_title(letters[p], titles[p])
    abline(h = 0.05, lty = 3, col = col_reference)
    for (m in seq_along(modes)) {
      d <- s2[s2$pattern == patterns[p] & s2$labels == modes[m], ]
      d <- d[order(d$beta), ]
      lines(d$beta, d$power_centroid, col = cols[m], lwd = 1.1, lty = 2)
      lines(d$beta, d$power_global, col = cols[m], lwd = lwd_main)
      points(d$beta, d$power_global, col = cols[m], pch = pchs[m], cex = 0.8)
    }
    box(lwd = lwd_box)
  }
  legend_bottom(legend = c(paste0(acc_label[modes], " (ARI ", fmt(ari, 2), ")"),
                           "pair statistic", "centroid statistic"),
                col = c(cols, "grey20", "grey20"),
                pch = c(pchs, NA, NA),
                lty = c(rep(1, length(modes)), 1, 2),
                lwd = c(rep(lwd_main, length(modes)), lwd_main, 1.1),
                ncol = 3)
  dev.off()
  cat("written:", file.path(dir_fig, "fig_power.pdf"), "\n")
}

## =============================================================================
## Table: typing errors (sim_03)
## =============================================================================
f3 <- file.path("results", "sim_03_misclass_summary.csv")
if (file.exists(f3)) {
  s3 <- read.csv(f3)
  lines <- c("\\begin{tabular}{rrrr}",
             "\\toprule",
             "$\\varepsilon$ & $\\gamma^2$ & Global test & Type-1 test \\\\",
             "\\midrule")
  for (r in 1:nrow(s3)) {
    lines <- c(lines, paste0(fmt(s3$eps[r], 1), " & ", fmt(s3$gamma2[r]), " & ",
                             fmt(s3$power_global[r]), " & ",
                             fmt(s3$power_type1[r]), " \\\\"))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(dir_tab, "tab_misclass.tex"))
}

## =============================================================================
## Table: confounding (sim_04)
## =============================================================================
f4 <- file.path("results", "sim_04_confounding_summary.csv")
if (file.exists(f4)) {
  s4 <- read.csv(f4)
  if (all(is.na(s4$rejection_rate))) {
    warning("results/sim_04_confounding_summary.csv contains only NA: ",
            "the run of sim_04 failed. Rerun sim_04_confounding.R.")
  }
  nice <- c("raw + unrestricted" = "Raw curves, unrestricted permutations",
            "raw + stratified on z" = "Raw curves, permutations within strata of $z$",
            "spectrum + unrestricted" = "Fourier modulus, unrestricted permutations")
  lines <- c("\\begin{tabular}{lrr}",
             "\\toprule",
             "Analysis & Rejection rate & ARI \\\\",
             "\\midrule")
  for (r in 1:nrow(s4)) {
    ari_r <- round(s4$mean_ari[r], 2)
    if (!is.na(ari_r) && ari_r == 0) ari_r <- 0   # avoid printing "-0.00"
    ari_txt <- fmt(ari_r, 2)
    lines <- c(lines, paste0(nice[[s4$analysis[r]]], " & ",
                             fmt(s4$rejection_rate[r]), " & ", ari_txt, " \\\\"))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(dir_tab, "tab_confounding.tex"))
}

## =============================================================================
## Competitors (sim_05): full table (appendix) and power figure (main text)
## =============================================================================
f5 <- file.path("results", "sim_05_competitors_summary.csv")
tests5 <- c("fcst", "centroid", "segregation", "mscan", "fscan", "markK", "dcov")
test_names5 <- c("FCST", "Centroid", "Segregation", "Mult. scan",
                 "Funct. scan", "Mark $K$", "dCov")
pattern_name <- c(one_bump = "one hotspot", two_bumps = "two hotspots",
                  band = "band")
if (file.exists(f5)) {
  s5 <- read.csv(f5)
  ## ---- table: the null cells do not depend on the pattern, show them once ----
  lines <- c("\\begin{tabular}{llr|rrrr|rrr}",
             "\\toprule",
             " & & & \\multicolumn{4}{c|}{Tests on estimated types} & \\multicolumn{3}{c}{Tests on curves} \\\\",
             paste0("Pattern & SNR & $\\beta$ & ",
                    paste(test_names5, collapse = " & "), " \\\\"),
             "\\midrule")
  for (snr in sort(unique(s5$snr), decreasing = TRUE)) {
    d0 <- s5[s5$beta == 0 & s5$snr == snr, ][1, ]
    lines <- c(lines, paste0("null & ", snr, " & 0 & ",
                             paste(fmt(unlist(d0[tests5])), collapse = " & "),
                             " \\\\"))
  }
  lines <- c(lines, "\\midrule")
  for (p in c("one_bump", "two_bumps", "band")) {
    for (snr in sort(unique(s5$snr), decreasing = TRUE)) {
      for (b in sort(unique(s5$beta[s5$beta > 0]))) {
        d <- s5[s5$pattern == p & s5$snr == snr & s5$beta == b, ]
        lines <- c(lines, paste0(pattern_name[[p]], " & ", snr, " & ", b, " & ",
                                 paste(fmt(unlist(d[tests5])), collapse = " & "),
                                 " \\\\"))
      }
    }
    if (p != "band") lines <- c(lines, "\\midrule")
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(dir_tab, "tab_competitors_full.tex"))

  ## ---- figure: rows = SNR, columns = pattern ---------------------------------
  open_pdf(file.path(dir_fig, "fig_competitors.pdf"), height = 5.3)
  set_panels(c(2, 3), mar = c(2.9, 3.1, 1.5, 0.4), oma = c(3.0, 0, 0, 0))
  panel <- 0
  for (snr in sort(unique(s5$snr), decreasing = TRUE)) {
    for (p in c("one_bump", "two_bumps", "band")) {
      panel <- panel + 1
      # the null cell (pattern "none") is the beta = 0 point of every curve
      d <- s5[(s5$pattern == p | s5$pattern == "none") & s5$snr == snr, ]
      d <- d[order(d$beta), ]
      plot(NA, xlim = range(s5$beta), ylim = c(0, 1), bty = "n",
           xlab = expression(beta), ylab = "Rejection rate")
      panel_title(letters[panel], paste0(c(one_bump = "One hotspot",
                                           two_bumps = "Two hotspots",
                                           band = "Elongated band")[[p]], ", SNR ", snr))
      abline(h = 0.05, lty = 3, col = col_reference)
      for (t in rev(tests_all)) {                 # FCST drawn last, on top
        lines(d$beta, d[[t]], col = col_tests[[t]], lty = lty_tests[[t]],
              lwd = lwd_tests[[t]])
        points(d$beta, d[[t]], col = col_tests[[t]], pch = pch_tests[[t]],
               cex = 0.75)
      }
      box(lwd = lwd_box)
    }
  }
  legend_bottom(legend = tests_name[tests_all],
                col = col_tests[tests_all], lty = lty_tests[tests_all],
                pch = pch_tests[tests_all], lwd = lwd_tests[tests_all],
                pt.cex = 0.8, ncol = 4)
  dev.off()
  cat("written:", file.path(dir_fig, "fig_competitors.pdf"), "\n")
}

## =============================================================================
## Several typings (sim_06)
## =============================================================================
f6 <- file.path("results", "sim_06_typings_summary.csv")
if (file.exists(f6)) {
  s6 <- read.csv(f6)
  algs <- sub("^ari_", "", grep("^ari_", names(s6), value = TRUE))
  alg_name <- c(fpca_kmeans = "$k$-means", fpca_tkmeans = "trimmed $k$-means",
                fpca_mclust = "Gaussian mixture", ward = "Ward", pam = "PAM",
                fpca_ward = "Ward", fpca_pam = "PAM",
                fpcac = "FPCAC", clusteff = "clustEff")
  lines <- c(paste0("\\begin{tabular}{lr", strrep("r", length(algs)), "|rrr}"),
             "\\toprule",
             paste0("Pattern & $\\beta$ & ",
                    paste(alg_name[algs], collapse = " & "),
                    " & Min $p$ & Bonferroni & Joint \\\\"),
             "\\midrule",
             paste0("\\multicolumn{2}{l}{mean ARI} & ",
                    paste(fmt(colMeans(s6[, paste0("ari_", algs)]), 2),
                          collapse = " & "), " & & & \\\\"),
             "\\midrule")
  for (r in 1:nrow(s6)) {
    vals <- c(unlist(s6[r, paste0("rej_", algs)]), s6$rej_naive_min[r],
              s6$rej_bonferroni[r], s6$rej_joint[r])
    lines <- c(lines, paste0(pattern_name[[s6$pattern[r]]], " & ", s6$beta[r],
                             " & ", paste(fmt(vals), collapse = " & "), " \\\\"))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(dir_tab, "tab_typings.tex"))
}

## =============================================================================
## Simulations sim_04b, sim_05b, sim_08
## =============================================================================

## ---- Confounding scenarios (sim_04b) --------------------------------------------
f4b <- file.path("results", "sim_04b_confounding_summary.csv")
if (file.exists(f4b)) {
  s4b <- read.csv(f4b)
  get4b <- function(scen, column) s4b[s4b$scenario == scen, column]
  row4b <- function(label, typing, expected_u, expected_s, scen, prefix, ari_col) {
    vals <- c(get4b(scen, ari_col), get4b(scen, paste0(prefix, "_unrestricted")),
              get4b(scen, paste0(prefix, "_strata8")),
              get4b(scen, paste0(prefix, "_strata16")))
    paste0(label, " & ", typing, " & ", expected_u, " & ", expected_s, " & ",
           fmt(vals[1], 2), " & ", paste(fmt(vals[2:4]), collapse = " & "), " \\\\")
  }
  lines <- c("\\begin{tabular}{@{}llllrrrr@{}}",
             "\\toprule",
             " & & \\multicolumn{2}{c}{Expected} & & \\multicolumn{3}{c}{Rejection rate} \\\\",
             "\\cmidrule(lr){3-4}\\cmidrule(lr){6-8}",
             "Scenario & Typing & Unrestr. & Stratified & ARI & Unrestr. & 8 strata & 16 strata \\\\",
             "\\midrule",
             row4b("(1) $z$ acts on amplitude", "raw curves", "rejects", "level",
                   "amplitude", "raw", "ari_raw"),
             row4b("", "Fourier modulus", "level", "level", "amplitude", "spec", "ari_spec"),
             row4b("(2) type depends on $z$", "Fourier modulus", "rejects", "level",
                   "type_z", "spec", "ari_spec"),
             row4b("(3) $z$ and hotspot", "Fourier modulus", "rejects", "rejects",
                   "type_z_hot", "spec", "ari_spec"),
             row4b("(4) $z$ and abscissa", "Fourier modulus", "rejects", "rejects",
                   "type_z_x", "spec", "ari_spec"),
             row4b("(5) $z$ acts on the SNR", "Fourier modulus", "may reject", "level",
                   "snr_z", "spec", "ari_spec"),
             "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(dir_tab, "tab_confounding_scenarios.tex"))
}

## ---- Segregation test with cross-validated bandwidth (sim_05b) ------------------
f5b <- file.path("results", "sim_05b_segregation_cv_summary.csv")
if (file.exists(f5b)) {
  s5b <- read.csv(f5b)
  geom_name <- c(none = "null", one_bump = "one hotspot",
                 two_bumps = "two hotspots", band = "band")
  lines <- c("\\begin{tabular}{@{}lrrrrrrrr@{}}",
             "\\toprule",
             " & & & & & \\multicolumn{3}{c}{Segregation test} & \\\\",
             "\\cmidrule(lr){6-8}",
             "Geometry & $\\beta$ & SNR & ARI & FCST & fixed & CV & CV, exact & bandwidth (fixed / CV) \\\\",
             "\\midrule")
  for (r in 1:nrow(s5b)) {
    lines <- c(lines, paste0(geom_name[[s5b$pattern[r]]], " & ", s5b$beta[r], " & ",
                             s5b$snr[r], " & ", fmt(s5b$mean_ari[r], 2), " & ",
                             paste(fmt(unlist(s5b[r, c("fcst", "seg_fixed", "seg_cv",
                                                        "seg_cvp")])), collapse = " & "),
                             " & ", fmt(s5b$median_bw_fixed[r], 2), " / ",
                             fmt(s5b$median_bw_cv[r], 2), " \\\\"))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(dir_tab, "tab_segregation_cv.tex"))
}

## ---- Robustness to ground process and hotspot profile (sim_08) ------------------
f8 <- file.path("results", "sim_08_robustness_summary.csv")
if (file.exists(f8)) {
  s8 <- read.csv(f8)
  proc8 <- c(hpp = "Homogeneous Poisson", ipp = "Inhomogeneous Poisson",
             thomas = "Thomas", thomas_overlap = "Thomas, overlapping")
  prof8 <- c(gauss = "Gaussian", t = "Student $t$")
  get8 <- function(proc, geom, prof, beta, lab, w, column) {
    sel <- s8$process == proc & s8$geometry == geom & s8$profile == prof &
      s8$beta == beta & s8$labels == lab & s8$weights == w
    return(s8[sel, column])
  }
  head8 <- c("\\toprule",
             " & & \\multicolumn{4}{c}{One hotspot} & \\multicolumn{4}{c}{Two hotspots} \\\\",
             "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
             " & & \\multicolumn{2}{c}{unit} & \\multicolumn{2}{c}{inverse} & \\multicolumn{2}{c}{unit} & \\multicolumn{2}{c}{inverse} \\\\",
             "Process & Profile & global & type 1 & global & type 1 & global & type 1 & global & type 1 \\\\",
             "\\midrule")
  rows8 <- function(beta, lab, first_col = TRUE) {
    out <- character(0)
    for (proc in names(proc8)) {
      for (prof in names(prof8)) {
        vals <- c()
        for (geom in c("one_bump", "two_bumps")) {
          for (w in c("none", "inverse_intensity")) {
            vals <- c(vals, get8(proc, geom, prof, beta, lab, w, "power_global"),
                      get8(proc, geom, prof, beta, lab, w, "holm_type1"))
          }
        }
        lab_proc <- ifelse(prof == "gauss", proc8[[proc]], "")
        out <- c(out, paste0(lab_proc, " & ", prof8[[prof]], " & ",
                             paste(fmt(vals, 2), collapse = " & "), " \\\\"))
      }
    }
    return(out)
  }
  ## main text: true types, beta = 6
  lines <- c("\\begin{tabular}{@{}llrrrrrrrr@{}}", head8, rows8(6, "oracle"),
             "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(dir_tab, "tab_robustness.tex"))

  ## appendix: all alternatives
  lines <- c("\\begin{tabular}{@{}llrrrrrrrr@{}}", head8)
  for (lab in c("oracle", "est_hard")) {
    for (beta in c(3, 6)) {
      lines <- c(lines, paste0("\\multicolumn{10}{@{}l}{\\emph{",
                               ifelse(lab == "oracle", "true types",
                                      "estimated types, SNR 0.05"),
                               ", $\\beta=", beta, "$}} \\\\"),
                 rows8(beta, lab))
    }
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(dir_tab, "tab_robustness_full.tex"))

  ## appendix: size under beta = 0
  lines <- c("\\begin{tabular}{@{}lrrrrrrrr@{}}",
             "\\toprule",
             " & \\multicolumn{4}{c}{One hotspot (geometry of the parents)} & \\multicolumn{4}{c}{Two hotspots} \\\\",
             "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
             " & \\multicolumn{2}{c}{true types} & \\multicolumn{2}{c}{estimated} & \\multicolumn{2}{c}{true types} & \\multicolumn{2}{c}{estimated} \\\\",
             "Process & unit & inverse & unit & inverse & unit & inverse & unit & inverse \\\\",
             "\\midrule")
  for (proc in names(proc8)) {
    vals <- c()
    for (geom in c("one_bump", "two_bumps")) {
      for (lab in c("oracle", "est_hard")) {
        for (w in c("none", "inverse_intensity")) {
          vals <- c(vals, get8(proc, geom, "none", 0, lab, w, "power_global"))
        }
      }
    }
    lines <- c(lines, paste0(proc8[[proc]], " & ", paste(fmt(vals), collapse = " & "), " \\\\"))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(dir_tab, "tab_robustness_size.tex"))
}

## ---- Full power results of sim_02 (Supplementary Material) ----------------------
## Longtable body (the longtable environment and caption are in the supplement).
f2s <- file.path("results", "sim_02_power_summary.csv")
if (file.exists(f2s)) {
  s2f <- read.csv(f2s)
  geom2 <- c(one_bump = "one hotspot", two_bumps = "two hotspots", band = "band")
  lab2 <- c(oracle = "true", est_easy = "SNR 1", est_hard = "SNR 0.05",
            est_very_hard = "SNR 0.03")
  lines <- character(0)
  for (r in 1:nrow(s2f)) {
    vals <- unlist(s2f[r, c("power_global", "power_maxT", "power_gauss",
                            "power_centroid", "holm_type1", "holm_other")])
    lines <- c(lines, paste0(geom2[[s2f$pattern[r]]], " & ", s2f$beta[r], " & ",
                             lab2[[s2f$labels[r]]], " & ", fmt(s2f$mean_ari[r], 2),
                             " & ", paste(fmt(vals), collapse = " & "), " \\\\"))
  }
  write_lines(lines, file.path(dir_tab, "tab_power_full_body.tex"))
}
