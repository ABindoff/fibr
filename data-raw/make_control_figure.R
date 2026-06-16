## data-raw/make_control_figure.R
## Negative-control comparison figure for the manuscript (Task 6).
##
## Reads : data-raw/control_pairs.rds
## Writes: data-raw/control_pairs_manuscript.png
## Run from package root: Rscript data-raw/make_control_figure.R

library(ggplot2)

out_dir <- "data-raw"
df      <- readRDS(file.path(out_dir, "control_pairs.rds"))

## ── Cell labels ordered high -> low pi_j ─────────────────────────────────────
cell_order <- c("nj3_s0.5", "sparse_benchmark", "nj3_s3.0", "nj100_s0.5")
cell_labs  <- c(
  "nj3_s0.5"         = "n[j]==3~','~sigma==0.5~~(high~pi[j])",
  "sparse_benchmark" = "sparse~benchmark~(n[j]==3~','~sigma==3)",
  "nj3_s3.0"         = "n[j]==3~','~sigma==3.0~~(moderate~pi[j])",
  "nj100_s0.5"       = "n[j]==100~','~sigma==0.5~~(low~pi[j])"
)
df$cell_label <- factor(df$cell,
  levels = cell_order,
  labels = cell_labs[cell_order]
)
df$group <- factor(df$group)

## ── Control range per cell x gap for the grey reference band ─────────────────
ctrl_summ <- do.call(rbind, Filter(Negate(is.null), lapply(
  split(df, list(df$cell_label, df$gap)),
  function(sub) {
    if (nrow(sub) == 0L) return(NULL)
    data.frame(
      cell_label = as.character(sub$cell_label[1L]),
      gap        = sub$gap[1L],
      ctrl_lo    = min(sub$h_control),
      ctrl_hi    = max(sub$h_control),
      ctrl_mean  = mean(sub$h_control)
    )
  }
)))
ctrl_summ$cell_label <- factor(ctrl_summ$cell_label,
  levels = levels(df$cell_label))
ctrl_summ <- ctrl_summ[order(ctrl_summ$cell_label, ctrl_summ$gap), ]

## ── Figure ────────────────────────────────────────────────────────────────────
fig <- ggplot(df, aes(x = gap, colour = group, group = group)) +
  ## Grey band = full range of control-pair rho_j across 8 groups
  geom_ribbon(
    data        = ctrl_summ,
    aes(x = gap, ymin = ctrl_lo, ymax = ctrl_hi),
    fill        = "grey75",
    alpha       = 0.50,
    inherit.aes = FALSE
  ) +
  ## Dashed line = control mean
  geom_line(
    data        = ctrl_summ,
    aes(x = gap, y = ctrl_mean),
    colour      = "grey40",
    linetype    = "dashed",
    linewidth   = 0.6,
    inherit.aes = FALSE
  ) +
  ## Dotted reference at zero
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
  ## True-loop rho_j per group: lines + points
  geom_line(aes(y = h_true),  linewidth = 0.55, alpha = 0.80) +
  geom_point(aes(y = h_true), size = 1.9,  shape = 16, alpha = 0.90) +
  ## Facets (2 x 2, ordered high -> low pi_j)
  facet_wrap(~ cell_label, ncol = 2, labeller = label_parsed) +
  ## Log-scale x so gaps 3/10/25/50 are evenly spaced
  scale_x_continuous(
    trans  = "log10",
    breaks = c(3, 10, 25, 50),
    labels = c("3", "10", "25", "50")
  ) +
  scale_colour_viridis_d(option = "turbo", name = "Group j") +
  labs(
    x = "Gap (log scale)",
    y = expression(hat(rho)[j] ~ "(loop-conditional dependence)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "right",
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "control_pairs_manuscript.png"),
  plot   = fig,
  width  = 9,
  height = 7,
  dpi    = 150
)
cat("Saved: data-raw/control_pairs_manuscript.png\n")
