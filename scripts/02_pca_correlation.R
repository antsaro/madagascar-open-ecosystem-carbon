# 02_pca_correlation.R
# PCA of environmental covariates and carbon pools, PERMANOVA, VIF and Spearman
# correlations. Writes figures and tables to outputs/pca.
# Run from the project root (the .here file marks it).

# Setup ----

library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)
library(vegan)
library(ragg)
library(systemfonts)

data_path <- here("data", "Datamapf_coordo.csv")
fig_dir <- here("outputs", "pca", "figures")
tab_dir <- here("outputs", "pca", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# kable in the rendered report, plain print when the script is sourced
tab <- function(x, caption = NULL, digits = 2) {
  if (isTRUE(getOption("knitr.in.progress")))
    return(knitr::kable(x, caption = caption, digits = digits))
  if (!is.null(caption)) cat("\n", caption, "\n", sep = "")
  if (inherits(x, "tbl")) print(x, n = 20) else print(x)
}

fams <- unique(systemfonts::system_fonts()$family)
fam <- intersect(c("Arial", "Liberation Sans", "Helvetica"), fams)[1]
if (is.na(fam)) fam <- "sans"

pal_site <- c(AF = "#E64B35", AT = "#4DBBD5", IV = "#00A087")
pal_lc <- c(Forest = "#00A087", Shrubland = "#F2A93B", Grassland = "#3C5488", Reforestation = "#8F5FA8")
shp_lc <- c(Forest = 16, Shrubland = 17, Grassland = 15, Reforestation = 18)
shp_site <- c(AF = 16, AT = 17, IV = 15)
col_sup <- "#B2182B"

theme_set(
  theme_classic(base_size = 17, base_family = fam) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.6),
          axis.text = element_text(colour = "black"),
          legend.position = "bottom", legend.box = "vertical",
          legend.box.spacing = unit(0.2, "cm"), legend.spacing.y = unit(0, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          legend.title = element_text(face = "bold"),
          legend.key.size = unit(0.6, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 17),
          plot.tag = element_text(face = "bold", size = 24),
          plot.margin = margin(8, 12, 8, 8))
)

save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         units = "in", dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h,
         units = "in", device = cairo_pdf)
  invisible(p)
}

# Data ----

dat <- read_delim(data_path, delim = ";", show_col_types = FALSE, trim_ws = TRUE) |>
  rename(SOC30 = ST30, SOC100 = ST100, SOC30_100 = ST30.100, AGC = AGCGIEC,
         BGC = BGCGIEC, DW = DWGIEC, Litter = LITERGIEC, TotC = TOTC,
         Land_cover = LU, Elevation = MNT, Fire = FF,
         Cos_aspect = Cosaspect, Sin_aspect = Sinaspect) |>
  mutate(Site = factor(Site, levels = c("AF", "AT", "IV")),
         Land_cover = factor(Land_cover, levels = c("Forest", "Shrubland", "Grassland", "Reforestation")))

lab_env <- c(Elevation = "Elevation", Slope = "Slope", Edge = "Edge", Roughness = "Roughness",
             TPI = "TPI", TRI = "TRI", MAP = "MAP", MAT = "MAT",
             Cos_aspect = "cos(Aspect)", Sin_aspect = "sin(Aspect)", Fire = "Fire frequency")
lab_pool <- c(SOC30 = "SOC 0\u201330 cm", SOC30_100 = "SOC 30\u2013100 cm", SOC100 = "SOC 0\u2013100 cm",
              AGC = "AGC", BGC = "BGC", DW = "Deadwood", Litter = "Litter", TotC = "Total C")

env_active <- setdiff(names(lab_env), "Fire")        # PCA input
env_cor <- names(lab_env)                            # correlations (fire frequency is ordinal)
pool_active <- c("SOC30", "SOC30_100", "AGC", "BGC", "DW", "Litter")
pool_cor <- names(lab_pool)

env_raw <- dat |> select(all_of(env_cor))
pool_raw <- dat |> select(all_of(pool_cor))

# edge distance and carbon pools are right-skewed
env_mat <- env_raw |> mutate(Edge = log1p(Edge))
pool_mat <- pool_raw |> mutate(across(everything(), log1p))

cc <- complete.cases(pool_mat[, pool_active])

tibble(`Plots` = nrow(dat), `Plots in carbon PCA` = sum(cc)) |> tab()
dat |> count(Site, Land_cover) |>
  pivot_wider(names_from = Land_cover, values_from = n, values_fill = 0) |>
  tab(caption = "Plots per site and land cover")

# Collinearity ----

vif_rank <- function(X) {
  R <- cor(X, method = "spearman", use = "pairwise.complete.obs")
  setNames(diag(solve(R)), colnames(X))
}

vif_tab <- tibble(Variable = lab_env[env_active], VIF = vif_rank(env_raw[env_active])) |>
  arrange(desc(VIF))
tab(vif_tab, caption = "Rank-based variance inflation factors (VIF > 10: strong collinearity)", digits = 1)
write_csv(vif_tab, file.path(tab_dir, "Table_VIF.csv"))

# PCA ----

pca_summary <- function(pca, vlab) {
  ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  eig <- tibble(PC = paste0("PC", seq_along(ve)), Eigenvalue = pca$sdev^2,
                `Variance (%)` = ve, `Cumulative (%)` = cumsum(ve))
  # variable-axis correlations
  load <- as_tibble(sweep(pca$rotation[, 1:4], 2, pca$sdev[1:4], "*"), rownames = "Variable") |>
    mutate(Variable = vlab[Variable])
  list(eig = eig, load = load)
}

pca_env <- prcomp(env_mat[, env_active], center = TRUE, scale. = TRUE)
pca_pool <- prcomp(pool_mat[cc, pool_active], center = TRUE, scale. = TRUE)
s_env <- pca_summary(pca_env, lab_env)
s_pool <- pca_summary(pca_pool, lab_pool)

tab(s_env$eig, "Eigenvalues, environmental PCA")
tab(s_env$load, "Variable-axis correlations, environmental PCA")
tab(s_pool$eig, "Eigenvalues, carbon-pool PCA")
tab(s_pool$load, "Variable-axis correlations, carbon-pool PCA")

write_csv(s_env$eig, file.path(tab_dir, "Table_PCA_env_eigenvalues.csv"))
write_csv(s_env$load, file.path(tab_dir, "Table_PCA_env_loadings.csv"))
write_csv(s_pool$eig, file.path(tab_dir, "Table_PCA_pool_eigenvalues.csv"))
write_csv(s_pool$load, file.path(tab_dir, "Table_PCA_pool_loadings.csv"))

# PERMANOVA ----

set.seed(2026)
meta <- dat |> select(Site, Land_cover)

perm_env <- adonis2(scale(as.matrix(env_mat[, env_active])) ~ Site + Land_cover, data = meta,
                    method = "euclidean", by = "margin", permutations = 9999)
perm_pool <- adonis2(scale(as.matrix(pool_mat[cc, pool_active])) ~ Site + Land_cover,
                     data = meta[cc, ], method = "euclidean", by = "margin", permutations = 9999)

perm_tab <- bind_rows(
  as_tibble(perm_env, rownames = "Term") |> mutate(Space = "Environmental covariates"),
  as_tibble(perm_pool, rownames = "Term") |> mutate(Space = "Carbon pools")
) |>
  filter(Term != "Total") |>
  select(Space, Term, Df, SumOfSqs, R2, `F`, `Pr(>F)`)

tab(perm_tab, caption = "Marginal PERMANOVA, 9999 permutations", digits = 3)
write_csv(perm_tab, file.path(tab_dir, "Table_PERMANOVA.csv"))

# Biplot function ----

pca_biplot <- function(pca, meta, col_var, shape_var, pal, shp, col_name, shp_name,
                       vlab, sup = NULL, slab = NULL, axes = c(1, 2), txt = 4.8) {
  ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  sc <- tibble(x = pca$x[, axes[1]], y = pca$x[, axes[2]]) |> bind_cols(meta)

  ld <- sweep(pca$rotation[, axes, drop = FALSE], 2, pca$sdev[axes], "*")
  ld <- tibble(var = rownames(ld), x = ld[, 1], y = ld[, 2], lab = vlab[rownames(ld)])

  # stretch the arrows to the size of the score cloud
  k <- 0.9 * min(max(abs(sc$x)) / max(abs(ld$x)), max(abs(sc$y)) / max(abs(ld$y)))

  p <- ggplot(sc, aes(x, y)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
    stat_ellipse(aes(colour = .data[[col_var]], fill = .data[[col_var]]), geom = "polygon",
                 type = "norm", level = 0.95, alpha = 0.10, linewidth = 0.8, show.legend = FALSE) +
    geom_point(aes(colour = .data[[col_var]], shape = .data[[shape_var]]), size = 3, alpha = 0.8) +
    geom_segment(data = ld, aes(x = 0, y = 0, xend = x * k, yend = y * k),
                 arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
                 colour = "grey10", linewidth = 0.75, inherit.aes = FALSE) +
    geom_text_repel(data = ld, aes(x * k, y * k, label = lab), family = fam, fontface = "bold",
                    size = txt, colour = "grey10", box.padding = 0.35, point.padding = 0.2,
                    segment.size = 0.3, min.segment.length = 0.3, max.overlaps = Inf,
                    inherit.aes = FALSE)

  if (!is.null(sup)) {
    sp <- cor(sup, as.matrix(sc[, c("x", "y")]), use = "pairwise.complete.obs")
    sp <- tibble(var = rownames(sp), x = sp[, 1] * k, y = sp[, 2] * k, lab = slab[rownames(sp)])
    p <- p +
      geom_segment(data = sp, aes(x = 0, y = 0, xend = x, yend = y),
                   arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
                   colour = col_sup, linetype = "dashed", linewidth = 0.75, inherit.aes = FALSE) +
      geom_text_repel(data = sp, aes(x, y, label = lab), family = fam, fontface = "italic",
                      size = txt, colour = col_sup, box.padding = 0.35, segment.size = 0.3,
                      min.segment.length = 0.3, max.overlaps = Inf, inherit.aes = FALSE)
  }

  p +
    scale_colour_manual(values = pal, name = col_name) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_shape_manual(values = shp, name = shp_name) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 4.5, alpha = 1)),
           shape = guide_legend(order = 2, override.aes = list(size = 4.5))) +
    labs(x = sprintf("PC%d (%.1f%%)", axes[1], ve[axes[1]]),
         y = sprintf("PC%d (%.1f%%)", axes[2], ve[axes[2]])) +
    coord_fixed()
}

# Environmental PCA ----

sup_env <- bind_cols(pool_mat, Fire = env_raw$Fire)
fig1 <- pca_biplot(pca_env, meta = select(dat, Site, Land_cover),
                   col_var = "Site", shape_var = "Land_cover",
                   pal = pal_site, shp = shp_lc, col_name = "Site", shp_name = "Land cover",
                   vlab = lab_env, sup = sup_env, slab = c(lab_pool, Fire = "Fire frequency"))
save_fig(fig1, "Fig_PCA_environment", 9.5, 9)
fig1

# Carbon-pool PCA ----

sup_pool <- env_mat[cc, c("MAT", "MAP", "Fire", "Slope", "Edge")]
fig2 <- pca_biplot(pca_pool, meta = select(dat[cc, ], Site, Land_cover),
                   col_var = "Land_cover", shape_var = "Site",
                   pal = pal_lc, shp = shp_site, col_name = "Land cover", shp_name = "Site",
                   vlab = lab_pool, sup = sup_pool, slab = lab_env)
save_fig(fig2, "Fig_PCA_carbon_pools", 9.5, 9)
fig2

# Correlation matrices ----

cor_long <- function(X, Y) {
  expand_grid(x = names(X), y = names(Y)) |>
    mutate(res = map2(x, y, \(a, b) {
      ok <- complete.cases(X[[a]], Y[[b]])
      ct <- suppressWarnings(cor.test(X[[a]][ok], Y[[b]][ok], method = "spearman", exact = FALSE))
      tibble(r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
    })) |>
    unnest(res)
}

stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

# covariates ordered by clustering on |rho|
ord_env <- env_cor[hclust(as.dist(1 - abs(cor(env_raw, method = "spearman"))), method = "average")$order]

heat_plot <- function(df, xlev, ylev, xlab, ylab, txt = 3.9, rot = 45) {
  ggplot(df, aes(factor(x, xlev), factor(y, rev(ylev)), fill = r)) +
    geom_tile(colour = "white", linewidth = 0.9) +
    geom_text(aes(label = paste0(sub("^-0\\.00$", "0.00", sprintf("%.2f", r)), stars(p_adj)),
                  colour = abs(r) > 0.6), family = fam, size = txt) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0,
                         limits = c(-1, 1), name = expression("Spearman " * rho),
                         guide = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                                 barwidth = unit(9, "cm"), barheight = unit(0.5, "cm"),
                                                 ticks.colour = "black")) +
    scale_x_discrete(labels = xlab, drop = FALSE) +
    scale_y_discrete(labels = ylab, drop = FALSE) +
    labs(x = NULL, y = NULL) +
    coord_fixed() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = rot, hjust = 1, vjust = 1),
          panel.grid = element_blank())
}

# (a) covariates, lower triangle only; (b) pools x covariates
cor_env <- cor_long(env_raw, env_raw) |>
  mutate(ix = match(x, ord_env), iy = match(y, ord_env)) |>
  filter(ix < iy) |>
  mutate(p_adj = p.adjust(p, "BH"))

cor_pe <- cor_long(env_raw, pool_raw) |>
  mutate(p_adj = p.adjust(p, "BH"))

write_csv(select(cor_env, Var1 = y, Var2 = x, rho = r, n, p, p_adj),
          file.path(tab_dir, "Table_corr_env_env.csv"))
write_csv(select(cor_pe, Covariate = x, Pool = y, rho = r, n, p, p_adj),
          file.path(tab_dir, "Table_corr_pool_env.csv"))

p3a <- heat_plot(cor_env, ord_env, ord_env, lab_env, lab_env)
p3b <- heat_plot(cor_pe, ord_env, pool_cor, lab_env, lab_pool)

fig3 <- (p3a / p3b) +
  plot_layout(heights = c(10, 7), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")
save_fig(fig3, "Fig_correlation_matrices", 12, 16)
fig3

# Scree plots and contributions ----

scree_plot <- function(pca, n = 8) {
  ev <- tibble(PC = factor(paste0("PC", seq_along(pca$sdev)), levels = paste0("PC", seq_along(pca$sdev))),
               ve = 100 * pca$sdev^2 / sum(pca$sdev^2)) |> slice_head(n = n)
  ggplot(ev, aes(PC, ve)) +
    geom_col(fill = "#3C5488", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f", ve)), vjust = -0.5, family = fam, size = 5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Variance explained (%)")
}

contrib_plot <- function(pca, vlab, pc = 1) {
  ct <- tibble(var = rownames(pca$rotation), contrib = 100 * pca$rotation[, pc]^2,
               dir = ifelse(pca$rotation[, pc] >= 0, "Positive", "Negative")) |>
    mutate(lab = vlab[var])
  ggplot(ct, aes(reorder(lab, contrib), contrib, fill = dir)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 100 / nrow(ct), linetype = "dashed", colour = "grey30") +
    scale_fill_manual(values = c(Positive = "#E64B35", Negative = "#3C5488"), name = "Loading",
                      drop = FALSE, limits = c("Negative", "Positive")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(x = NULL, y = sprintf("Contribution to PC%d (%%)", pc))
}

figS1 <- (scree_plot(pca_env) | contrib_plot(pca_env, lab_env, 1) | contrib_plot(pca_env, lab_env, 2)) /
  (scree_plot(pca_pool) | contrib_plot(pca_pool, lab_pool, 1) | contrib_plot(pca_pool, lab_pool, 2)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")
save_fig(figS1, "FigS_scree_contributions", 17, 10)
figS1

# Correlations by site ----

cor_site <- dat |>
  group_by(Site) |>
  group_modify(\(d, key) cor_long(d[env_cor], d[pool_cor])) |>
  mutate(p_adj = p.adjust(p, "BH")) |>
  ungroup()

write_csv(select(cor_site, Site, Covariate = x, Pool = y, rho = r, n, p, p_adj),
          file.path(tab_dir, "Table_corr_pool_env_by_site.csv"))

figS2 <- heat_plot(cor_site, ord_env, pool_cor, lab_env, lab_pool, txt = 3.2, rot = 60) +
  facet_wrap(~Site, nrow = 1)
save_fig(figS2, "FigS_corr_pool_env_by_site", 20, 7.5)
figS2

# Session information ----

sessionInfo()
