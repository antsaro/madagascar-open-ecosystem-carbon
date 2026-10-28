# 01_statistical_tests.R
# Assumption checks and one-way tests (ANOVA / Welch / Kruskal-Wallis) of carbon pools
# between land covers and sites. Writes tables to outputs/stat_tests.
# Run from the project root (the .here file marks it).

# Setup ----

library(tidyverse)
library(here)
library(car)
library(rstatix)
library(openxlsx)
library(flextable)
library(officer)
library(systemfonts)
library(ragg)

data_path <- here("data", "Datamapf_coordo.csv")
out_dir <- here("outputs", "stat_tests")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

alpha_assump <- 0.05
alpha_test <- 0.05
set.seed(2026)

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

theme_set(
  theme_classic(base_size = 14, base_family = fam) +
    theme(axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
)

# Data ----

sites <- c("AF", "AT", "IV")
lc_levels <- c("Forest", "Shrubland", "Grassland", "Reforestation")

dat <- read_delim(data_path, delim = ";", show_col_types = FALSE, trim_ws = TRUE) |>
  rename(SOC30 = ST30, SOC100 = ST100, SOC30_100 = ST30.100, AGC = AGCGIEC,
         BGC = BGCGIEC, DW = DWGIEC, Litter = LITERGIEC, TotC = TOTC, Land_cover = LU) |>
  mutate(Site = factor(Site, levels = sites),
         Land_cover = factor(Land_cover, levels = lc_levels))

pools <- c(SOC30 = "SOC 0\u201330 cm", SOC30_100 = "SOC 30\u2013100 cm", SOC100 = "SOC 0\u2013100 cm",
           AGC = "AGC", BGC = "BGC", DW = "Deadwood", Litter = "Litter", TotC = "Total C")

# land cover (all sites), site (all land covers), land cover within each site
design <- bind_rows(
  tibble(Factor = "Land cover", Subset = "All sites", group = "Land_cover"),
  tibble(Factor = "Site", Subset = "All land covers", group = "Site"),
  tibble(Factor = "Land cover", Subset = paste("Site", sites), group = "Land_cover")
)

dat |> count(Site, Land_cover) |>
  pivot_wider(names_from = Land_cover, values_from = n, values_fill = 0) |>
  tab(caption = "Plots per site and land cover")

# Assumptions and test choice ----

fmt_p <- function(p) ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))

check_assump <- function(d, y, g, transf = c("none", "log1p")) {
  transf <- match.arg(transf)
  yy <- if (transf == "log1p") log1p(d[[y]]) else d[[y]]
  dd <- tibble(y = yy, g = d[[g]])
  res <- residuals(aov(y ~ g, data = dd))
  sw <- shapiro.test(res)
  lv <- car::leveneTest(y ~ g, data = dd, center = median)
  by_group <- dd |> group_by(g) |> filter(n() >= 3) |>
    summarise(p = shapiro.test(y)$p.value, .groups = "drop")
  tibble(Transformation = ifelse(transf == "log1p", "log(1 + x)", "none"),
         W = unname(sw$statistic), p_shapiro = sw$p.value,
         n_groups_nonnormal = sum(by_group$p < alpha_assump), n_groups = nrow(by_group),
         F_levene = lv$`F value`[1], df1_levene = lv$Df[1], df2_levene = lv$Df[2],
         p_levene = lv$`Pr(>F)`[1],
         normal = sw$p.value > alpha_assump, homosc = lv$`Pr(>F)`[1] > alpha_assump)
}

choose_test <- function(raw, lg) {
  if (raw$normal && raw$homosc) return(list(test = "ANOVA", transf = "none"))
  if (lg$normal && lg$homosc) return(list(test = "ANOVA", transf = "log1p"))
  if (raw$normal && !raw$homosc) return(list(test = "Welch ANOVA", transf = "none"))
  if (lg$normal && !lg$homosc) return(list(test = "Welch ANOVA", transf = "log1p"))
  list(test = "Kruskal-Wallis", transf = "none")
}

run_test <- function(d, y, g, test, transf) {
  yy <- if (transf == "log1p") log1p(d[[y]]) else d[[y]]
  dd <- tibble(y = yy, g = d[[g]])

  if (test == "ANOVA") {
    fit <- aov(y ~ g, data = dd)
    s <- summary(fit)[[1]]
    om <- tibble(Statistic = "F", Value = s$`F value`[1], df1 = s$Df[1], df2 = s$Df[2],
                 p = s$`Pr(>F)`[1], Effect_size = s$`Sum Sq`[1] / sum(s$`Sum Sq`),
                 Effect_metric = "eta-squared")
    ph <- NULL
    if (om$p < alpha_test) {
      tk <- TukeyHSD(fit)$g
      ph <- tibble(Comparison = rownames(tk), Estimate = tk[, "diff"], Lower95 = tk[, "lwr"],
                   Upper95 = tk[, "upr"], p_adj = tk[, "p adj"], Method = "Tukey HSD")
    }
  } else if (test == "Welch ANOVA") {
    w <- oneway.test(y ~ g, data = dd, var.equal = FALSE)
    ss <- summary(aov(y ~ g, data = dd))[[1]]$`Sum Sq`
    om <- tibble(Statistic = "F (Welch)", Value = unname(w$statistic), df1 = unname(w$parameter[1]),
                 df2 = unname(w$parameter[2]), p = w$p.value,
                 Effect_size = ss[1] / sum(ss), Effect_metric = "eta-squared (OLS)")
    ph <- NULL
    if (om$p < alpha_test) {
      gh <- rstatix::games_howell_test(dd, y ~ g)
      ph <- tibble(Comparison = paste(gh$group2, gh$group1, sep = "-"), Estimate = gh$estimate,
                   Lower95 = gh$conf.low, Upper95 = gh$conf.high, p_adj = gh$p.adj,
                   Method = "Games-Howell")
    }
  } else {
    kw <- kruskal.test(y ~ g, data = dd)
    om <- tibble(Statistic = "H", Value = unname(kw$statistic), df1 = unname(kw$parameter),
                 df2 = NA_real_, p = kw$p.value,
                 Effect_size = rstatix::kruskal_effsize(dd, y ~ g)$effsize,
                 Effect_metric = "epsilon-squared (H-based)")
    ph <- NULL
    if (om$p < alpha_test) {
      dn <- rstatix::dunn_test(dd, y ~ g, p.adjust.method = "holm")
      ph <- tibble(Comparison = paste(dn$group2, dn$group1, sep = "-"), Estimate = dn$statistic,
                   Lower95 = NA_real_, Upper95 = NA_real_, p_adj = dn$p.adj,
                   Method = "Dunn (Holm), estimate = z")
    }
  }
  list(omnibus = bind_cols(tibble(n = nrow(dd), k = nlevels(dd$g)), om), posthoc = ph)
}

res_assump <- list(); res_omni <- list(); res_ph <- list()

for (y in names(pools)) {
  for (i in seq_len(nrow(design))) {
    g <- design$group[i]
    lab <- design$Subset[i]
    d <- dat
    if (startsWith(lab, "Site ")) d <- filter(d, Site == sub("Site ", "", lab))
    d <- d |> filter(!is.na(.data[[y]]), !is.na(.data[[g]])) |> droplevels()
    if (nlevels(d[[g]]) < 2 || min(table(d[[g]])) < 3) next

    a_raw <- check_assump(d, y, g, "none")
    a_log <- check_assump(d, y, g, "log1p")
    ch <- choose_test(a_raw, a_log)
    rt <- run_test(d, y, g, ch$test, ch$transf)
    id <- tibble(Variable = pools[[y]], Factor = design$Factor[i], Subset = lab)

    res_assump[[length(res_assump) + 1]] <- bind_cols(id, n = nrow(d), bind_rows(a_raw, a_log))
    res_omni[[length(res_omni) + 1]] <- bind_cols(
      id, Test = ch$test, Transformation = ifelse(ch$transf == "log1p", "log(1 + x)", "none"),
      rt$omnibus)
    if (!is.null(rt$posthoc)) res_ph[[length(res_ph) + 1]] <- bind_cols(id, Test = ch$test, rt$posthoc)
  }
}

assump <- bind_rows(res_assump)
omni <- bind_rows(res_omni) |>
  group_by(Factor, Subset) |> mutate(p_BH = p.adjust(p, "BH")) |> ungroup()
posthoc <- bind_rows(res_ph)

test_summary <- omni |>
  select(Variable, Factor, Subset, Test, Transformation) |>
  pivot_wider(names_from = c(Factor, Subset), values_from = c(Test, Transformation), names_sep = " | ")

omni |> count(Test) |> tab(caption = "Number of comparisons per test")

omni |>
  mutate(p = fmt_p(p), p_BH = fmt_p(p_BH)) |>
  select(Variable, Factor, Subset, Test, Transformation, n, Statistic, Value, df1, df2, p, p_BH, Effect_size) |>
  tab(caption = "Omnibus tests")

# Descriptive statistics ----

desc_fun <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x); m <- mean(x); s <- sd(x); se <- s / sqrt(n)
  tq <- qt(0.975, n - 1)
  tibble(n = n, Mean = m, SD = s, SE = se, CI95_low = m - tq * se, CI95_high = m + tq * se,
         Median = median(x), Q1 = unname(quantile(x, .25)), Q3 = unname(quantile(x, .75)))
}

long <- dat |>
  pivot_longer(all_of(names(pools)), names_to = "var", values_to = "value") |>
  mutate(Variable = factor(pools[var], levels = unname(pools)))

desc_site <- long |> group_by(Variable, Site, Land_cover) |> reframe(desc_fun(value)) |>
  mutate(Site = as.character(Site))
desc_all <- long |> group_by(Variable, Land_cover) |> reframe(desc_fun(value)) |>
  mutate(Site = "All sites")
desc <- bind_rows(desc_site, desc_all) |> arrange(Variable, Site, Land_cover)
head(desc)

# QQ-plots of the residuals ----

qq_dat <- map_dfr(names(pools), \(y) {
  d <- filter(dat, !is.na(.data[[y]]))
  tibble(Variable = pools[[y]], res = residuals(aov(d[[y]] ~ d$Land_cover)))
}) |> mutate(Variable = factor(Variable, levels = unname(pools)))

fig_qq <- ggplot(qq_dat, aes(sample = res)) +
  stat_qq(alpha = 0.5, size = 1.2, colour = "#3C5488") +
  stat_qq_line(colour = "#E64B35") +
  facet_wrap(~Variable, scales = "free_y", nrow = 2) +
  labs(x = "Theoretical quantiles", y = "Residuals (one-way ANOVA on land cover)")

ggsave(file.path(out_dir, "Fig_QQ_residuals.png"), fig_qq, width = 14, height = 8,
       dpi = 600, device = ragg::agg_png, bg = "white")
fig_qq

# Tables ----

assump_out <- assump |>
  transmute(Variable, Factor, Subset, n, Transformation,
            `Shapiro-Wilk W` = round(W, 3), `Shapiro-Wilk p` = fmt_p(p_shapiro),
            `Groups non-normal` = paste0(n_groups_nonnormal, "/", n_groups),
            `Levene F` = round(F_levene, 2), `Levene df` = paste0(df1_levene, ", ", df2_levene),
            `Levene p` = fmt_p(p_levene),
            `Normal residuals` = ifelse(normal, "yes", "no"),
            `Homogeneous variances` = ifelse(homosc, "yes", "no"))

omni_out <- omni |>
  transmute(Variable, Factor, Subset, n, k, Test, Transformation, Statistic,
            Value = round(Value, 2), df1 = round(df1, 2), df2 = round(df2, 2),
            p = fmt_p(p), `p (BH-adjusted)` = fmt_p(p_BH),
            `Effect size` = round(Effect_size, 3), `Effect metric` = Effect_metric)

ph_out <- posthoc |>
  transmute(Variable, Factor, Subset, Test, Comparison, Method,
            Estimate = round(Estimate, 2), `CI95 low` = round(Lower95, 2), `CI95 high` = round(Upper95, 2),
            `p (adjusted)` = fmt_p(p_adj), p_num = p_adj)

desc_out <- desc |>
  transmute(Variable, Site, `Land cover` = Land_cover, n, Mean = round(Mean, 2), SD = round(SD, 2),
            SE = round(SE, 2), `CI95 low` = round(CI95_low, 2), `CI95 high` = round(CI95_high, 2),
            Median = round(Median, 2), Q1 = round(Q1, 2), Q3 = round(Q3, 2))

readme <- tibble(
  Sheet = c("Assumptions", "Omnibus_tests", "PostHoc", "Descriptives", "Test_used_summary"),
  Content = c("Shapiro-Wilk and Levene tests, raw and log(1 + x)",
              "Selected test, statistic, df, p, BH-adjusted p, effect size",
              "Tukey HSD, Games-Howell or Dunn-Holm comparisons (omnibus p < 0.05 only)",
              "n, mean, SD, SE, 95% CI, median and quartiles by variable, site and land cover",
              "Test used for each variable"))

sheets <- list(README = readme, Assumptions = assump_out, Omnibus_tests = omni_out,
               PostHoc = select(ph_out, -p_num), Descriptives = desc_out,
               Test_used_summary = test_summary)

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", border = "bottom", fontName = "Times New Roman", fontSize = 11)
bs <- createStyle(fontName = "Times New Roman", fontSize = 11)
for (nm in names(sheets)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, sheets[[nm]], headerStyle = hs)
  addStyle(wb, nm, bs, rows = 2:(nrow(sheets[[nm]]) + 1), cols = 1:ncol(sheets[[nm]]), gridExpand = TRUE)
  setColWidths(wb, nm, cols = 1:ncol(sheets[[nm]]), widths = "auto")
  freezePane(wb, nm, firstRow = TRUE)
}
saveWorkbook(wb, file.path(out_dir, "Table_stat_tests.xlsx"), overwrite = TRUE)
walk2(sheets, names(sheets), \(d, n) write_csv(d, file.path(out_dir, paste0("Table_", n, ".csv"))))

ft_style <- function(d, caption) {
  flextable(d) |>
    font(fontname = "Times New Roman", part = "all") |> fontsize(size = 8, part = "all") |>
    theme_booktabs(bold_header = FALSE) |> italic(part = "header") |>
    align(align = "center", part = "all") |> align(j = 1, align = "left", part = "all") |>
    autofit() |> set_caption(caption)
}
make_cap <- function(n, txt) {
  as_paragraph(
    as_chunk(paste0("Table S", n, ". "), props = fp_text(bold = TRUE, font.family = "Times New Roman", font.size = 11)),
    as_chunk(txt, props = fp_text(italic = TRUE, font.family = "Times New Roman", font.size = 11)))
}
add_tab <- function(doc, d, n, txt, note) {
  doc <- body_add_flextable(doc, ft_style(d, make_cap(n, txt)))
  doc |>
    body_add_fpar(fpar(ftext(paste0("Note. ", note),
                             fp_text(font.family = "Times New Roman", font.size = 9)))) |>
    body_add_break()
}

doc <- read_docx()
doc <- add_tab(doc, assump_out, 1, "Assumption checks of the one-way ANOVA for each carbon pool and grouping factor",
  "Shapiro-Wilk test on the ANOVA residuals and Levene (median-centred) test for homogeneity of variances, for raw and log(1 + x) data. Groups non-normal: groups with Shapiro-Wilk p < 0.05 / number of groups. SOC, soil organic carbon; AGC, aboveground carbon; BGC, belowground carbon.")
doc <- add_tab(doc, omni_out, 2, "Omnibus test used for each carbon pool and grouping factor",
  "ANOVA when residuals were normal and variances homogeneous (raw or log(1 + x)), Welch ANOVA when only the variances differed, Kruskal-Wallis otherwise. Effect size is eta-squared (ANOVA, Welch) or epsilon-squared (Kruskal-Wallis). BH, Benjamini-Hochberg adjustment within each comparison family.")
doc <- add_tab(doc, ph_out |> filter(p_num < 0.05) |> select(-p_num), 3, "Significant pairwise comparisons",
  "Tukey HSD after ANOVA, Games-Howell after Welch ANOVA, Dunn's test with Holm adjustment after Kruskal-Wallis (estimate = z). Only adjusted p < 0.05 shown; the full table is in the Excel file.")
doc <- add_tab(doc, desc_out, 4, "Carbon pools (Mg C ha-1) by site and land cover",
  "n, number of plots; SD, standard deviation; SE, standard error; CI95, 95% confidence interval of the mean (t distribution); Q1 and Q3, first and third quartile.")
doc <- body_end_section_landscape(doc)
print(doc, target = file.path(out_dir, "Table_stat_tests.docx"))

# Tests used ----

omni |>
  mutate(Used = ifelse(Test == "ANOVA" & Transformation != "none", "ANOVA (log)", Test)) |>
  count(Factor, Subset, Used) |>
  pivot_wider(names_from = Used, values_from = n, values_fill = 0) |>
  tab(caption = "Number of pools analysed with each test")

# Session information ----

sessionInfo()
