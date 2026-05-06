##########################################################################
# 一次性安装（如已装可跳过）
#install.packages(c("tidyverse","data.table","lme4","lmerTest","glmmTMB","emmeans","phyloseq","vegan"))
## ==============================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(vegan)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

## ---- paths & outdir ----
otu_path  <- "otu.txt"
meta_path <- "metadata.txt"
tree_path <- if (file.exists("TREE.nwk")) "TREE.nwk" else if (file.exists("tree.txt")) "tree.txt" else NA
outdir <- "LMM_GLMM_outputs"; dir.create(outdir, showWarnings = FALSE)

## ---- 1) 读入 + 对齐（样本名交集 keep）----
otu_df <- read.table(otu_path, header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE, quote = "", comment.char = "")
otu_mat <- as.matrix(otu_df)

meta <- read.table(meta_path, header = TRUE, sep = "\t",
                   check.names = FALSE, quote = "", comment.char = "")
stopifnot(all(c("SampleID","Group") %in% names(meta)))
rownames(meta) <- meta$SampleID

keep <- intersect(colnames(otu_mat), rownames(meta))
stopifnot(length(keep) >= 2)
otu_mat <- otu_mat[, keep, drop = FALSE]
meta    <- meta[keep, , drop = FALSE]

meta$Group <- factor(meta$Group)
if ("Control" %in% levels(meta$Group)) meta$Group <- relevel(meta$Group, ref = "Control")

## ---- 2) α 多样性（不含 Chao1/ACE；可选 PD）----
X <- t(otu_mat)  # 行=样本，列=OTU

# Good's coverage
goods_cov <- function(v){
  v <- as.numeric(v)
  N <- sum(v)
  if (N == 0) return(NA_real_)
  F1 <- sum(v == 1)
  1 - F1 / N
}

alpha_tbl <- tibble::tibble(
  SampleID   = rownames(X),
  Observed   = rowSums(X > 0),
  Shannon    = vegan::diversity(X, index = "shannon"),
  Simpson    = vegan::diversity(X, index = "simpson"),
  InvSimpson = vegan::diversity(X, index = "invsimpson"),
  Fisher     = apply(X, 1, vegan::fisher.alpha),
  Coverage   = apply(X, 1, goods_cov)
)

# Faith's PD（可选）
alpha_tbl$PD <- NA_real_
if (!is.na(tree_path) &&
    requireNamespace("picante", quietly = TRUE) &&
    requireNamespace("ape", quietly = TRUE)) {
  
  tree_txt <- readLines(tree_path, warn = FALSE)
  tree_str <- paste(tree_txt, collapse = "\n")
  if (!grepl(";$", tree_str)) tree_str <- paste0(tree_str, ";")
  
  tr <- try(ape::read.tree(text = tree_str), silent = TRUE)
  
  if (!inherits(tr, "try-error")) {
    otus_in_tree <- intersect(colnames(X), tr$tip.label)
    if (length(otus_in_tree) >= 2) {
      X_pd <- X[, otus_in_tree, drop = FALSE]
      trm  <- ape::keep.tip(tr, otus_in_tree)
      PDres <- picante::pd(X_pd, trm, include.root = FALSE)
      alpha_tbl$PD[match(rownames(X_pd), alpha_tbl$SampleID)] <- PDres$PD
    } else {
      message("PD 跳过：树与 OTU 交集 < 2。")
    }
  } else {
    message("树读取失败，PD 跳过。")
  }
} else {
  message("未检测到树或 picante/ape，PD 跳过。")
}

# 保存每样本指数
readr::write_tsv(alpha_tbl, file.path(outdir, "alpha_indices_per_sample.tsv"))

# 合并元数据
alpha_df_all <- alpha_tbl %>%
  dplyr::left_join(meta, by = c("SampleID" = "SampleID"))

## ---- 3) 逐指标建模（LM，取实验组 vs Control 的对比）----
metrics <- intersect(
  c("Observed","Shannon","Simpson","InvSimpson","Fisher","Coverage","PD"),
  colnames(alpha_df_all)
)

# 选择对比（优先 eDNA - Control）
.choice_contrast <- function(fctr){
  lv <- levels(fctr)
  if (length(lv) < 2) return(NULL)
  ctrl <- if ("Control" %in% lv) "Control" else lv[1]
  edna <- lv[grep("^eDNA$", lv, ignore.case = TRUE)]
  if (length(edna) == 0) edna <- setdiff(lv, ctrl)[1]
  list(edna = edna, ctrl = ctrl)
}

fit_one_metric <- function(y){
  dat <- alpha_df_all %>%
    dplyr::select(
      SampleID,
      Group,
      dplyr::all_of(y)
    ) %>%
    stats::na.omit()
  
  if (nrow(dat) < 4 || dplyr::n_distinct(dat$Group) < 2) return(NULL)
  
  mod <- stats::lm(
    stats::as.formula(paste0(y, " ~ Group")),
    data = dat
  )
  
  model_type <- "LM"
  
  cc <- .choice_contrast(dat$Group)
  if (is.null(cc)) return(NULL)
  
  emm <- emmeans::emmeans(mod, ~ Group)
  contr_all <- pairs(emm)
  contr_df  <- as.data.frame(contr_all)
  ci_df     <- as.data.frame(confint(contr_all, level = 0.95))
  
  want  <- paste0(cc$edna, " - ", cc$ctrl)
  want2 <- paste0(cc$ctrl, " - ", cc$edna)
  
  pick <- which(tolower(contr_df$contrast) == tolower(want))
  flip <- FALSE
  
  if (length(pick) == 0) {
    pick <- which(tolower(contr_df$contrast) == tolower(want2))
    flip <- TRUE
  }
  if (length(pick) == 0) return(NULL)
  
  est <- contr_df$estimate[pick]
  se  <- contr_df$SE[pick]
  p   <- contr_df$p.value[pick]
  df_resid <- if ("df" %in% colnames(contr_df)) contr_df$df[pick] else NA_real_
  lo  <- ci_df$lower.CL[pick]
  hi  <- ci_df$upper.CL[pick]
  
  if (flip) {
    est_new <- -est
    lo_new  <- -hi
    hi_new  <- -lo
    est <- est_new
    lo  <- lo_new
    hi  <- hi_new
  }
  
  tibble::tibble(
    metric = y,
    model = model_type,
    contrast = paste0(cc$edna, " - ", cc$ctrl),
    estimate = est,
    std.error = se,
    df_resid = df_resid,
    p.value = p,
    conf.low = lo,
    conf.high = hi
  )
}

res_alpha_models <- metrics %>%
  lapply(fit_one_metric) %>%
  Filter(Negate(is.null), .) %>%
  dplyr::bind_rows() %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(padj = stats::p.adjust(p.value, method = "BH")) %>%
  dplyr::ungroup()

readr::write_tsv(
  res_alpha_models,
  file.path(outdir, "alpha_models_by_metric.tsv")
)

## ---- 4) 阿尔法多样性森林图：P值替代星号 ----
library(readr)
library(dplyr)
library(ggplot2)

df <- readr::read_tsv(
  file.path(outdir, "alpha_models_by_metric.tsv"),
  show_col_types = FALSE
)

# 统一计算 CI 和 P 标签
df <- df %>%
  dplyr::mutate(
    sig_used = dplyr::coalesce(padj, p.value),
    tcrit = ifelse(is.na(df_resid), 1.96, stats::qt(0.975, df_resid)),
    lower = estimate - tcrit * std.error,
    upper = estimate + tcrit * std.error,
    ci_nonzero = lower * upper > 0,
    p_label = dplyr::case_when(
      is.na(sig_used) ~ "",
      sig_used >= 0.05 ~ "",
      sig_used < 0.001 & ci_nonzero ~ paste0("p=", formatC(sig_used, format = "e", digits = 1)),
      ci_nonzero ~ paste0("p=", sprintf("%.3f", sig_used)),
      TRUE ~ ""
    )
  )

# 排序
df <- df %>%
  dplyr::arrange(estimate) %>%
  dplyr::mutate(metric = factor(metric, levels = metric))

# P值标签位置
df <- df %>%
  dplyr::mutate(
    label_x = ifelse(estimate >= 0, upper, lower)
  )

# 作图
p <- ggplot(df, aes(x = estimate, y = metric)) +
  
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.22,
    size = 1.1,
    color = "black"
  ) +
  
  geom_point(
    size = 4.2,
    shape = 21,
    fill = "white",
    stroke = 1.2
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    color = "grey40",
    size = 1
  ) +
  
  geom_text(
    aes(x = label_x, label = p_label),
    nudge_y = 0.18,
    size = 4.5
  ) +
  
  theme_classic(base_size = 18) +
  theme(
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(10, 30, 10, 20)
  ) +
  
  labs(
    title = "α-diversity effects (eDNA vs Control)",
    x = "Effect size (model estimate)"
  ) +
  coord_cartesian(clip = "off")

ggsave(
  file.path(outdir, "alpha_forest_vertical_compact_pvalue.pdf"),
  p,
  width = 4.6,
  height = 5.8
)

ggsave(
  file.path(outdir, "alpha_forest_vertical_compact_pvalue.png"),
  p,
  width = 4.6,
  height = 5.8,
  dpi = 300
)





# ================== 差异类群 GLM/GLMM（负二项 + offset）==================
suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
  library(readr)
  library(data.table)
})

# ====== 参数设置 ======
tax_level <- "Family"    # 可改成 "Genus" / "Species" / "Order" 等
outdir <- "LMM_GLMM_outputs"
dir.create(outdir, showWarnings = FALSE)

# ====== 读取 taxonomy.txt，兼容中文列名 ======
tax <- tryCatch(
  readr::read_tsv(
    "taxonomy.txt",
    col_types = readr::cols(.default = readr::col_character()),
    locale = readr::locale(encoding = "UTF-8"),
    progress = FALSE
  ),
  error = function(e) {
    data.table::fread(
      "taxonomy.txt",
      sep = "\t",
      header = TRUE,
      encoding = "UTF-8",
      data.table = FALSE
    )
  }
)

names(tax) <- trimws(names(tax))
tax <- dplyr::mutate_all(tax, ~trimws(.))

# 中文列名映射
cn2en <- c(
  "界" = "Kingdom",
  "门" = "Phylum",
  "纲" = "Class",
  "目" = "Order",
  "科" = "Family",
  "属" = "Genus",
  "种" = "Species",
  "OTU" = "OTUID"
)
names(tax) <- dplyr::recode(names(tax), !!!cn2en, .default = names(tax))

# 检查关键列
req <- c("OTUID", tax_level)
missing <- setdiff(req, names(tax))
if (length(missing)) {
  stop(sprintf(
    "在 taxonomy 里找不到列：%s\n现有列：%s",
    paste(missing, collapse = ", "),
    paste(names(tax), collapse = ", ")
  ))
}

# ====== 假设以下对象已在前文存在 ======
# otu: 第一列 OTUID，其余列为样本计数
# meta: 包含 SampleID 和 Group 列
# libsize: data.frame(SampleID, libsize)

# 如果没有 libsize，这里自动从 otu 生成
if (!exists("libsize")) {
  libsize <- data.frame(
    SampleID = colnames(otu)[-1],
    libsize = colSums(as.matrix(otu[, -1, drop = FALSE])),
    stringsAsFactors = FALSE
  )
}

# 确保 Group 参考水平是 Control
if (!is.factor(meta$Group)) meta$Group <- factor(meta$Group)
if ("Control" %in% levels(meta$Group)) {
  meta$Group <- relevel(meta$Group, ref = "Control")
}

# 自动识别实验组名称（除 Control 外的另一个水平）
group_levels <- levels(meta$Group)
treat_group <- setdiff(group_levels, "Control")[1]
if (is.na(treat_group) || length(treat_group) == 0) {
  stop("无法识别实验组名称，请检查 meta$Group。")
}

# ====== 数据处理：OTU + taxonomy + meta + libsize ======
tax_lvl <- tax %>%
  dplyr::select(OTUID, !!rlang::sym(tax_level))

otu_long <- otu %>%
  tidyr::pivot_longer(
    cols = -OTUID,
    names_to = "SampleID",
    values_to = "count"
  ) %>%
  dplyr::left_join(meta, by = "SampleID") %>%
  dplyr::left_join(libsize, by = "SampleID") %>%
  dplyr::left_join(tax_lvl, by = "OTUID") %>%
  dplyr::rename(Level = !!rlang::sym(tax_level)) %>%
  dplyr::mutate(
    Level = dplyr::coalesce(Level, "Unassigned"),
    count = as.numeric(count)
  )

# 聚合到目标分类层级
otu_level <- otu_long %>%
  dplyr::group_by(Level, SampleID, Group, libsize) %>%
  dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

# ====== 逐类群拟合模型 ======
run_one_level <- function(x) {
  dat <- dplyr::filter(otu_level, Level == x)
  
  # 过滤极低丰度类群
  if (sum(dat$count, na.rm = TRUE) < 10L) return(NULL)
  if (dplyr::n_distinct(dat$Group) < 2L) return(NULL)
  
  # 负二项模型
  fit <- try(
    glmmTMB::glmmTMB(
      count ~ Group + offset(log(libsize)),
      family = glmmTMB::nbinom2,
      data = dat
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NULL)
  
  co <- summary(fit)$coefficients$cond
  
  # 找到实验组那一行系数
  group_row <- grep("^Group", rownames(co), value = TRUE)
  group_row <- setdiff(group_row, "GroupControl")
  if (length(group_row) == 0) return(NULL)
  
  group_row <- group_row[1]
  
  tibble::tibble(
    Level = x,
    group_term = group_row,
    logFC = unname(co[group_row, "Estimate"]),
    SE    = unname(co[group_row, "Std. Error"]),
    z     = unname(co[group_row, "z value"]),
    pval  = unname(co[group_row, "Pr(>|z|)"])
  )
}

res_glmm_lvl <- purrr::map_df(unique(otu_level$Level), run_one_level) %>%
  dplyr::mutate(
    padj = p.adjust(pval, method = "BH"),
    ci_low = logFC - 1.96 * SE,
    ci_high = logFC + 1.96 * SE,
    direction = dplyr::case_when(
      logFC > 0 ~ paste0("Higher in ", treat_group),
      logFC < 0 ~ "Higher in Control",
      TRUE ~ "No change"
    ),
    fold_change = exp(logFC),
    fc_low = exp(ci_low),
    fc_high = exp(ci_high),
    p_label = dplyr::case_when(
      is.na(padj) ~ "",
      padj < 0.001 ~ paste0("p=", formatC(padj, format = "e", digits = 1)),
      TRUE ~ paste0("p=", sprintf("%.3f", padj))
    )
  ) %>%
  dplyr::arrange(padj, dplyr::desc(logFC))

lvl_tag <- tolower(tax_level)

# 保存完整结果
readr::write_tsv(
  res_glmm_lvl,
  file.path(outdir, paste0(lvl_tag, "_diff_abundance_glmm_full.tsv"))
)

# ====== 只保留实验组更高的类群 ======
res_treat_high <- res_glmm_lvl %>%
  dplyr::filter(
    !is.na(padj),
    logFC > 0
  ) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(logFC)
  )

readr::write_tsv(
  res_treat_high,
  file.path(outdir, paste0(lvl_tag, "_higher_in_", treat_group, ".tsv"))
)

# ====== 取实验组中更高的前6个类群 ======
res_top6 <- res_treat_high %>%
  dplyr::slice_head(n = 6)

if (nrow(res_top6) == 0) {
  stop("没有检测到实验组中更高的差异类群，无法绘图。")
}

# 排序：从小到大，便于森林图展示
res_top6 <- res_top6 %>%
  dplyr::arrange(logFC) %>%
  dplyr::mutate(
    Level = factor(Level, levels = Level),
    label_x = ci_high
  )

readr::write_tsv(
  res_top6,
  file.path(outdir, paste0(lvl_tag, "_top6_higher_in_", treat_group, ".tsv"))
)

# ====== 作图：P值替代星号 ======
p_top6 <- ggplot(res_top6, aes(x = logFC, y = Level)) +
  
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0.22,
    size = 1.1,
    color = "black"
  ) +
  
  geom_point(
    fill = "#BF7FB9",
    size = 4.2,
    shape = 21,
    stroke = 1.2,
    color = "black"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    color = "grey40",
    size = 1
  ) +
  
  geom_text(
    aes(x = label_x, label = p_label),
    nudge_y = 0.18,
    hjust = 0,
    size = 4.2
  ) +
  
  theme_classic(base_size = 18) +
  theme(
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 30, 10, 20)
  ) +
  labs(
    title = paste0("Top 6 ", tax_level, " enriched in ", treat_group),
    x = paste0("Log fold change (", treat_group, " vs Control)")
  ) +
  coord_cartesian(clip = "off")

# ====== 输出图片 ======
ggsave(
  file.path(outdir, paste0(lvl_tag, "_top6_higher_in_", treat_group, "_pvalue.pdf")),
  p_top6,
  width = 4.8,
  height = 5.8
)

ggsave(
  file.path(outdir, paste0(lvl_tag, "_top6_higher_in_", treat_group, "_pvalue.png")),
  p_top6,
  width = 4.8,
  height = 5.8,
  dpi = 300
)


