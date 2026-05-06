suppressPackageStartupMessages({
  library(vegan)
  library(tidyverse)
  library(ggpubr)
})

# ==========================================
# 0. 参数设置
# ==========================================
set.seed(123)

otu_file  <- "otu.txt"
meta_file <- "metadata.txt"
outdir <- "Niche_Analysis_Output"
if (!dir.exists(outdir)) dir.create(outdir)

# 组颜色
group_cols <- c(
  "Treatment" = "#97328E",   # 紫色
  "Control"   = "#3c9696"    # 浅绿
)

# ==========================================
# 1. 数据读取与预处理
# ==========================================
otu_raw <- read.table(
  otu_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.table(
  meta_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

common_samples <- intersect(colnames(otu_raw), rownames(metadata))
otu_raw <- otu_raw[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# 这里只需要 Group
stopifnot("Group" %in% colnames(metadata))

# 调整因子顺序
metadata$Group <- factor(metadata$Group, levels = c("Treatment", "Control"))

# ==========================================
# 2. 核心计算：社区生态位宽度 (基于 OTU)
# ==========================================
otu_data <- t(otu_raw)  # 行=样本，列=OTU

# 稀释到最小测序深度
otu_rare <- vegan::rrarefy(otu_data, sample = min(rowSums(otu_data)))

# 标准化生态位宽度
community_niche <- data.frame(
  SampleID = rownames(otu_rare),
  Community_B = vegan::diversity(otu_rare, index = "invsimpson") / ncol(otu_rare)
) %>%
  dplyr::left_join(
    tibble::rownames_to_column(metadata, var = "SampleID"),
    by = "SampleID"
  )

# 保存每个样本原始数值
readr::write_csv(
  community_niche,
  file.path(outdir, "Community_Niche_Values.csv")
)

# ==========================================
# 3. 统计学差异：只比较 Treatment vs Control
# ==========================================
p_value <- wilcox.test(Community_B ~ Group, data = community_niche)$p.value

stats_results <- tibble::tibble(
  Comparison = "Treatment_vs_Control",
  Method = "Wilcoxon rank sum test",
  P_value = p_value,
  Significance = dplyr::case_when(
    P_value < 0.001 ~ "***",
    P_value < 0.01  ~ "**",
    P_value < 0.05  ~ "*",
    TRUE ~ "ns"
  ),
  P_label = dplyr::case_when(
    P_value < 0.001 ~ paste0("p = ", format(P_value, scientific = TRUE, digits = 2)),
    TRUE ~ paste0("p = ", signif(P_value, 3))
  )
)

print("--- 社区生态位宽度统计差异 (P-values) ---")
print(stats_results)

readr::write_csv(
  stats_results,
  file.path(outdir, "Niche_Stats_Pvalues.csv")
)

# ==========================================
# 4. 输出误差分析和范围表
# ==========================================
niche_summary <- community_niche %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(Community_B, na.rm = TRUE),
    median = median(Community_B, na.rm = TRUE),
    sd = sd(Community_B, na.rm = TRUE),
    se = sd / sqrt(n),
    ci95_low = mean - qt(0.975, df = n - 1) * se,
    ci95_high = mean + qt(0.975, df = n - 1) * se,
    min = min(Community_B, na.rm = TRUE),
    max = max(Community_B, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  niche_summary,
  file.path(outdir, "Community_Niche_Summary_Stats.csv")
)

# 宽表
niche_summary_wide <- niche_summary %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = c(n, mean, median, sd, se, ci95_low, ci95_high, min, max)
  )

readr::write_csv(
  niche_summary_wide,
  file.path(outdir, "Community_Niche_Summary_Stats_Wide.csv")
)

# ==========================================
# 5. 统一 Fig.2 风格的单图
# ==========================================
y_max <- max(community_niche$Community_B, na.rm = TRUE)
y_min <- min(community_niche$Community_B, na.rm = TRUE)
y_rng <- y_max - y_min
if (y_rng == 0) y_rng <- 0.02
y_pad <- y_rng * 0.18

p_niche_fig2 <- ggplot(community_niche, aes(x = Group, y = Community_B, fill = Group)) +
  
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    linewidth = 0.9,
    color = "black"
  ) +
  
  geom_jitter(
    width = 0.10,
    size = 2.1,
    alpha = 0.95,
    color = "black"
  ) +
  
  scale_fill_manual(values = c(
    "Treatment" = "#97328E",
    "Control" = "#3c9696"
  )) +
  
  annotate(
    "text",
    x = 1.5,
    y = y_max + y_pad * 0.55,
    label = stats_results$P_label,
    size = 4.8,
    fontface = "plain"
  ) +
  
  coord_cartesian(ylim = c(y_min, y_max + y_pad)) +
  
  labs(
    title = "Community niche breadth",
    x = NULL,
    y = "Standardized Niche Index (B)"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15, face = "plain", color = "black"),
    axis.text.x = element_text(size = 15, face = "plain", color = "black"),
    axis.text.y = element_text(size = 14, face = "plain", color = "black"),
    axis.line = element_line(linewidth = 0.9, color = "black"),
    axis.ticks = element_line(linewidth = 0.9, color = "black"),
    axis.ticks.length = unit(0.14, "cm"),
    legend.position = "none",
    plot.margin = margin(6, 6, 6, 6)
  )

# ==========================================
# 6. 输出图片
# ==========================================
ggsave(
  file.path(outdir, "Community_Niche_GroupOnly_Fig2Style.pdf"),
  p_niche_fig2,
  width = 3.2, height = 3.6
)

ggsave(
  file.path(outdir, "Community_Niche_GroupOnly_Fig2Style.png"),
  p_niche_fig2,
  width = 3.2, height = 3.6, dpi = 300
)

cat(
  "\n[分析完成]\n",
  "1. 每样本数值: ", file.path(outdir, "Community_Niche_Values.csv"), "\n",
  "2. Wilcoxon结果表: ", file.path(outdir, "Niche_Stats_Pvalues.csv"), "\n",
  "3. 误差分析表: ", file.path(outdir, "Community_Niche_Summary_Stats.csv"), "\n",
  "4. 宽格式统计表: ", file.path(outdir, "Community_Niche_Summary_Stats_Wide.csv"), "\n",
  "5. 图片: ", file.path(outdir, "Community_Niche_GroupOnly_Fig2Style.pdf"), "\n"
)




