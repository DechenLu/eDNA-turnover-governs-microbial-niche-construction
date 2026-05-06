# ==========================================
# 1. 准备工作与数据读取
# ==========================================
library(vegan)
library(ggplot2)
library(ggpubr)
library(patchwork) # 用于组合图片

beta_dir <- "Beta_Diversity_Analysis"
if(!dir.exists(beta_dir)) dir.create(beta_dir)

# 读取数据
otu_raw <- read.table("otu.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, quote = "")
meta    <- read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# 对齐样本
common_samples <- intersect(colnames(otu_raw), rownames(meta))
otu_raw <- otu_raw[, common_samples]
meta    <- meta[common_samples, , drop = FALSE]

# ==========================================
# 2. 数据抽平 (Rarefaction)
# ==========================================
otu_data <- t(otu_raw)
set.seed(123)
min_depth <- min(rowSums(otu_data))
otu_rare <- rrarefy(otu_data, sample = min_depth)

# ==========================================
# 3. 核心函数：执行特定条件下的 Beta 多样性分析
# ==========================================
analyze_beta <- function(sub_otu, sub_meta, title_suffix) {
  
  # 1. 计算 Bray-Curtis 距离
  dist_matrix <- vegdist(sub_otu, method = "bray")
  
  # 2. PERMANOVA 检验 (Control vs Treatment)
  set.seed(123)
  adonis_res <- adonis2(dist_matrix ~ Group, data = sub_meta)
  r2_val <- round(adonis_res$R2[1], 3)
  p_val  <- adonis_res$`Pr(>F)`[1]
  
  # 3. PCoA 计算
  pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  points <- as.data.frame(pcoa_res$points)
  colnames(points) <- c("PCoA1", "PCoA2")
  eig <- pcoa_res$eig
  pc1_exp <- round(eig[1] / sum(eig) * 100, 1)
  pc2_exp <- round(eig[2] / sum(eig) * 100, 1)
  
  plot_df <- cbind(points, sub_meta)
  plot_df$Group <- factor(plot_df$Group, levels = c("Control", "Treatment"))
  
  # 4. 绘图
  p <- ggplot(plot_df, aes(x = PCoA1, y = PCoA2, color = Group, shape = Group)) +
    stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.1, linetype = 0) +
    geom_point(size = 3, alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray80") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray80") +
    labs(x = paste0("PCoA1 (", pc1_exp, "%)"),
         y = paste0("PCoA2 (", pc2_exp, "%)"),
         title = paste0("Beta Diversity: ", title_suffix),
         subtitle = paste0("PERMANOVA: R² = ", r2_val, ", p = ", p_val)) +
    scale_color_manual(values = c("Control" = "#3c9696", "Treatment" = "#97328E")) +
    scale_fill_manual(values = c("Control" = "#3c9696", "Treatment" = "#97328E")) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "right")
  
  return(list(plot = p, stats = adonis_res))
}

# ==========================================
# 4. 分层执行分析 (Surface vs Liquid)
# ==========================================

# --- 4.1 表面样品 (Surface) ---
meta_surf <- subset(meta, Conditions == "Surface")
otu_surf  <- otu_rare[rownames(meta_surf), ]
res_surf  <- analyze_beta(otu_surf, meta_surf, "Surface Samples")

# --- 4.2 沉积物样品 (Sediment) ---
meta_Sedi <- subset(meta, Conditions == "Sediment")
otu_Sedi  <- otu_rare[rownames(meta_Sedi), ]
res_Sedi  <- analyze_beta(otu_Sedi, meta_Sedi, "Sediment Samples")

# ==========================================
# 5. 结果保存与打印
# ==========================================

# 组合图片并保存
final_plot <- res_surf$plot / res_Sedi$plot
ggsave(file.path(beta_dir, "PCoA_Combined_by_Conditions.pdf"), final_plot, width = 8, height = 10)
ggsave(file.path(beta_dir, "PCoA_Combined_by_Conditions.png"), final_plot, width = 8, height = 10, dpi = 300)

# 打印统计结果
cat("\n--- Surface 组 PERMANOVA 结果 ---\n")
print(res_surf$stats)
cat("\n--- Sediment 组 PERMANOVA 结果 ---\n")
print(res_liq$stats)

cat("\n分析完成！图片已存至 Beta_Diversity_Analysis 文件夹。")



