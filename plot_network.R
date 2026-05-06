##################网络图########################################
suppressPackageStartupMessages({
  library(igraph)
  library(tidyverse)
})

# =========================
# 1. 参数
# =========================
group_name <- "Treatment"   # 或 "Control，Treatment" 
base_dir <- file.path("SpiecEasi_network_results", group_name)

graph_file <- file.path(base_dir, "network.graphml")
node_file  <- file.path(base_dir, "network.node_list.tsv")
edge_file  <- file.path(base_dir, "network.edge_list.tsv")

out_pdf <- file.path(base_dir, paste0(group_name, "_network_plot_nature.pdf"))
out_png <- file.path(base_dir, paste0(group_name, "_network_plot_nature.png"))

# =========================
# 2. 读取
# =========================
g <- read_graph(graph_file, format = "graphml")
node_df <- read.table(node_file, header = TRUE, sep = "\t", check.names = FALSE)
edge_df <- read.table(edge_file, header = TRUE, sep = "\t", check.names = FALSE)

node_df$node_id <- as.character(node_df$node_id)
node_df <- node_df[match(V(g)$name, node_df$node_id), , drop = FALSE]

# =========================
# 3. 节点大小：degree
# =========================
deg_vals <- if ("degree" %in% colnames(node_df)) node_df$degree else degree(g)
deg_vals[is.na(deg_vals)] <- 0
V(g)$size <- scales::rescale(deg_vals, to = c(3.5, 11))

# =========================
# 4. 节点颜色：按 Phylum，先清理再映射
# =========================
if ("Phylum" %in% colnames(node_df)) {
  phylum <- as.character(node_df$Phylum)
  
  # 去空格
  phylum <- trimws(phylum)
  
  # 去掉常见前缀：p__, p_, P_, k__Bacteria|p__XXX 这种中的门名前缀部分
  phylum <- gsub("^p__", "", phylum)
  phylum <- gsub("^p_", "", phylum)
  phylum <- gsub("^P__", "", phylum)
  phylum <- gsub("^P_", "", phylum)
  
  # 若有类似 "k__Bacteria|p__Firmicutes" 这种格式，取最后一个层级
  phylum <- gsub(".*\\|p__", "", phylum)
  phylum <- gsub(".*\\|p_", "", phylum)
  
  # 空值统一
  phylum[is.na(phylum) | phylum == ""] <- "Unassigned"
  
  # 同义名统一
  phylum <- dplyr::recode(
    phylum,
    "Bacteroidetes" = "Bacteroidota",
    "Actinobacteria" = "Actinobacteriota",
    "Acidobacteria" = "Acidobacteriota",
    .default = phylum
  )
  
  # 只保留前6个最常见门，其余并入 Other
  phy_tab <- sort(table(phylum), decreasing = TRUE)
  top_phy <- names(phy_tab)[seq_len(min(6, length(phy_tab)))]
  phylum_plot <- ifelse(phylum %in% top_phy, phylum, "Other")
  phylum_plot <- factor(phylum_plot, levels = c(top_phy, "Other"))
  
  # Nature 风格、柔和配色
  phy_cols <- c(
    "Proteobacteria"   = "#4C78A8",
    "Bacteroidota"     = "#97328E",
    "Firmicutes"       = "#3C9696",
    "Actinobacteriota" = "#E45756",
    "Acidobacteriota"  = "#72B7B2",
    "Cyanobacteria"    = "#B279A2",
    "Other"            = "#BDBDBD",
    "Unassigned"       = "#BDBDBD"
  )
  
  # 如果有没预设颜色的门，自动补色
  miss_phy <- setdiff(levels(phylum_plot), names(phy_cols))
  if (length(miss_phy) > 0) {
    extra_cols <- setNames(grDevices::hcl.colors(length(miss_phy), "Set 2"), miss_phy)
    phy_cols <- c(phy_cols, extra_cols)
  }
  
  # 严格按清理后的门名映射颜色
  V(g)$node_group <- as.character(phylum_plot)
  V(g)$color <- unname(phy_cols[as.character(phylum_plot)])
  
  # 轻微透明，柔和一点
  V(g)$color <- scales::alpha(V(g)$color, 0.92)
  
} else if ("module" %in% colnames(node_df)) {
  md <- as.character(node_df$module)
  md[is.na(md)] <- "NA"
  V(g)$node_group <- md
  md_cols <- setNames(grDevices::hcl.colors(length(unique(md)), "Set 2"), unique(md))
  V(g)$color <- unname(md_cols[md])
  V(g)$color <- scales::alpha(V(g)$color, 0.92)
} else {
  V(g)$node_group <- "Node"
  V(g)$color <- scales::alpha("#7F7F7F", 0.92)
}




# =========================
# 5. 关键节点：黑色细边框
# =========================
if ("Node_role" %in% colnames(node_df)) {
  nr <- as.character(node_df$Node_role)
  nr[is.na(nr)] <- "Peripheral"
  key_idx <- nr %in% c("Connector", "Module hub", "Network hub")
  V(g)$frame.color <- ifelse(key_idx, "black", NA)
  V(g)$frame.width <- ifelse(key_idx, 1.0, 0)
  V(g)$size[key_idx] <- V(g)$size[key_idx] * 1.25
} else {
  V(g)$frame.color <- NA
  V(g)$frame.width <- 0
}

# =========================
# 6. 边：正负关系 + 曲线
# =========================
if ("edge_type" %in% colnames(edge_df)) {
  edge_type <- as.character(edge_df$edge_type)
  edge_type[is.na(edge_type)] <- "Unknown"
  
  # 柔和配色：正相关暖灰橙，负相关蓝灰
  edge_cols <- c(
    "Positive" = "#C98A4B",
    "Negative" = "#6A9FB5",
    "Unknown"  = "#C7C7C7"
  )
  
  E(g)$color <- edge_cols[edge_type]
} else {
  E(g)$color <- "#C7C7C7"
}

if ("weight" %in% colnames(edge_df)) {
  wt <- edge_df$weight
  if (all(is.na(wt)) || length(unique(na.omit(wt))) <= 1) {
    E(g)$width <- 0.8
  } else {
    wt[is.na(wt)] <- min(wt, na.rm = TRUE)
    E(g)$width <- scales::rescale(wt, to = c(0.4, 1.6))
  }
} else {
  E(g)$width <- 0.8
}

# 曲线强度
E(g)$curved <- 0.18

# =========================
# 7. 布局
# =========================
set.seed(123)
lay <- layout_with_fr(g, niter = 2000)

# =========================
# 8. 作图
# =========================
pdf(out_pdf, width = 8.2, height = 7.2)
par(mar = c(1, 1, 3, 1))

plot(
  g,
  layout = lay,
  vertex.label = NA,
  vertex.color = V(g)$color,
  vertex.size = V(g)$size,
  vertex.frame.color = V(g)$frame.color,
  vertex.frame.width = V(g)$frame.width,
  edge.color = E(g)$color,
  edge.width = E(g)$width,
  edge.curved = E(g)$curved,
  main = paste0(group_name, " network")
)

legend(
  "topleft",
  legend = names(table(V(g)$node_group)),
  pch = 21,
  pt.bg = unique(V(g)$color[match(names(table(V(g)$node_group)), V(g)$node_group)]),
  pt.cex = 1.4,
  bty = "n",
  cex = 0.8,
  title = "Node color"
)

legend(
  "bottomleft",
  legend = c("Positive edge", "Negative edge", "Key node"),
  lty = c(1, 1, NA),
  lwd = c(1.5, 1.5, NA),
  col = c("#C98A4B", "#6A9FB5", NA),
  pch = c(NA, NA, 21),
  pt.bg = c(NA, NA, "white"),
  pt.cex = c(NA, NA, 1.3),
  pt.lwd = c(NA, NA, 1),
  bty = "n",
  cex = 0.85
)
dev.off()

png(out_png, width = 2600, height = 2300, res = 300)
par(mar = c(1, 1, 3, 1))

plot(
  g,
  layout = lay,
  vertex.label = NA,
  vertex.color = V(g)$color,
  vertex.size = V(g)$size,
  vertex.frame.color = V(g)$frame.color,
  vertex.frame.width = V(g)$frame.width,
  edge.color = E(g)$color,
  edge.width = E(g)$width,
  edge.curved = E(g)$curved,
  main = paste0(group_name, " network")
)

legend(
  "topleft",
  legend = names(table(V(g)$node_group)),
  pch = 21,
  pt.bg = unique(V(g)$color[match(names(table(V(g)$node_group)), V(g)$node_group)]),
  pt.cex = 1.4,
  bty = "n",
  cex = 0.9,
  title = "Node color"
)

legend(
  "bottomleft",
  legend = c("Positive edge", "Negative edge", "Key node"),
  lty = c(1, 1, NA),
  lwd = c(1.5, 1.5, NA),
  col = c("#C98A4B", "#6A9FB5", NA),
  pch = c(NA, NA, 21),
  pt.bg = c(NA, NA, "white"),
  pt.cex = c(NA, NA, 1.3),
  pt.lwd = c(NA, NA, 1),
  bty = "n",
  cex = 0.9
)
dev.off()







##############Zi–Pi 二维关键节点图######################################
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

# =========================================================
# 1. 参数
# =========================================================
group_name <- "Control"   # 改成 "Control Treatment" 也可以
base_dir <- file.path("SpiecEasi_network_results", group_name)

node_file <- file.path(base_dir, "network.node_list.tsv")

out_pdf <- file.path(base_dir, paste0(group_name, "_ZiPi_plot_phylumcolor.pdf"))
out_png <- file.path(base_dir, paste0(group_name, "_ZiPi_plot_phylumcolor.png"))

# =========================================================
# 2. 读入
# =========================================================
plot_df <- read.table(node_file, header = TRUE, sep = "\t", check.names = FALSE)

# 如果没有 Node_role，就自动分类
if (!("Node_role" %in% colnames(plot_df))) {
  plot_df <- plot_df %>%
    mutate(
      Node_role = case_when(
        Zi <= 2.5 & Pi <= 0.62 ~ "Peripheral",
        Zi <= 2.5 & Pi > 0.62  ~ "Connector",
        Zi > 2.5  & Pi <= 0.62 ~ "Module hub",
        Zi > 2.5  & Pi > 0.62  ~ "Network hub",
        TRUE ~ "Unclassified"
      )
    )
}

# =========================================================
# 3. 清理门名，并建立和网络图一致的颜色
# =========================================================
if ("Phylum" %in% colnames(plot_df)) {
  phylum <- as.character(plot_df$Phylum)
  phylum <- trimws(phylum)
  phylum <- gsub("^p__", "", phylum)
  phylum <- gsub("^p_", "", phylum)
  phylum <- gsub("^P__", "", phylum)
  phylum <- gsub("^P_", "", phylum)
  phylum <- gsub(".*\\|p__", "", phylum)
  phylum <- gsub(".*\\|p_", "", phylum)
  
  phylum[is.na(phylum) | phylum == ""] <- "Unassigned"
  
  phylum <- dplyr::recode(
    phylum,
    "Bacteroidetes" = "Bacteroidota",
    "Actinobacteria" = "Actinobacteriota",
    "Acidobacteria" = "Acidobacteriota",
    .default = phylum
  )
  
  # 只保留前6个最常见门，其余并到 Other
  phy_tab <- sort(table(phylum), decreasing = TRUE)
  top_phy <- names(phy_tab)[seq_len(min(6, length(phy_tab)))]
  phylum_plot <- ifelse(phylum %in% top_phy, phylum, "Other")
  
  phy_cols <- c(
    "Proteobacteria"   = "#4C78A8",
    "Bacteroidota"     = "#97328E",
    "Firmicutes"       = "#3C9696",
    "Actinobacteriota" = "#E45756",
    "Acidobacteriota"  = "#72B7B2",
    "Cyanobacteria"    = "#B279A2",
    "Other"            = "#BDBDBD",
    "Unassigned"       = "#BDBDBD"
  )
  
  miss_phy <- setdiff(unique(phylum_plot), names(phy_cols))
  if (length(miss_phy) > 0) {
    extra_cols <- setNames(grDevices::hcl.colors(length(miss_phy), "Set 2"), miss_phy)
    phy_cols <- c(phy_cols, extra_cols)
  }
  
  plot_df$Phylum_clean <- phylum_plot
} else {
  plot_df$Phylum_clean <- "Other"
  phy_cols <- c("Other" = "#BDBDBD")
}

# =========================================================
# 4. 颜色与标签逻辑
# Peripheral 统一灰色；关键节点按门着色
# =========================================================
plot_df <- plot_df %>%
  mutate(
    color_group = ifelse(
      Node_role %in% c("Connector", "Module hub", "Network hub"),
      Phylum_clean,
      "Peripheral"
    ),
    label = ifelse(
      Node_role %in% c("Connector", "Module hub", "Network hub"),
      node_id,
      ""
    )
  )

# 如果有 Genus，用 Genus 当标签更友好
if ("Genus" %in% colnames(plot_df)) {
  plot_df$label <- ifelse(
    plot_df$Node_role %in% c("Connector", "Module hub", "Network hub") &
      !is.na(plot_df$Genus) & plot_df$Genus != "",
    plot_df$Genus,
    plot_df$label
  )
}

# 加上 Peripheral 颜色
all_cols <- c("Peripheral" = "grey75", phy_cols)

# 用形状区分节点角色
plot_df$Node_role2 <- factor(
  plot_df$Node_role,
  levels = c("Peripheral", "Connector", "Module hub", "Network hub", "Unclassified")
)

shape_vals <- c(
  "Peripheral" = 16,
  "Connector" = 16,
  "Module hub" = 17,
  "Network hub" = 15,
  "Unclassified" = 16
)

# =========================================================
# 5. 作图
# =========================================================
p <- ggplot(plot_df, aes(x = Pi, y = Zi)) +
  geom_point(
    aes(color = color_group, shape = Node_role2),
    size = 2.8,
    alpha = 0.9
  ) +
  geom_vline(xintercept = 0.62, linetype = 2, color = "black") +
  geom_hline(yintercept = 2.5, linetype = 2, color = "black") +
  ggrepel::geom_text_repel(
    data = subset(plot_df, Node_role %in% c("Connector", "Module hub", "Network hub")),
    aes(label = label, color = color_group),
    size = 3.1,
    max.overlaps = 30,
    box.padding = 0.35,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  scale_color_manual(values = all_cols) +
  scale_shape_manual(values = shape_vals) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 10.5)
  ) +
  labs(
    title = paste0(group_name, " Zi–Pi plot"),
    x = "Among-module connectivity (Pi)",
    y = "Within-module connectivity (Zi)"
  )

ggsave(out_pdf, p, width = 6.4, height = 5.4)
ggsave(out_png, p, width = 6.4, height = 5.4, dpi = 300)













