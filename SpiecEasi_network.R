suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(SpiecEasi)
})

# =========================================================
# 0. 参数设置
# =========================================================
otu_file <- "otu.txt"            # 行=OTU, 列=样本
meta_file <- "metadata.txt"      # 包含 SampleID, Group
tax_file <- "taxonomy.txt"       # 包含 OTUID 及分类列
outdir <- "SpiecEasi_network_results"
dir.create(outdir, showWarnings = FALSE)

# 过滤阈值
prev_cutoff <- 0.2
abund_cutoff <- 0.0005
min_sample_sum <- 1000

# SPIEC-EASI 参数
spiec_method <- "mb"
nlambda <- 20
lambda_min_ratio <- 1e-2
pulsar_rep <- 50
set.seed(123)

# =========================================================
# 1. 读取数据
# =========================================================
otu_df <- read.table(
  otu_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

meta <- read.table(
  meta_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

tax <- read.table(
  tax_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

otu_df <- as.data.frame(otu_df)
otu_df[] <- lapply(otu_df, as.numeric)

stopifnot(all(c("SampleID", "Group") %in% colnames(meta)))

if (!("OTUID" %in% colnames(tax))) {
  colnames(tax)[1] <- "OTUID"
}
tax$OTUID <- as.character(tax$OTUID)

# =========================================================
# 2. 样本对齐
# =========================================================
common_samples <- intersect(colnames(otu_df), meta$SampleID)
if (length(common_samples) < 4) {
  stop("OTU表和metadata交集样本太少。")
}

otu_df <- otu_df[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$SampleID), , drop = FALSE]

# 样本 × OTU
otu_mat_samples <- t(as.matrix(otu_df))

# =========================================================
# 3. 工具函数
# =========================================================
calc_geodesic_efficiency <- function(g) {
  if (vcount(g) < 2) return(NA_real_)
  d <- distances(g, mode = "all")
  d[is.infinite(d)] <- NA
  diag(d) <- NA
  inv_d <- 1 / d
  inv_d[is.na(inv_d)] <- 0
  sum(inv_d) / (vcount(g) * (vcount(g) - 1))
}

calc_harmonic_geodesic_distance <- function(g) {
  if (vcount(g) < 2) return(NA_real_)
  d <- distances(g, mode = "all")
  d[is.infinite(d)] <- NA
  diag(d) <- NA
  mean(d, na.rm = TRUE)
}

calc_connectedness <- function(g) {
  if (vcount(g) < 2) return(NA_real_)
  comp <- components(g)
  max(comp$csize) / vcount(g)
}

calc_eigen_centralization <- function(g) {
  if (vcount(g) < 2 || ecount(g) == 0) return(NA_real_)
  ec <- eigen_centrality(g, directed = FALSE, scale = TRUE)$vector
  max_ec <- max(ec, na.rm = TRUE)
  sum(max_ec - ec, na.rm = TRUE)
}

prepare_taxonomy <- function(tax_df, otu_ids) {
  tax_df <- as.data.frame(tax_df, stringsAsFactors = FALSE)
  tax_df$OTUID <- as.character(tax_df$OTUID)
  tax_sub <- tax_df[match(otu_ids, tax_df$OTUID), , drop = FALSE]
  if (nrow(tax_sub) != length(otu_ids)) {
    tax_sub <- data.frame(OTUID = otu_ids, stringsAsFactors = FALSE)
  } else {
    rownames(tax_sub) <- otu_ids
  }
  tax_sub
}

calc_zipi <- function(g, module_membership) {
  nodes <- V(g)$name
  mods <- module_membership[nodes]
  
  adj_mat <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
  
  Zi <- rep(NA_real_, length(nodes))
  Pi <- rep(NA_real_, length(nodes))
  
  for (i in seq_along(nodes)) {
    node <- nodes[i]
    node_mod <- mods[node]
    
    ki_all <- sum(adj_mat[node, ], na.rm = TRUE)
    
    if (ki_all == 0) {
      Zi[i] <- 0
      Pi[i] <- 0
      next
    }
    
    same_mod_nodes <- names(mods)[mods == node_mod]
    ki_in <- sum(adj_mat[node, same_mod_nodes], na.rm = TRUE)
    
    mod_levels <- unique(mods)
    ki_s <- sapply(mod_levels, function(m) {
      mod_nodes <- names(mods)[mods == m]
      sum(adj_mat[node, mod_nodes], na.rm = TRUE)
    })
    Pi[i] <- 1 - sum((ki_s / ki_all)^2, na.rm = TRUE)
    
    kin_mod <- sapply(same_mod_nodes, function(nm) {
      sum(adj_mat[nm, same_mod_nodes], na.rm = TRUE)
    })
    kin_mean <- mean(kin_mod, na.rm = TRUE)
    kin_sd <- sd(kin_mod, na.rm = TRUE)
    
    if (is.na(kin_sd) || kin_sd == 0) {
      Zi[i] <- 0
    } else {
      Zi[i] <- (ki_in - kin_mean) / kin_sd
    }
  }
  
  data.frame(
    node_id = nodes,
    Zi = Zi,
    Pi = Pi,
    stringsAsFactors = FALSE
  )
}

# =========================================================
# 4. 单组建网函数
# =========================================================
build_network_one_group <- function(group_name, otu_mat_samples, meta, tax, outdir) {
  
  message("======== Building network for: ", group_name, " ========")
  
  sample_ids <- meta$SampleID[meta$Group == group_name]
  group_mat <- otu_mat_samples[rownames(otu_mat_samples) %in% sample_ids, , drop = FALSE]
  
  if (nrow(group_mat) < 5) {
    stop("组 ", group_name, " 的样本数过少（<5），不适合建网。")
  }
  
  group_mat <- group_mat[rowSums(group_mat, na.rm = TRUE) >= min_sample_sum, , drop = FALSE]
  
  if (nrow(group_mat) < 5) {
    stop("组 ", group_name, " 过滤低测序深度后样本数过少。")
  }
  
  group_mat <- group_mat[, colSums(group_mat, na.rm = TRUE) > 0, drop = FALSE]
  
  prev_prop <- apply(group_mat, 2, function(x) sum(x > 0, na.rm = TRUE) / nrow(group_mat))
  total_relative_abundance <- colSums(group_mat, na.rm = TRUE) /
    sum(colSums(group_mat, na.rm = TRUE))
  
  keep_otu <- prev_prop >= prev_cutoff & total_relative_abundance >= abund_cutoff
  group_mat_filt <- group_mat[, keep_otu, drop = FALSE]
  
  if (ncol(group_mat_filt) < 5) {
    stop("组 ", group_name, " 过滤后 OTU 数过少（<5），无法稳定建网。")
  }
  
  otu_ids <- colnames(group_mat_filt)
  
  # --------------------------
  # SPIEC-EASI
  # --------------------------
  se <- spiec.easi(
    group_mat_filt,
    method = spiec_method,
    lambda.min.ratio = lambda_min_ratio,
    nlambda = nlambda,
    pulsar.params = list(rep.num = pulsar_rep)
  )
  
  adj <- as.matrix(getRefit(se))
  if (is.null(rownames(adj))) rownames(adj) <- otu_ids
  if (is.null(colnames(adj))) colnames(adj) <- otu_ids
  
  if (sum(adj) == 0) {
    warning("组 ", group_name, " 未检测到网络边。")
    g <- make_empty_graph(n = ncol(group_mat_filt), directed = FALSE)
    V(g)$name <- otu_ids
  } else {
    g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
    if (is.null(V(g)$name) || length(V(g)$name) == 0) {
      V(g)$name <- otu_ids
    }
  }
  
  # --------------------------
  # 边权
  # --------------------------
  edge_df <- as_data_frame(g, what = "edges")
  beta_mat <- try(symBeta(getOptBeta(se), mode = "maxabs"), silent = TRUE)
  
  if (!inherits(beta_mat, "try-error")) {
    beta_mat <- as.matrix(beta_mat)
    if (is.null(rownames(beta_mat))) rownames(beta_mat) <- otu_ids
    if (is.null(colnames(beta_mat))) colnames(beta_mat) <- otu_ids
    
    if (nrow(edge_df) > 0) {
      edge_df$weight_signed <- purrr::map2_dbl(
        edge_df$from, edge_df$to,
        function(.x, .y) {
          if (.x %in% rownames(beta_mat) && .y %in% colnames(beta_mat)) {
            beta_mat[.x, .y]
          } else {
            NA_real_
          }
        }
      )
      edge_df$weight <- abs(edge_df$weight_signed)
      edge_df$edge_type <- ifelse(edge_df$weight_signed > 0, "Positive", "Negative")
    } else {
      edge_df$weight_signed <- numeric(0)
      edge_df$weight <- numeric(0)
      edge_df$edge_type <- character(0)
    }
  } else {
    edge_df$weight_signed <- rep(NA_real_, nrow(edge_df))
    edge_df$weight <- rep(NA_real_, nrow(edge_df))
    edge_df$edge_type <- rep(NA_character_, nrow(edge_df))
  }
  
  if (ecount(g) > 0) {
    E(g)$weight_signed <- edge_df$weight_signed
    E(g)$weight <- edge_df$weight
    E(g)$edge_type <- edge_df$edge_type
  }
  
  # --------------------------
  # taxonomy 注释
  # --------------------------
  tax_sub <- prepare_taxonomy(tax, V(g)$name)
  
  rank_names <- intersect(
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    colnames(tax_sub)
  )
  
  for (rk in rank_names) {
    vals <- tax_sub[[rk]]
    if (is.null(vals)) vals <- rep(NA_character_, vcount(g))
    vals <- as.character(vals)
    if (length(vals) != vcount(g)) vals <- rep(NA_character_, vcount(g))
    igraph::vertex_attr(g, rk) <- vals
  }
  
  # --------------------------
  # 节点表
  # --------------------------
  node_ids <- V(g)$name
  if (is.null(node_ids) || length(node_ids) == 0) {
    node_ids <- otu_ids
  }
  
  deg_vals <- degree(g)
  bet_vals <- betweenness(g, directed = FALSE, normalized = TRUE)
  close_vals <- suppressWarnings(closeness(g, normalized = TRUE))
  eig_vals <- if (ecount(g) > 0) eigen_centrality(g, directed = FALSE)$vector else rep(NA_real_, vcount(g))
  
  if (length(node_ids) != vcount(g)) node_ids <- paste0("Node_", seq_len(vcount(g)))
  if (length(deg_vals) != vcount(g)) deg_vals <- rep(NA_real_, vcount(g))
  if (length(bet_vals) != vcount(g)) bet_vals <- rep(NA_real_, vcount(g))
  if (length(close_vals) != vcount(g)) close_vals <- rep(NA_real_, vcount(g))
  if (length(eig_vals) != vcount(g)) eig_vals <- rep(NA_real_, vcount(g))
  
  node_list <- data.frame(
    node_id = node_ids,
    degree = deg_vals,
    betweenness = bet_vals,
    closeness = close_vals,
    eigenvector = eig_vals,
    stringsAsFactors = FALSE
  )
  
  keep_tax_cols <- setdiff(colnames(tax_sub), "OTUID")
  if (length(keep_tax_cols) > 0) {
    node_list <- cbind(node_list, tax_sub[, keep_tax_cols, drop = FALSE])
  }
  
  # --------------------------
  # 模块
  # --------------------------
  if (ecount(g) > 0) {
    cl <- cluster_louvain(g)
    modules_num <- length(unique(membership(cl)))
    modularity_value <- modularity(cl)
    
    module_df <- data.frame(
      node_id = names(membership(cl)),
      module = as.integer(membership(cl)),
      stringsAsFactors = FALSE
    )
    
    zipi_df <- calc_zipi(g, membership(cl)) %>%
      dplyr::mutate(
        Node_role = dplyr::case_when(
          Zi <= 2.5 & Pi <= 0.62 ~ "Peripheral",
          Zi <= 2.5 & Pi > 0.62  ~ "Connector",
          Zi > 2.5  & Pi <= 0.62 ~ "Module hub",
          Zi > 2.5  & Pi > 0.62  ~ "Network hub",
          TRUE ~ "Unclassified"
        )
      )
  } else {
    modules_num <- 0
    modularity_value <- NA_real_
    
    module_df <- data.frame(
      node_id = V(g)$name,
      module = NA_integer_,
      stringsAsFactors = FALSE
    )
    
    zipi_df <- data.frame(
      node_id = V(g)$name,
      Zi = NA_real_,
      Pi = NA_real_,
      Node_role = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  node_list <- node_list %>%
    dplyr::left_join(module_df, by = "node_id") %>%
    dplyr::left_join(zipi_df, by = "node_id")
  
  # --------------------------
  # 边表
  # --------------------------
  edge_list <- if (nrow(edge_df) > 0) {
    data.frame(
      edge_id = paste0("edge_", seq_len(nrow(edge_df))),
      source = edge_df$from,
      target = edge_df$to,
      weight_signed = edge_df$weight_signed,
      weight = edge_df$weight,
      edge_type = edge_df$edge_type,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      edge_id = character(0),
      source = character(0),
      target = character(0),
      weight_signed = numeric(0),
      weight = numeric(0),
      edge_type = character(0)
    )
  }
  
  # --------------------------
  # 网络参数
  # --------------------------
  nodes_num <- vcount(g)
  edges_num <- ecount(g)
  
  positive_num <- if (ecount(g) > 0) sum(E(g)$edge_type == "Positive", na.rm = TRUE) else 0
  negative_num <- if (ecount(g) > 0) sum(E(g)$edge_type == "Negative", na.rm = TRUE) else 0
  positive_prop <- if (edges_num > 0) positive_num / edges_num else NA_real_
  negative_prop <- if (edges_num > 0) negative_num / edges_num else NA_real_
  
  comp <- components(g)
  largest_comp_id <- if (length(comp$csize) > 0) which.max(comp$csize) else 1
  g_lcc <- induced_subgraph(g, vids = V(g)[comp$membership == largest_comp_id])
  
  average_degree <- mean(degree(g))
  avgCC <- transitivity(g, type = "average")
  avg_path <- if (ecount(g_lcc) > 0 && vcount(g_lcc) > 1) mean_distance(g_lcc, directed = FALSE) else NA_real_
  geo_eff <- calc_geodesic_efficiency(g)
  harm_geo <- calc_harmonic_geodesic_distance(g)
  centralization_eigen <- calc_eigen_centralization(g)
  dens <- edge_density(g, loops = FALSE)
  trans <- transitivity(g, type = "global")
  con <- calc_connectedness(g)
  
  if (nodes_num > 5) {
    deg_tab <- table(degree(g))
    deg_vals2 <- as.numeric(names(deg_tab))
    freq_vals <- as.numeric(deg_tab)
    valid <- deg_vals2 > 0 & freq_vals > 0
    
    if (sum(valid) >= 2) {
      fit_df <- data.frame(
        x = log10(deg_vals2[valid]),
        y = log10(freq_vals[valid])
      )
      lm_fit <- lm(y ~ x, data = fit_df)
      r2_powerlaw <- summary(lm_fit)$r.squared
    } else {
      r2_powerlaw <- NA_real_
    }
  } else {
    r2_powerlaw <- NA_real_
  }
  
  network_parameter <- data.frame(
    Group = group_name,
    Numbers_of_OTUs = ncol(group_mat_filt),
    Total_nodes = nodes_num,
    Total_links = edges_num,
    Modules = modules_num,
    Modularity = modularity_value,
    R_square_of_power_law = r2_powerlaw,
    Average_degree_avgK = average_degree,
    Average_clustering_coefficient_avgCC = avgCC,
    Average_path_distance_GD = avg_path,
    Geodesic_efficiency_E = geo_eff,
    Harmonic_geodesic_distance_HD = harm_geo,
    Centralization_of_eigenvector_centrality_CE = centralization_eigen,
    Density_D = dens,
    Transitivity_Trans = trans,
    Connectedness_Con = con,
    Positive_edges = positive_num,
    Positive_edge_prop = positive_prop,
    Negative_edges = negative_num,
    Negative_edge_prop = negative_prop,
    Largest_component_nodes = vcount(g_lcc),
    stringsAsFactors = FALSE
  )
  
  # --------------------------
  # 关键节点表
  # --------------------------
  top_degree <- node_list %>%
    dplyr::arrange(dplyr::desc(degree)) %>%
    dplyr::slice_head(n = 20)
  
  top_betweenness <- node_list %>%
    dplyr::arrange(dplyr::desc(betweenness)) %>%
    dplyr::slice_head(n = 20)
  
  top_eigenvector <- node_list %>%
    dplyr::arrange(dplyr::desc(eigenvector)) %>%
    dplyr::slice_head(n = 20)
  
  key_nodes <- node_list %>%
    dplyr::filter(Node_role %in% c("Connector", "Module hub", "Network hub")) %>%
    dplyr::arrange(
      dplyr::desc(Node_role),
      dplyr::desc(degree),
      dplyr::desc(betweenness)
    )
  
  # --------------------------
  # 输出
  # --------------------------
  group_dir <- file.path(outdir, group_name)
  dir.create(group_dir, showWarnings = FALSE)
  
  write.table(
    node_list,
    file = file.path(group_dir, "network.node_list.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    edge_list,
    file = file.path(group_dir, "network.edge_list.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    network_parameter,
    file = file.path(group_dir, "network_parameter.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    module_df,
    file = file.path(group_dir, "network_modules.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    top_degree,
    file = file.path(group_dir, "top20_degree_nodes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    top_betweenness,
    file = file.path(group_dir, "top20_betweenness_nodes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    top_eigenvector,
    file = file.path(group_dir, "top20_eigenvector_nodes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    key_nodes,
    file = file.path(group_dir, "keystone_nodes_by_ZiPi.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write_graph(
    g,
    file.path(group_dir, "network.graphml"),
    format = "graphml"
  )
  
  return(network_parameter)
}

# =========================================================
# 5. Treatment / Control 分别建网
# =========================================================
target_groups <- c("Treatment", "Control")
target_groups <- intersect(target_groups, unique(meta$Group))

if (length(target_groups) < 2) {
  stop("metadata中的 Group 不包含 Treatment 和 Control。")
}

all_params <- lapply(
  target_groups,
  build_network_one_group,
  otu_mat_samples = otu_mat_samples,
  meta = meta,
  tax = tax,
  outdir = outdir
) %>%
  dplyr::bind_rows()

write.table(
  all_params,
  file = file.path(outdir, "network_parameter_all_groups.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat(
  "\n[完成]\n",
  "输出目录: ", outdir, "\n",
  "Treatment 和 Control 已分别建网，并输出关键节点结果。\n"
)




