rm(list = ls())

library(R.matlab)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(RColorBrewer)
library(scales)

# ============================================================
# 0. 设置本地文件路径
# ============================================================

mat_dir <- "/Users/cj/Desktop/WM_prediction/code/barplot/barplot_pre"

fig_dir <- file.path(mat_dir, "NetworkPair_bar_plots_R_NetAvg_only_topWeighted_noZero_equalBar")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# ============================================================
# 1. 文件名设置
#    每个 .mat 文件中都只需要包含变量 Net_avg
# ============================================================

files <- list(
  age = list(
    label = "Age",
    GG = file.path(mat_dir, "datasetsAveraged_Haufe_FCmatrix_gg_networkAveraged.mat"),
    GW = file.path(mat_dir, "datasetsAveraged_Haufe_FCmatrix_gw_networkAveraged.mat"),
    WW = file.path(mat_dir, "datasetsAveraged_Haufe_FCmatrix_ww_networkAveraged.mat")
  ),
  cog = list(
    label = "Cognition",
    GG = file.path(mat_dir, "ABCDcog_Haufe_FCmatrix_gg_networkAveraged.mat"),
    GW = file.path(mat_dir, "ABCDcog_Haufe_FCmatrix_gw_networkAveraged.mat"),
    WW = file.path(mat_dir, "ABCDcog_Haufe_FCmatrix_ww_networkAveraged.mat")
  ),
  pfactor = list(
    label = "p-factor / ADHD",
    GG = file.path(mat_dir, "ABCDpfactor_Haufe_FCmatrix_gg_networkAveraged.mat"),
    GW = file.path(mat_dir, "ABCDpfactor_Haufe_FCmatrix_gw_networkAveraged.mat"),
    WW = file.path(mat_dir, "ABCDpfactor_Haufe_FCmatrix_ww_networkAveraged.mat")
  )
)

all_files <- c(
  files$age$GG, files$age$GW, files$age$WW,
  files$cog$GG, files$cog$GW, files$cog$WW,
  files$pfactor$GG, files$pfactor$GW, files$pfactor$WW
)

missing_files <- all_files[!file.exists(all_files)]
if (length(missing_files) > 0) {
  stop(paste0(
    "下面这些文件不存在，请检查 mat_dir 或文件名：\n",
    paste(missing_files, collapse = "\n")
  ))
}

# ============================================================
# 2. 网络名称和颜色
#    颜色保持不变
# ============================================================

GM_net_names <- c("VIS", "SMN", "DAN", "VAN", "LIM", "CON", "DMN")
WM_net_names <- c("BS", "COMM", "ASSOC", "PROJ", "SWM")

GM_cols <- c(
  VIS = "#A84AA7",
  SMN = "#75B4D6",
  DAN = "#008C3A",
  VAN = "#E26AF2",
  LIM = "#B8D96A",
  CON = "#FDC06B",
  DMN = "#FF7F8E"
)

WM_cols <- c(
  BS    = "#8A8A8A",
  COMM  = "#75C1D0",
  ASSOC = "#7FBF90",
  PROJ  = "#F5D1B0",
  SWM   = "#B39BC8"
)

all_cols <- c(GM_cols, WM_cols)

GM_legend_labels <- c(
  VIS = "VS",
  SMN = "SM",
  DAN = "DA",
  VAN = "VA",
  LIM = "LIM",
  CON = "FP",
  DMN = "DM"
)

WM_legend_labels <- c(
  BS = "BS",
  COMM = "COMM",
  ASSOC = "ASSOC",
  PROJ = "PROJ",
  SWM = "SWM"
)

# ============================================================
# 3. 读取 .mat 中的 Net_avg
# ============================================================

get_mat_var <- function(S, varname) {
  var_dot <- gsub("_", ".", varname)
  
  if (varname %in% names(S)) {
    return(S[[varname]])
  }
  
  if (var_dot %in% names(S)) {
    return(S[[var_dot]])
  }
  
  stop(paste0(
    "Cannot find variable: ", varname,
    "\nAvailable variables are:\n",
    paste(names(S), collapse = ", ")
  ))
}

read_net_avg <- function(file) {
  S <- readMat(file)
  Net_avg <- get_mat_var(S, "Net_avg")
  Net_avg <- as.matrix(Net_avg)
  return(Net_avg)
}

read_one_phenotype_results <- function(file_GG, file_GW, file_WW) {
  list(
    GG = list(value = read_net_avg(file_GG)),
    GW = list(value = read_net_avg(file_GW)),
    WW = list(value = read_net_avg(file_WW))
  )
}

age_res <- read_one_phenotype_results(
  file_GG = files$age$GG,
  file_GW = files$age$GW,
  file_WW = files$age$WW
)

cog_res <- read_one_phenotype_results(
  file_GG = files$cog$GG,
  file_GW = files$cog$GW,
  file_WW = files$cog$WW
)

pfactor_res <- read_one_phenotype_results(
  file_GG = files$pfactor$GG,
  file_GW = files$pfactor$GW,
  file_WW = files$pfactor$WW
)

# ============================================================
# 4. matrix -> network-pair data frame
# ============================================================

make_symmetric_pair_df <- function(mat_value,
                                   net_names,
                                   phenotype,
                                   feature_type,
                                   net1_type,
                                   net2_type,
                                   keep_upper = TRUE) {
  
  n_net <- length(net_names)
  
  if (!all(dim(mat_value) == c(n_net, n_net))) {
    stop(paste0(
      phenotype, " ", feature_type,
      ": matrix dimension does not match network names. Matrix dim = ",
      paste(dim(mat_value), collapse = " x "),
      ", expected = ", n_net, " x ", n_net
    ))
  }
  
  out <- list()
  
  for (i in seq_len(n_net)) {
    for (j in seq_len(n_net)) {
      if (keep_upper && j < i) next
      
      out[[length(out) + 1]] <- data.frame(
        phenotype = phenotype,
        feature_type = feature_type,
        net1_type = net1_type,
        net2_type = net2_type,
        net1_id = i,
        net2_id = j,
        net1_name = net_names[i],
        net2_name = net_names[j],
        pair_label = paste0(net_names[i], "-", net_names[j]),
        haufe_weight = mat_value[i, j],
        stringsAsFactors = FALSE
      )
    }
  }
  
  bind_rows(out)
}

make_rectangular_pair_df <- function(mat_value,
                                     row_names,
                                     col_names,
                                     phenotype,
                                     feature_type,
                                     net1_type,
                                     net2_type) {
  
  n_row <- length(row_names)
  n_col <- length(col_names)
  
  if (!all(dim(mat_value) == c(n_row, n_col))) {
    stop(paste0(
      phenotype, " ", feature_type,
      ": matrix dimension does not match network names. Matrix dim = ",
      paste(dim(mat_value), collapse = " x "),
      ", expected = ", n_row, " x ", n_col
    ))
  }
  
  out <- list()
  
  for (i in seq_len(n_row)) {
    for (j in seq_len(n_col)) {
      out[[length(out) + 1]] <- data.frame(
        phenotype = phenotype,
        feature_type = feature_type,
        net1_type = net1_type,
        net2_type = net2_type,
        net1_id = i,
        net2_id = j,
        net1_name = row_names[i],
        net2_name = col_names[j],
        pair_label = paste0(row_names[i], "-", col_names[j]),
        haufe_weight = mat_value[i, j],
        stringsAsFactors = FALSE
      )
    }
  }
  
  bind_rows(out)
}

make_one_phenotype_df <- function(res, phenotype) {
  
  df_GG <- make_symmetric_pair_df(
    mat_value = res$GG$value,
    net_names = GM_net_names,
    phenotype = phenotype,
    feature_type = "G-G",
    net1_type = "GM",
    net2_type = "GM",
    keep_upper = TRUE
  )
  
  df_GW <- make_rectangular_pair_df(
    mat_value = res$GW$value,
    row_names = GM_net_names,
    col_names = WM_net_names,
    phenotype = phenotype,
    feature_type = "G-W",
    net1_type = "GM",
    net2_type = "WM"
  )
  
  df_WW <- make_symmetric_pair_df(
    mat_value = res$WW$value,
    net_names = WM_net_names,
    phenotype = phenotype,
    feature_type = "W-W",
    net1_type = "WM",
    net2_type = "WM",
    keep_upper = TRUE
  )
  
  bind_rows(df_GG, df_GW, df_WW)
}

df_age <- make_one_phenotype_df(age_res, "Age")
df_cog <- make_one_phenotype_df(cog_res, "Cognition")
df_pfactor <- make_one_phenotype_df(pfactor_res, "p-factor / ADHD")

df_all <- bind_rows(df_age, df_cog, df_pfactor)

write.csv(
  df_all,
  file.path(fig_dir, "Age_Cognition_pfactor_GG_GW_WW_networkPair_barInput_NetAvgOnly.csv"),
  row.names = FALSE
)

# ============================================================
# 5. 只选择 top positive 和 top negative network-pairs
#    0 不占位；接近 0 保留
# ============================================================

select_top_pairs <- function(df,
                             n_pos = 12,
                             n_neg = 12) {
  
  df_clean <- df %>%
    filter(
      !is.na(haufe_weight),
      haufe_weight != 0
    )
  
  df_pos <- df_clean %>%
    filter(haufe_weight > 0) %>%
    arrange(desc(haufe_weight)) %>%
    slice_head(n = n_pos)
  
  df_neg <- df_clean %>%
    filter(haufe_weight < 0) %>%
    arrange(haufe_weight) %>%
    slice_head(n = n_neg) %>%
    arrange(desc(haufe_weight))
  
  bind_rows(df_pos, df_neg)
}

make_bar_segment_df <- function(df,
                                sort_by = c("signed", "absolute"),
                                n_pos = 12,
                                n_neg = 12,
                                gap_size = 0,
                                bar_half_width = 0.4,
                                total_slots = 24) {
  
  sort_by <- match.arg(sort_by)
  
  df <- select_top_pairs(
    df = df,
    n_pos = n_pos,
    n_neg = n_neg
  )
  
  if (nrow(df) == 0) {
    stop("No non-zero bars remain after top-pair selection.")
  }
  
  if (sort_by == "signed") {
    
    df_pos <- df %>%
      filter(haufe_weight > 0) %>%
      arrange(desc(haufe_weight))
    
    df_neg <- df %>%
      filter(haufe_weight < 0) %>%
      arrange(desc(haufe_weight))
    
    n_pos_real <- nrow(df_pos)
    
    if (nrow(df_pos) > 0) {
      df_pos <- df_pos %>%
        mutate(x = row_number())
    }
    
    if (nrow(df_neg) > 0) {
      df_neg <- df_neg %>%
        mutate(x = n_pos_real + gap_size + row_number())
    }
    
    df <- bind_rows(df_pos, df_neg)
    
  } else {
    
    df <- df %>%
      arrange(desc(abs(haufe_weight))) %>%
      mutate(x = row_number())
  }
  
  df <- df %>%
    mutate(
      pair_label = factor(pair_label, levels = pair_label)
    )
  
  seg1 <- df %>%
    transmute(
      phenotype,
      feature_type,
      x,
      pair_label,
      net_name = net1_name,
      net_type = net1_type,
      part = "net1",
      haufe_weight,
      ymin = 0,
      ymax = haufe_weight / 2
    )
  
  seg2 <- df %>%
    transmute(
      phenotype,
      feature_type,
      x,
      pair_label,
      net_name = net2_name,
      net_type = net2_type,
      part = "net2",
      haufe_weight,
      ymin = haufe_weight / 2,
      ymax = haufe_weight
    )
  
  bar_df <- bind_rows(seg1, seg2) %>%
    mutate(
      xmin = x - bar_half_width,
      xmax = x + bar_half_width
    )
  
  list(
    pair_df = df,
    bar_df = bar_df,
    total_slots = total_slots
  )
}

# ============================================================
# 6. WM 斜线 hatch
# ============================================================

make_hatch_df <- function(bar_df, n_lines = 3) {
  
  wm_df <- bar_df %>% filter(net_type == "WM")
  
  if (nrow(wm_df) == 0) {
    return(data.frame())
  }
  
  hatch_list <- list()
  
  for (r in seq_len(nrow(wm_df))) {
    row <- wm_df[r, ]
    
    xmin <- row$xmin
    xmax <- row$xmax
    ymin <- min(row$ymin, row$ymax)
    ymax <- max(row$ymin, row$ymax)
    
    if (!is.finite(ymin) || !is.finite(ymax) || ymin == ymax) {
      next
    }
    
    dx <- xmax - xmin
    dy <- ymax - ymin
    
    starts <- seq(xmin - dx * 0.8, xmax, length.out = n_lines)
    
    for (s in starts) {
      x1 <- max(xmin, s)
      x2 <- min(xmax, s + dx * 0.8)
      
      if (x2 <= x1) next
      
      y1 <- ymin + (x1 - s) / (dx * 0.8) * dy
      y2 <- ymin + (x2 - s) / (dx * 0.8) * dy
      
      hatch_list[[length(hatch_list) + 1]] <- data.frame(
        x = x1,
        xend = x2,
        y = y1,
        yend = y2
      )
    }
  }
  
  bind_rows(hatch_list)
}

# ============================================================
# 7. 主画图函数
# ============================================================

plot_network_pair_bar <- function(df,
                                  title_text,
                                  sort_by = "signed",
                                  y_label = "Feature weights",
                                  show_title = TRUE,
                                  hatch_color = "white",
                                  n_pos = 12,
                                  n_neg = 12,
                                  gap_size = 0,
                                  bar_half_width = 0.4,
                                  total_slots = 24,
                                  x_positions = NULL) {
  
  res <- make_bar_segment_df(
    df = df,
    sort_by = sort_by,
    n_pos = n_pos,
    n_neg = n_neg,
    gap_size = gap_size,
    bar_half_width = bar_half_width,
    total_slots = total_slots
  )
  
  pair_df <- res$pair_df
  bar_df <- res$bar_df
  total_slots <- res$total_slots
  
  # ------------------------------------------------------------
  # Optional custom x positions
  # 只改变 bar 的中心位置，不改变 bar_half_width。
  # 因此 bar 宽度保持一致，只改变 bar 之间的距离。
  # 主要用于 Cognition 的 W-W panel。
  # ------------------------------------------------------------
  if (!is.null(x_positions)) {
    
    n_bar <- nrow(pair_df)
    
    if (length(x_positions) < n_bar) {
      stop("Length of x_positions is smaller than the number of bars.")
    }
    
    x_map <- data.frame(
      x_old = pair_df$x,
      x_new = x_positions[seq_len(n_bar)]
    )
    
    pair_df <- pair_df %>%
      left_join(x_map, by = c("x" = "x_old")) %>%
      mutate(x = x_new) %>%
      select(-x_new)
    
    bar_df <- bar_df %>%
      left_join(x_map, by = c("x" = "x_old")) %>%
      mutate(
        x = x_new,
        xmin = x - bar_half_width,
        xmax = x + bar_half_width
      ) %>%
      select(-x_new)
    
    total_slots <- max(total_slots, max(x_positions, na.rm = TRUE))
  }
  
  hatch_df <- make_hatch_df(bar_df, n_lines = 3)
  
  max_val <- max(pair_df$haufe_weight, na.rm = TRUE)
  min_val <- min(pair_df$haufe_weight, na.rm = TRUE)
  yrange <- max_val - min_val
  
  if (!is.finite(yrange) || yrange == 0) {
    yrange <- max(abs(c(max_val, min_val)), na.rm = TRUE)
    if (!is.finite(yrange) || yrange == 0) yrange <- 1
  }
  
  y_top <- max_val + 0.06 * yrange
  y_bottom <- min_val - 0.06 * yrange
  
  if (min_val >= 0) {
    y_bottom <- -0.04 * max_val
  }
  
  if (max_val <= 0) {
    y_top <- 0.04 * abs(min_val)
  }
  
  p <- ggplot() +
    geom_rect(
      data = bar_df,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = net_name
      ),
      color = NA
    ) +
    scale_fill_manual(values = all_cols, drop = FALSE) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      color = "black"
    ) +
    scale_x_continuous(
      limits = c(0.5, total_slots + 0.5),
      breaks = NULL,
      labels = NULL,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = function(x) pretty(x, n = 4),
      labels = label_number(accuracy = 0.001)
    ) +
    coord_cartesian(
      ylim = c(y_bottom, y_top),
      clip = "off"
    ) +
    labs(
      x = NULL,
      y = y_label,
      title = ifelse(show_title, title_text, "")
    ) +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 15, hjust = 0.5, face = "plain"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(size = 10.5, color = "black"),
      axis.title.y = element_text(size = 12.5, color = "black"),
      axis.line.y = element_line(linewidth = 0.4, color = "black"),
      axis.ticks.y = element_line(linewidth = 0.35, color = "black"),
      axis.ticks.length.y = unit(2.5, "pt"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 8, 5, 6)
    )
  
  if (nrow(hatch_df) > 0) {
    p <- p +
      geom_segment(
        data = hatch_df,
        aes(x = x, xend = xend, y = y, yend = yend),
        color = hatch_color,
        linewidth = 0.18,
        inherit.aes = FALSE
      )
  }
  
  return(p)
}

# ============================================================
# 8. 生成 9 张图
# ============================================================

make_plot_for <- function(phenotype_name, feature_name) {
  
  df <- df_all %>%
    filter(
      phenotype == phenotype_name,
      feature_type == feature_name
    )
  
  # 所有 feature type 都使用最多 12 positive + 12 negative。
  # 0 不画，不补 fake bar；接近 0 的非零值保留。
  n_pos_use <- 12
  n_neg_use <- 12
  
  # 默认不改变 x 位置。
  x_positions_use <- NULL
  
  # ------------------------------------------------------------
  # 只调整 Cognition 的 W-W：
  # bar 宽度保持一致，只增加 bar 和 bar 之间的距离。
  # 将实际存在的 bar 均匀 spread 到 24 个 slots，
  # 使其总占据位置接近 Age 的 W-W panel。
  # ------------------------------------------------------------
  if (phenotype_name == "Cognition" && feature_name == "W-W") {
    
    df_tmp <- select_top_pairs(
      df = df,
      n_pos = n_pos_use,
      n_neg = n_neg_use
    )
    
    n_bar_tmp <- nrow(df_tmp)
    
    if (n_bar_tmp > 1) {
      x_positions_use <- seq(
        from = 1,
        to = 24,
        length.out = n_bar_tmp
      )
    } else if (n_bar_tmp == 1) {
      x_positions_use <- 12
    }
  }
  
  plot_network_pair_bar(
    df,
    title_text = paste0(phenotype_name, " | ", feature_name),
    sort_by = "signed",
    y_label = "Feature weights",
    hatch_color = "white",
    n_pos = n_pos_use,
    n_neg = n_neg_use,
    gap_size = 0,
    bar_half_width = 0.4,
    total_slots = 24,
    x_positions = x_positions_use
  )
}

p_age_GG <- make_plot_for("Age", "G-G")
p_age_GW <- make_plot_for("Age", "G-W")
p_age_WW <- make_plot_for("Age", "W-W")

p_cog_GG <- make_plot_for("Cognition", "G-G")
p_cog_GW <- make_plot_for("Cognition", "G-W")
p_cog_WW <- make_plot_for("Cognition", "W-W")

p_pfactor_GG <- make_plot_for("p-factor / ADHD", "G-G")
p_pfactor_GW <- make_plot_for("p-factor / ADHD", "G-W")
p_pfactor_WW <- make_plot_for("p-factor / ADHD", "W-W")

# ============================================================
# 9. 保存单张图
# ============================================================

plot_list <- list(
  Age_GG = p_age_GG,
  Age_GW = p_age_GW,
  Age_WW = p_age_WW,
  Cognition_GG = p_cog_GG,
  Cognition_GW = p_cog_GW,
  Cognition_WW = p_cog_WW,
  pfactor_ADHD_GG = p_pfactor_GG,
  pfactor_ADHD_GW = p_pfactor_GW,
  pfactor_ADHD_WW = p_pfactor_WW
)

# 所有单图宽度一致，避免 W-W 被单独压窄或放大
get_plot_width <- function(feature_name) {
  4.9
}

for (nm in names(plot_list)) {
  
  p <- plot_list[[nm]]
  
  feature_name <- dplyr::case_when(
    grepl("_GG$", nm) ~ "G-G",
    grepl("_GW$", nm) ~ "G-W",
    grepl("_WW$", nm) ~ "W-W",
    TRUE ~ "G-G"
  )
  
  w <- get_plot_width(feature_name)
  
  ggsave(
    filename = file.path(fig_dir, paste0(nm, "_networkPair_bar_sorted_topWeighted_noZero_equalBar.pdf")),
    plot = p,
    width = w,
    height = 4,
    device = cairo_pdf
  )
  
  ggsave(
    filename = file.path(fig_dir, paste0(nm, "_networkPair_bar_sorted_topWeighted_noZero_equalBar.png")),
    plot = p,
    width = w,
    height = 4,
    dpi = 600
  )
}

# ============================================================
# 10. 合并为 3 × 3 总图
#     三列等宽；每个 panel 内部 x 轴统一为 24 slots
# ============================================================

combined_3x3 <-
  (
    (p_age_GG | p_age_GW | p_age_WW) +
      plot_layout(widths = c(1, 1, 1))
  ) /
  (
    (p_cog_GG | p_cog_GW | p_cog_WW) +
      plot_layout(widths = c(1, 1, 1))
  ) /
  (
    (p_pfactor_GG | p_pfactor_GW | p_pfactor_WW) +
      plot_layout(widths = c(1, 1, 1))
  ) +
  plot_annotation(tag_levels = "A")

ggsave(
  filename = file.path(fig_dir, "Age_Cognition_pfactor_GG_GW_WW_networkPair_bar_3x3_topWeighted_noZero_equalBar.pdf"),
  plot = combined_3x3,
  width = 15.0,
  height = 12,
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_dir, "Age_Cognition_pfactor_GG_GW_WW_networkPair_bar_3x3_topWeighted_noZero_equalBar.png"),
  plot = combined_3x3,
  width = 15.0,
  height = 12,
  dpi = 600
)

# ============================================================
# 11. 每个 phenotype 的 1 × 3 图
# ============================================================

age_1x3 <- (p_age_GG | p_age_GW | p_age_WW) +
  plot_layout(widths = c(1, 1, 1))

cog_1x3 <- (p_cog_GG | p_cog_GW | p_cog_WW) +
  plot_layout(widths = c(1, 1, 1))

pfactor_1x3 <- (p_pfactor_GG | p_pfactor_GW | p_pfactor_WW) +
  plot_layout(widths = c(1, 1, 1))

ggsave(
  filename = file.path(fig_dir, "Age_GG_GW_WW_networkPair_bar_1x3_topWeighted_noZero_equalBar.pdf"),
  plot = age_1x3,
  width = 15.0,
  height = 4.2,
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_dir, "Age_GG_GW_WW_networkPair_bar_1x3_topWeighted_noZero_equalBar.png"),
  plot = age_1x3,
  width = 15.0,
  height = 4.2,
  dpi = 600
)

ggsave(
  filename = file.path(fig_dir, "Cognition_GG_GW_WW_networkPair_bar_1x3_topWeighted_noZero_equalBar.pdf"),
  plot = cog_1x3,
  width = 15.0,
  height = 4.2,
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_dir, "Cognition_GG_GW_WW_networkPair_bar_1x3_topWeighted_noZero_equalBar.png"),
  plot = cog_1x3,
  width = 15.0,
  height = 4.2,
  dpi = 600
)

ggsave(
  filename = file.path(fig_dir, "pfactor_ADHD_GG_GW_WW_networkPair_bar_1x3_topWeighted_noZero_equalBar.pdf"),
  plot = pfactor_1x3,
  width = 15.0,
  height = 4.2,
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_dir, "pfactor_ADHD_GG_GW_WW_networkPair_bar_1x3_topWeighted_noZero_equalBar.png"),
  plot = pfactor_1x3,
  width = 15.0,
  height = 4.2,
  dpi = 600
)

# ============================================================
# 12. 单独生成颜色图例
# ============================================================

make_color_legend_plot <- function() {
  
  gm_legend_df <- data.frame(
    group = "GM",
    name = names(GM_cols),
    label = GM_legend_labels[names(GM_cols)],
    color = unname(GM_cols),
    x = seq_along(GM_cols),
    y = 2,
    stringsAsFactors = FALSE
  )
  
  wm_legend_df <- data.frame(
    group = "WM",
    name = names(WM_cols),
    label = WM_legend_labels[names(WM_cols)],
    color = unname(WM_cols),
    x = seq_along(WM_cols),
    y = 1,
    stringsAsFactors = FALSE
  )
  
  legend_df <- bind_rows(gm_legend_df, wm_legend_df)
  
  p_leg <- ggplot(legend_df) +
    geom_tile(
      aes(x = x, y = y, fill = name),
      width = 0.42,
      height = 0.42,
      color = NA
    ) +
    geom_text(
      aes(x = x + 0.35, y = y, label = label),
      hjust = 0,
      size = 4.5,
      family = "Arial"
    ) +
    scale_fill_manual(values = all_cols, drop = FALSE) +
    scale_y_continuous(
      limits = c(0.4, 2.6),
      breaks = c(1, 2),
      labels = c("WM", "GM")
    ) +
    scale_x_continuous(
      limits = c(0.5, max(length(GM_cols), length(WM_cols)) + 1.2)
    ) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = "Arial") +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 12, family = "Arial"),
      plot.margin = margin(5, 20, 5, 5)
    )
  
  return(p_leg)
}

p_color_legend <- make_color_legend_plot()

ggsave(
  filename = file.path(fig_dir, "GM_WM_color_legend.pdf"),
  plot = p_color_legend,
  width = 8.5,
  height = 1.4,
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_dir, "GM_WM_color_legend.png"),
  plot = p_color_legend,
  width = 8.5,
  height = 1.4,
  dpi = 600
)

message("All figures saved to: ", fig_dir)