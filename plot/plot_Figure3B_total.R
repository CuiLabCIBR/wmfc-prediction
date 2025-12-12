# 5-fold Predicted vs True scatter 
rm(list = ls())

library(R.matlab)
library(ggplot2)
library(visreg)
library(grid)   # unit()

fold_cols <- c("#7F7F7F", "#9A9A9A", "#8A8A8A", "#000000", "#595959")
pt_shape  <- 16
pt_alpha  <- 0.30
pt_size   <- 2
jit_w     <- 0.12
rib_alpha <- 0.22
line_w    <- 1.4

# ========= datapath =========
ProjectFolder <- "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/ABCD/code/5th_prediction/PLS_prediction/networkFC/cognition_nosmooth_motion2runFD"
targetStr     <- "nihtbx_totalcomp"
BaseFolder    <- file.path(ProjectFolder, targetStr, "RegressCovariates_RandomCV")
OutDir        <- "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/Figure/cog_scatterplot"
if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)

# FC 
fc_map <- c("GGFC" = "G-G", "GWFC" = "G-W", "WWFC" = "W-W")

for (fc_name in names(fc_map)) {
  message("==== ", fc_name, " (", fc_map[[fc_name]], ") ====")
  
  # ---- 1) find the median r: Time_*  ----
  corr_vec <- rep(NA_real_, 101)
  for (t in 0:100) {
    f_res <- file.path(BaseFolder, paste0("Time_", t), fc_name, "Res_NFold.mat")
    if (file.exists(f_res)) {
      R <- readMat(f_res); nms <- names(R)
      if      ("Mean.Corr" %in% nms) corr_vec[t+1] <- as.numeric(R$Mean.Corr)
      else if ("Mean_Corr" %in% nms) corr_vec[t+1] <- as.numeric(R$Mean_Corr)
      else if ("Corr"      %in% nms) corr_vec[t+1] <- as.numeric(R$Corr)
    }
  }
  ord        <- order(corr_vec, decreasing = TRUE, na.last = TRUE)
  MedianTime <- as.integer(ord[51] - 1L)
  cat(sprintf("[%s] Use median Time_%d (median corr=%.3f)\n",
              fc_name, MedianTime, median(corr_vec, na.rm=TRUE)))
  
  # ---- 2) Read the True/Pred value at 5fold at this 'Time_' ----
  fold_data_list <- vector("list", 5)
  for (k in 0:4) {
    f_sc <- file.path(BaseFolder, paste0("Time_", MedianTime), fc_name, paste0("Fold_", k, "_Score.mat"))
    if (!file.exists(f_sc)) stop("Missing file: ", f_sc)
    S <- readMat(f_sc); nms <- names(S)
    
    if      ("Test_Score" %in% nms) TrueScore <- as.numeric(as.vector(S$Test_Score))
    else if ("Test.Score" %in% nms) TrueScore <- as.numeric(as.vector(S$Test.Score))
    else stop("Missing field Test_Score/Test.Score in: ", f_sc)
    
    if      ("Predict_Score" %in% nms) PredScore <- as.numeric(as.vector(S$Predict_Score))
    else if ("Predict.Score" %in% nms) PredScore <- as.numeric(as.vector(S$Predict.Score))
    else stop("Missing field Predict_Score/Predict.Score in: ", f_sc)
    
    fold_data_list[[k+1]] <- data.frame(TrueScore = TrueScore,
                                        PredScore = PredScore,
                                        FoldID    = factor(k))
  }
  
  # ---- 3) Fitting per fold: lm(PredAge ~ TrueAge)----
  curve_list  <- vector("list", 5)
  points_list <- vector("list", 5)
  r_by_fold   <- rep(NA_real_, 5)
  
  for (k in 1:5) {
    dfk  <- fold_data_list[[k]]
    fitk <- lm(PredScore ~ TrueScore, data = dfk)
    vr   <- visreg(fitk, "TrueScore", type = "conditional", scale = "linear", plot = FALSE)
    
    curve_list[[k]] <- data.frame(
      TrueScore = vr$fit[["TrueScore"]],
      line      = vr$fit$visregFit,
      lower     = vr$fit$visregLwr,
      upper     = vr$fit$visregUpr,
      FoldID    = factor(k-1)
    )
    points_list[[k]] <- data.frame(
      TrueScore = vr$res[["TrueScore"]],
      PredScore = vr$res$visregRes,
      FoldID    = factor(k-1)
    )
    
    r_by_fold[k] <- suppressWarnings(cor(dfk$TrueScore, dfk$PredScore, use = "complete.obs"))
  }
  
  # ---- 4) y-axis range： PredScore(min-1 / max+1)----
  pooled_pts <- do.call(rbind, points_list)
  yr_raw <- range(pooled_pts$PredScore, na.rm = TRUE)
  
  y_min <- floor(yr_raw[1]) - 1
  y_max <- ceiling(yr_raw[2]) + 1
  
  b_start <- ceiling(y_min / 5) * 5
  b_end   <- floor(y_max / 5) * 5
  y_breaks <- if (b_start <= b_end) seq(b_start, b_end, by = 5) else NULL  
  
  annot_x <- 45 + 0.05 * (120 - 45)
  annot_y <- y_max - 0.08 * (y_max - y_min)
  
  
  # ---- 5) plot ----
  x_min <- 45; x_max <- 120
  
  g <- ggplot()
  for (k in 1:5) {
    col_k <- fold_cols[k]
    g <- g +
      geom_ribbon(data = curve_list[[k]],
                  aes(x = TrueScore, ymin = lower, ymax = upper),
                  fill = col_k, alpha = rib_alpha) +
      geom_line(data = curve_list[[k]],
                aes(x = TrueScore, y = line),
                colour = col_k, linewidth = line_w) +
      geom_jitter(data = points_list[[k]],
                  aes(x = TrueScore, y = PredScore),
                  width = jit_w, colour = col_k, alpha = pt_alpha, size = pt_size, shape = pt_shape)
  }
  
  mean_r    <- mean(r_by_fold, na.rm = TRUE)
  label_str <- paste0("Mean~italic(r)==", formatC(mean_r, format = "f", digits = 2))
  
  g <- g +
    theme_classic(base_size = 13) +
    theme(
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text  = element_text(size = 16, color = "black"),
      axis.ticks.length = unit(2.5, "pt"),
      plot.margin = margin(4, 4, 4, 4),
      legend.position = "none"
    ) +
    labs(title = fc_map[[fc_name]],
         x = "Actual Total Cognition Score",
         y = "Predicted Total Cognition Score") +
    scale_x_continuous(breaks = seq(50, 110, by = 10),
                       minor_breaks = seq(45, 120, by = 5)) +
    scale_y_continuous(breaks = y_breaks) +
    coord_cartesian(xlim = c(45, 120), ylim = c(y_min, y_max), expand = FALSE) +
    annotate("text",
             x = annot_x, y = annot_y,
             label = label_str, parse = TRUE,
             hjust = 0, vjust = 1, size = 6, color = "grey20")
  
  # ---- 6) save ----
  stub <- file.path(OutDir, paste0(fc_name, "_Time_", MedianTime, "_raw"))
  ggsave(paste0(stub, ".tif"), g, width = 6, height = 6, units = "in", dpi = 300, bg = "transparent")
  ggsave(paste0(stub, ".svg"), g, width = 6, height = 6, units = "in", bg = "transparent")
  ggsave(paste0(stub, ".pdf"), g, width = 6, height = 6, units = "in", bg = "transparent")
  print(g)
  cat("Saved:", paste0(stub, ".tif/.svg"), "\n")
}