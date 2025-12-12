# 5-fold Predicted vs True scatter with per-fold line + CI 
rm(list = ls())

library(R.matlab)
library(ggplot2)
library(visreg)

# ---- data path ----
ProjectFolder <- "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/PNC/code/4th_prediction/nosmooth/age_withHaufe_Schaefer100"
targetStr     <- "PNC_age"
BaseFolder    <- file.path(ProjectFolder, targetStr, "RegressCovariates_RandomCV")
OutDir        <- "/ibmgpfs/cuizaixu_lab/congjing/WM_prediction/Figure/Age_scatterplot/PNC"
if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)


# ---- FC  ----
fc_map <- c("GGFC" = "G-G", "GWFC" = "G-W", "WWFC" = "W-W")

# ---- color ----
fold_cols   <- c("#7F7F7F", 
                 "#000000", 
                 "#9A9A9A",
                 "#595959",  
                 "#8A8A8A")
pt_shape  <- 16
pt_alpha  <- 0.6

for (fc_name in names(fc_map)) {
  message("==== ", fc_name, " (", fc_map[[fc_name]], ") ====")
  
  # ---------- 1) find the median r: Time_* ----------
  corr_vec <- rep(NA_real_, 101)
  for (t in 0:100) {
    f_res <- file.path(BaseFolder, paste0("Time_", t), fc_name, "Res_NFold.mat")
    if (file.exists(f_res)) {
      R   <- readMat(f_res); nms <- names(R)
      if      ("Mean.Corr" %in% nms) corr_vec[t+1] <- as.numeric(R$Mean.Corr)
      else if ("Mean_Corr" %in% nms) corr_vec[t+1] <- as.numeric(R$Mean_Corr)
      else if ("Corr"      %in% nms) corr_vec[t+1] <- as.numeric(R$Corr)
    }
  }
  ord        <- order(corr_vec, decreasing = TRUE, na.last = TRUE)
  MedianTime <- as.integer(ord[51] - 1L)  # 0-based
  cat(sprintf("[%s] Use median Time_%d (median corr=%.3f)\n",
              fc_name, MedianTime, median(corr_vec, na.rm=TRUE)))
  
  # ---------- 2) Read the True/Pred value at 5fold at this 'Time_' ----------
  fold_data_list <- vector("list", 5)  #TrueAge, PredAge, FoldID
  for (k in 0:4) {
    f_sc <- file.path(BaseFolder, paste0("Time_", MedianTime), fc_name, paste0("Fold_", k, "_Score.mat"))
    if (!file.exists(f_sc)) stop("Missing file: ", f_sc)
    
    S   <- readMat(f_sc); nms <- names(S)
    
    # TrueAge（Test_*）
    if      ("Test_Score" %in% nms) TrueAge <- as.numeric(as.vector(S$Test_Score))
    else if ("Test.Score" %in% nms) TrueAge <- as.numeric(as.vector(S$Test.Score))
    else stop("Missing field Test_Score/Test.Score in: ", f_sc)
    
    # PredAge（Predict_*）
    if      ("Predict_Score" %in% nms) PredAge <- as.numeric(as.vector(S$Predict_Score))
    else if ("Predict.Score" %in% nms) PredAge <- as.numeric(as.vector(S$Predict.Score))
    else stop("Missing field Predict_Score/Predict.Score in: ", f_sc)
    
    fold_data_list[[k+1]] <- data.frame(TrueAge = TrueAge,
                                        PredAge = PredAge,
                                        FoldID  = factor(k))  # 0..4
  }
  pooled <- do.call(rbind, fold_data_list)
  
  # ---------- 3) Fitting per fold: lm(PredAge ~ TrueAge) ----------
  curve_list  <- vector("list", 5)  # TrueAge, line, lower, upper, FoldID
  points_list <- vector("list", 5)  # TrueAge, PredAge, FoldID
  r_by_fold   <- rep(NA_real_, 5)
  
  for (k in 1:5) {
    dfk <- fold_data_list[[k]]   
    
    fit_k <- lm(PredAge ~ TrueAge, data = dfk)
    vr_k  <- visreg(fit_k, "TrueAge", type = "conditional", scale = "linear", plot = FALSE)
    
    curve_list[[k]] <- data.frame(
      TrueAge = vr_k$fit[["TrueAge"]],
      line    = vr_k$fit$visregFit,
      lower   = vr_k$fit$visregLwr,
      upper   = vr_k$fit$visregUpr,
      FoldID  = factor(k-1)
    )
    

    points_list[[k]] <- data.frame(
      TrueAge = vr_k$res[["TrueAge"]],
      PredAge = vr_k$res$visregRes,
      FoldID  = factor(k-1)
    )
    
    r_by_fold[k] <- suppressWarnings(cor(dfk$TrueAge, dfk$PredAge, use = "complete.obs"))
  }
  
  mean_r    <- mean(r_by_fold, na.rm = TRUE)
  label_str <- paste0("Mean~italic(r)==", formatC(mean_r, format = "f", digits = 2))
  
  # ---------- 4) range ----------
  # x : 7–25；y: min(PredAge) -- 25
  y_min_src <- pooled$PredAge
  x_min <- 7;  x_max <- 25
  y_min <- floor(min(y_min_src, na.rm = TRUE));  y_max <- 25
  
  # ---------- 5) plot----------
  g <- ggplot()
  for (k in 1:5) {
    col_k <- fold_cols[k]
    g <- g +
      geom_ribbon(data = curve_list[[k]],
                  aes(x = TrueAge, ymin = lower, ymax = upper),
                  fill = col_k, alpha = 0.25) +
      geom_line(data = curve_list[[k]],
                aes(x = TrueAge, y = line),
                colour = col_k, linewidth = 1.4) +
      geom_jitter(data = points_list[[k]],
                  aes(x = TrueAge, y = PredAge),
                  width = 0.2, colour = col_k, alpha = pt_alpha, size = 2, shape = pt_shape)
  }
  
  g <- g +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.grid       = element_blank(),
      plot.title = element_text(face = "bold", size = 16),
      axis.title = element_text(size = 15),
      axis.text  = element_text(size = 13, color = "black")
    ) +
    labs(title = fc_map[[fc_name]],
         x = "Chronological Age (years)",
         y = "Predicted Age (years)") +
    scale_x_continuous(breaks = c(10, 15, 20, 25)) +
    scale_y_continuous(breaks = c(5, 10, 15, 20, 25)) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE) +
    annotate("text",
             x = x_min + 0.58*(x_max - x_min),
             y = y_min + 0.12*(y_max - y_min),
             label = label_str, parse = TRUE,
             hjust = 0, vjust = 0, size = 5, color = "grey15")
  
  # ---------- 6) save ----------
  stub <- file.path(OutDir, paste0(fc_name, "_Time_", MedianTime))
  ggsave(paste0(stub, ".tif"), g, width = 7, height = 6.0, units = "in", dpi = 300, bg = "transparent")
  ggsave(paste0(stub, ".svg"), g, width = 7, height = 6.0, units = "in", bg = "transparent")
  ggsave(paste0(stub, ".pdf"), g, width = 6, height = 6, units = "in", bg = "transparent")
  print(g)
  cat("Saved:", paste0(stub, ".tif/.svg/.pdf"), "\n")
}