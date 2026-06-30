run_fabia_grid <- function(
    X,
    k = 10,
    scale_features = TRUE,
    grid = NULL,               # data.frame with columns alpha, spz, spl
    alphas = c(0.01, 0.02, 0.05),
    spz_vals = c(0.20),
    spl_vals = c(0.10, 0.20, 0.30),
    cyc = 1000,
    lap = 1.0,
    nL = 2,
    seed = 21,
    tol = 1e-8,
    verbose = TRUE
) {
  .stop_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is required.")
  .as_list_views <- function(X) if (is.list(X)) X else list(omic1 = X)
  .stack_views <- function(lst) {
    do.call(rbind, lapply(names(lst), function(v) {
      M <- lst[[v]]
      rn <- rownames(M); if (is.null(rn)) rn <- paste0("f", seq_len(nrow(M)))
      rownames(M) <- paste0(v, "::", rn)
      M
    }))
  }
  
  # ---- input checks & harmonize ----
  if (!(is.matrix(X) || is.list(X))) stop("X must be a matrix or a list of matrices.")
  V <- .as_list_views(X)
  sample_names <- colnames(V[[1]])
  if (is.null(sample_names)) stop("All matrices must have column (sample) names.")
  for (nm in names(V)) {
    v <- V[[nm]]
    if (is.null(colnames(v))) stop("All matrices must have column (sample) names.")
    if (!setequal(colnames(v), sample_names)) stop("All views must contain the same sample names.")
  }
  V <- lapply(V, function(m) m[, sample_names, drop = FALSE])
  
  # optional scaling (feature-wise)
  if (scale_features) {
    V <- lapply(V, function(m) {
      m[] <- scale(t(m)) |> t()
      m[is.na(m)] <- 0
      m
    })
  }
  
  # build FABIA matrix (features x samples)
  M_fabia <- .stack_views(V)
  p <- min(nrow(M_fabia), k)
  
  # ---- build grid (if not provided) ----
  if (is.null(grid)) {
    grid <- expand.grid(
      alpha = alphas,
      spz   = spz_vals,
      spl   = spl_vals,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    req <- c("alpha", "spz", "spl")
    if (!all(req %in% colnames(grid))) stop("grid must contain columns: ", paste(req, collapse = ", "))
    grid <- as.data.frame(grid)
  }
  
  .stop_pkg("fabia")
  
  # ---- helper: extract L and Z robustly ----
  .extract_LZ <- function(fit, M_fabia) {
    L <- tryCatch(fit@L, error = function(e) NULL)
    Z <- tryCatch(fit@Z, error = function(e) NULL)
    
    # loadings: features x k
    if (!is.null(L)) {
      if (nrow(L) != nrow(M_fabia) && ncol(L) == nrow(M_fabia)) L <- t(L)
      rownames(L) <- rownames(M_fabia)
      if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
    } else {
      L <- matrix(0, nrow = nrow(M_fabia), ncol = 0, dimnames = list(rownames(M_fabia), NULL))
    }
    
    # scores: samples x k
    if (!is.null(Z)) {
      if (nrow(Z) != ncol(M_fabia) && ncol(Z) == ncol(M_fabia)) Z <- t(Z)
      rownames(Z) <- colnames(M_fabia)
      if (ncol(Z) > 0) colnames(Z) <- paste0("F", seq_len(ncol(Z)))
    } else {
      Z <- matrix(0, nrow = ncol(M_fabia), ncol = 0, dimnames = list(colnames(M_fabia), NULL))
    }
    
    list(L = L, Z = Z)
  }
  
  # ---- helper: summarize sparsity ----
  .summarize_sparsity <- function(L, Z, tol) {
    nf <- nrow(L); kL <- ncol(L)
    ns <- nrow(Z); kZ <- ncol(Z)
    
    # entry-wise nonzeros
    L_nz <- if (kL == 0) 0L else sum(abs(L) > tol, na.rm = TRUE)
    Z_nz <- if (kZ == 0) 0L else sum(abs(Z) > tol, na.rm = TRUE)
    
    L_tot <- nf * kL
    Z_tot <- ns * kZ
    
    # row-wise "selected if any factor nonzero"
    feat_selected <- if (kL == 0) 0L else sum(apply(abs(L) > tol, 1, any))
    samp_selected <- if (kZ == 0) 0L else sum(apply(abs(Z) > tol, 1, any))
    
    list(
      n_features = nf,
      n_samples  = ns,
      n_factors_L = kL,
      n_factors_Z = kZ,
      
      nz_entries_features = L_nz,
      nz_entries_samples  = Z_nz,
      prop_nz_entries_features = if (L_tot == 0) NA_real_ else L_nz / L_tot,
      prop_nz_entries_samples  = if (Z_tot == 0) NA_real_ else Z_nz / Z_tot,
      
      nz_features_any = feat_selected,
      nz_samples_any  = samp_selected,
      prop_nz_features_any = if (nf == 0) NA_real_ else feat_selected / nf,
      prop_nz_samples_any  = if (ns == 0) NA_real_ else samp_selected / ns
    )
  }
  
  # ---- run grid ----
  out_tbl <- vector("list", nrow(grid))
  fits <- vector("list", nrow(grid))  # optional: keep fits; set to NULL if you don't want
  
  for (i in seq_len(nrow(grid))) {
    a   <- grid$alpha[i]
    spz <- grid$spz[i]
    spl <- grid$spl[i]
    
    if (verbose) message(sprintf("FABIA run %d/%d: alpha=%.5g spz=%.5g spl=%.5g", i, nrow(grid), a, spz, spl))
    
    set.seed(seed + i)  # different randomness per run while keeping reproducibility
    fit <- fabia::fabia(
      X = M_fabia,
      p = p,
      alpha = a,
      cyc = cyc,
      spz = spz,
      spl = spl,
      lap = lap,
      nL = nL
    )
    
    LZ <- .extract_LZ(fit, M_fabia)
    stats <- .summarize_sparsity(LZ$L, LZ$Z, tol)
    
    out_tbl[[i]] <- data.frame(
      run = i,
      alpha = a,
      spz = spz,
      spl = spl,
      stats,
      stringsAsFactors = FALSE
    )
    
    fits[[i]] <- list(fit = fit, loadings = LZ$L, scores = LZ$Z)
  }
  
  summary_table <- do.call(rbind, out_tbl)
  
  # rank helper: you can change objective
  summary_table$rank_sparse_entries <- rank(-summary_table$prop_nz_entries_features - summary_table$prop_nz_entries_samples, ties.method = "min")
  
  list(
    summary = summary_table,
    fits = fits,
    meta = list(
      k = k,
      p = p,
      scale_features = scale_features,
      tol = tol,
      seed_base = seed,
      grid = grid
    )
  )
}

set.seed(1234)

# Example 1: keep spz constant, vary alpha + spl
# res1 <- run_fabia_grid(
#   X = list(mRNA = as.matrix(mRNA), proteomics = as.matrix(proteomics)),
#   k = 1,
#   scale_features = TRUE,
#   spz_vals = 0.20,
#   alphas = c(0.005, 0.01, 0.02, 0.05),
#   spl_vals = c(0.05, 0.10, 0.20, 0.30)
# )
# # res1$summary[order(res1$summary$rank_sparse_entries), ]
# alphas_grid <- unique(round(c(
#   0.005,
#   0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0
# ), 6))
# 
# spl_grid <- unique(round(seq(0, 1, by = 0.1), 2))   # 11 values
# spz_grid <- unique(round(seq(0, 1, by = 0.1), 2))   # 11 values
# 
# res1 <- run_fabia_grid(
#   X = list(mRNA = as.matrix(mRNA), proteomics = as.matrix(proteomics)),
#   k = 1,
#   scale_features = TRUE,
#   spz_vals = 0.20,
#   alphas = alphas_grid,
#   spl_vals = spl_grid
# )
# 
# 
# # Example 2: keep alpha constant, vary spz + spl
# res2 <- run_fabia_grid(
#   X = list(mRNA = as.matrix(mRNA), proteomics = as.matrix(proteomics)),
#   k = 1,
#   scale_features = TRUE,
#   alphas = 0.02,
#   spz_vals = c(0.05, 0.10, 0.20, 0.30),
#   spl_vals = c(0.05, 0.10, 0.20, 0.30)
# )
# 
# res2$summary

# Example 3: custom grid (full control)
# custom_grid <- data.frame(
#   alpha = c(0.005,0.01, 0.02, 0.025, 0.03, 0.35, 0.04, 0.05, 0.1),
#   spz   = c(0.01, 0.01, 0.02,0.025, 0.03, 0.05, 0.10, 0.20, 0.30, 0.5,0.6, 0.8, 0.9, 1.0),
#   spl   = c(0.01, 0.01, 0.02,0.025, 0.03, 0.05, 0.10, 0.20, 0.30, 0.5,0.6, 0.8, 0.9, 1.0)
# )
# 
# 
# res3 <- run_fabia_grid(
#   X = list(prot = as.matrix(prot), mRNA = as.matrix(mRNAation)),
#   k = 1,
#   scale_features = TRUE,
#   grid = custom_grid
# )
# 
# # res_rad_alpha_spl_spz
# custom_grid <- data.frame(
#   alpha = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1),
#   spz   = c(0.01, 0.05, 0.10, 0.20, 0.30, 0.5),
#   spl   = c(0.01, 0.05, 0.10, 0.20, 0.30, 0.5)
# )
run_fabia_grid <- function(
    X,
    k = 10,
    scale_features = TRUE,
    grid = NULL,               # data.frame with columns alpha, spz, spl
    alphas = c(0.01, 0.02, 0.05),
    spz_vals = c(0.20),
    spl_vals = c(0.10, 0.20, 0.30),
    cyc = 1000,
    lap = 1.0,
    nL = 2,
    seed = 21,
    tol = 1e-8,
    verbose = TRUE
) {
  .stop_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is required.")
  .as_list_views <- function(X) if (is.list(X)) X else list(omic1 = X)
  .stack_views <- function(lst) {
    do.call(rbind, lapply(names(lst), function(v) {
      M <- lst[[v]]
      rn <- rownames(M); if (is.null(rn)) rn <- paste0("f", seq_len(nrow(M)))
      rownames(M) <- paste0(v, "::", rn)
      M
    }))
  }
  
  # ---- input checks & harmonize ----
  if (!(is.matrix(X) || is.list(X))) stop("X must be a matrix or a list of matrices.")
  V <- .as_list_views(X)
  sample_names <- colnames(V[[1]])
  if (is.null(sample_names)) stop("All matrices must have column (sample) names.")
  for (nm in names(V)) {
    v <- V[[nm]]
    if (is.null(colnames(v))) stop("All matrices must have column (sample) names.")
    if (!setequal(colnames(v), sample_names)) stop("All views must contain the same sample names.")
  }
  V <- lapply(V, function(m) m[, sample_names, drop = FALSE])
  
  # optional scaling (feature-wise)
  if (scale_features) {
    V <- lapply(V, function(m) {
      m[] <- scale(t(m)) |> t()
      m[is.na(m)] <- 0
      m
    })
  }
  
  # build FABIA matrix (features x samples)
  M_fabia <- .stack_views(V)
  p <- min(nrow(M_fabia), k)
  
  # ---- build grid (if not provided) ----
  if (is.null(grid)) {
    grid <- expand.grid(
      alpha = alphas,
      spz   = spz_vals,
      spl   = spl_vals,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    req <- c("alpha", "spz", "spl")
    if (!all(req %in% colnames(grid))) stop("grid must contain columns: ", paste(req, collapse = ", "))
    grid <- as.data.frame(grid)
  }
  
  .stop_pkg("fabia")
  
  # ---- helper: extract L and Z robustly ----
  .extract_LZ <- function(fit, M_fabia) {
    L <- tryCatch(fit@L, error = function(e) NULL)
    Z <- tryCatch(fit@Z, error = function(e) NULL)
    
    # loadings: features x k
    if (!is.null(L)) {
      if (nrow(L) != nrow(M_fabia) && ncol(L) == nrow(M_fabia)) L <- t(L)
      rownames(L) <- rownames(M_fabia)
      if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
    } else {
      L <- matrix(0, nrow = nrow(M_fabia), ncol = 0, dimnames = list(rownames(M_fabia), NULL))
    }
    
    # scores: samples x k
    if (!is.null(Z)) {
      if (nrow(Z) != ncol(M_fabia) && ncol(Z) == ncol(M_fabia)) Z <- t(Z)
      rownames(Z) <- colnames(M_fabia)
      if (ncol(Z) > 0) colnames(Z) <- paste0("F", seq_len(ncol(Z)))
    } else {
      Z <- matrix(0, nrow = ncol(M_fabia), ncol = 0, dimnames = list(colnames(M_fabia), NULL))
    }
    
    list(L = L, Z = Z)
  }
  
  # ---- helper: summarize sparsity ----
  # .summarize_sparsity <- function(L, Z, tol) {
  #   nf <- nrow(L); kL <- ncol(L)
  #   ns <- nrow(Z); kZ <- ncol(Z)
  #   
  #   # entry-wise nonzeros
  #   L_nz <- if (kL == 0) 0L else sum(abs(L) > tol, na.rm = TRUE)
  #   Z_nz <- if (kZ == 0) 0L else sum(abs(Z) > tol, na.rm = TRUE)
  #   
  #   L_tot <- nf * kL
  #   Z_tot <- ns * kZ
  #   
  #   # row-wise "selected if any factor nonzero"
  #   feat_selected <- if (kL == 0) 0L else sum(apply(abs(L) > tol, 1, any))
  #   samp_selected <- if (kZ == 0) 0L else sum(apply(abs(Z) > tol, 1, any))
  #   
  #   list(
  #     n_features = nf,
  #     n_samples  = ns,
  #     n_factors_L = kL,
  #     n_factors_Z = kZ,
  #     
  #     nz_entries_features = L_nz,
  #     nz_entries_samples  = Z_nz,
  #     prop_nz_entries_features = if (L_tot == 0) NA_real_ else L_nz / L_tot,
  #     prop_nz_entries_samples  = if (Z_tot == 0) NA_real_ else Z_nz / Z_tot,
  #     
  #     nz_features_any = feat_selected,
  #     nz_samples_any  = samp_selected,
  #     prop_nz_features_any = if (nf == 0) NA_real_ else feat_selected / nf,
  #     prop_nz_samples_any  = if (ns == 0) NA_real_ else samp_selected / ns
  #   )
  # }
  .summarize_sparsity <- function(L, Z, tol, view_names = NULL) {
    nf <- nrow(L); kL <- ncol(L)
    ns <- nrow(Z); kZ <- ncol(Z)
    
    # entry-wise nonzeros
    L_nz <- if (kL == 0) 0L else sum(abs(L) > tol, na.rm = TRUE)
    Z_nz <- if (kZ == 0) 0L else sum(abs(Z) > tol, na.rm = TRUE)
    
    L_tot <- nf * kL
    Z_tot <- ns * kZ
    
    # row-wise "selected if any factor nonzero"
    feat_selected <- if (kL == 0) 0L else sum(apply(abs(L) > tol, 1, any))
    samp_selected <- if (kZ == 0) 0L else sum(apply(abs(Z) > tol, 1, any))
    
    out <- list(
      n_features = nf,
      n_samples  = ns,
      n_factors_L = kL,
      n_factors_Z = kZ,
      
      nz_entries_features = L_nz,
      nz_entries_samples  = Z_nz,
      prop_nz_entries_features = if (L_tot == 0) NA_real_ else L_nz / L_tot,
      prop_nz_entries_samples  = if (Z_tot == 0) NA_real_ else Z_nz / Z_tot,
      
      nz_features_any = feat_selected,
      nz_samples_any  = samp_selected,
      prop_nz_features_any = if (nf == 0) NA_real_ else feat_selected / nf,
      prop_nz_samples_any  = if (ns == 0) NA_real_ else samp_selected / ns
    )
    
    # ---- per-view feature selection (requires rownames with "view::") ----
    rn <- rownames(L)
    if (!is.null(rn) && kL > 0) {
      view_of_row <- sub("::.*$", "", rn)
      
      # If caller provides view_names, keep a stable set/order; otherwise infer from L.
      if (is.null(view_names)) view_names <- sort(unique(view_of_row))
      
      for (v in view_names) {
        idx <- view_of_row == v
        n_v <- sum(idx)
        if (n_v == 0) {
          out[[paste0("nz_features_any_", v)]] <- NA_integer_
          out[[paste0("prop_nz_features_any_", v)]] <- NA_real_
        } else {
          selected_v <- sum(apply(abs(L[idx, , drop = FALSE]) > tol, 1, any))
          out[[paste0("nz_features_any__", v)]] <- selected_v
          out[[paste0("prop_nz_features_any_", v)]] <- selected_v / n_v
        }
      }
    }
    
    out
  }
  
  # ---- run grid ----
  out_tbl <- vector("list", nrow(grid))
  fits <- vector("list", nrow(grid))  # optional: keep fits; set to NULL if you don't want
  
  for (i in seq_len(nrow(grid))) {
    a   <- grid$alpha[i]
    spz <- grid$spz[i]
    spl <- grid$spl[i]
    
    if (verbose) message(sprintf("FABIA run %d/%d: alpha=%.5g spz=%.5g spl=%.5g", i, nrow(grid), a, spz, spl))
    
    set.seed(seed + i)  # different randomness per run while keeping reproducibility
    fit <- fabia::fabia(
      X = M_fabia,
      p = p,
      alpha = a,
      cyc = cyc,
      spz = spz,
      spl = spl,
      lap = lap,
      nL = nL
    )
    
    LZ <- .extract_LZ(fit, M_fabia)
    stats <- .summarize_sparsity(LZ$L, LZ$Z, tol, view_names = names(V)) # .summarize_sparsity(LZ$L, LZ$Z, tol)
    
    out_tbl[[i]] <- data.frame(
      run = i,
      alpha = a,
      spz = spz,
      spl = spl,
      stats,
      stringsAsFactors = FALSE
    )
    
    fits[[i]] <- list(fit = fit, loadings = LZ$L, scores = LZ$Z)
  }
  
  summary_table <- do.call(rbind, out_tbl)
  
  # rank helper: you can change objective
  summary_table$rank_sparse_entries <- rank(-summary_table$prop_nz_entries_features - summary_table$prop_nz_entries_samples, ties.method = "min")
  
  list(
    summary = summary_table,
    fits = fits,
    meta = list(
      k = k,
      p = p,
      scale_features = scale_features,
      tol = tol,
      seed_base = seed,
      grid = grid
    )
  )
}

#writexl::write_xlsx(mRNA, path = "C:/Users/bosangir/OneDrive - Studiecentrum voor Kernenergie/Desktop/Ausan - DELL/Output - Paper II/R project/paper_II/mRNA.xlsx")
#writexl::write_xlsx(proteomics, path = "C:/Users/bosangir/OneDrive - Studiecentrum voor Kernenergie/Desktop/Ausan - DELL/Output - Paper II/R project/paper_II/proteomics.xlsx")

library(readxl)

#mRNA <- read_xlsx("C:/Users/bosangir/OneDrive - Studiecentrum voor Kernenergie/Desktop/Ausan - DELL/Output - Paper II/R project/paper_II//mRNA.xlsx")
#proteomics <- read_xlsx("C:/Users/bosangir/OneDrive - Studiecentrum voor Kernenergie/Desktop/Ausan - DELL/Output - Paper II/R project/paper_II//proteomics.xlsx")

mRNA <- read.csv("D17_mRNA_Hippocampus_normalized_tmm.csv", header = TRUE, check.names = FALSE, row.names = 1)
proteomics <- read.csv("D17_Protein_Hippocampus_normalized_quant.csv", header = TRUE, check.names = FALSE, row.names = 1)

set.seed(123)
library(fabia)
res_rad_alpha_spl_spzX <- run_fabia_grid(
  X = list(mRNA = as.matrix(mRNA), proteomics = as.matrix(proteomics)),
  k = 2, #2
  scale_features = TRUE,
  seed = 123,
  spz_vals = c(0.01, 0.05, 0.10, 0.20, 0.30, 0.5),
  alphas = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1),
  spl_vals = c(0.01, 0.05, 0.10, 0.20, 0.30, 0.5)
)


res_rad_alpha_spl_spzX$summary

tab_mRNA <- res_rad_alpha_spl_spzX$summary[, c("alpha","spz","spl","prop_nz_features_any_mRNA")]
names(tab_mRNA)[4] <- "proportion_features"

tab_prot <- res_rad_alpha_spl_spzX$summary[, c("alpha","spz","spl","prop_nz_features_any_proteomics")]
names(tab_prot)[4] <- "proportion_features"

res_rad_alpha_spl_spzX$summary

library(ggplot2)

df = res_rad_alpha_spl_spzX
df$alpha <- factor(df$alpha)
df$spz   <- factor(df$spz, levels = sort(unique(df$spz)))
df$spl   <- factor(df$spl, levels = sort(unique(df$spl)))

res_rad_alpha_spl_spzX$summary

# Table from res_rad_alpha_spl_spz (samples proportion = row-wise "any factor" non-zero)
tab_samples_rad <- df$summary[, c("alpha", "spz", "spl", "prop_nz_samples_any")]
names(tab_samples_rad)[4] <- "proportion_samples"

# Order like your screenshot
tab_samples_rad <- tab_samples_rad[order(tab_samples_rad$alpha), ]

# Optional: round for display (keep more digits if you want)
tab_samples_rad$proportion_samples <- round(tab_samples_rad$proportion_samples, 3)

tab_samples_rad

write.csv(tab_samples_rad, "fabia_samples_table_rad.csv", row.names = FALSE)
write.csv(tab_mRNA, "tab_mRNA_rad.csv", row.names = FALSE)


library(dplyr)
# library(ggplot2)
# ggplot(tab_prot,
#        aes(x = alpha,
#            y = proportion_samples,
#            color = factor(spz),
#            group = spz)) +
#   geom_line(linewidth = 0.7) +
#   geom_point(size = 1.2) +
#   facet_wrap(~ spl, ncol = 3) +
#   scale_y_continuous(limits = c(0, 1)) +
#   scale_color_brewer(palette = "Blues") +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.title = element_blank(),
#     panel.grid.minor = element_blank()
#   )
library(ggplot2)

best_row_samples_rad <- subset(
  tab_samples_rad,
  spz == 0.1 & spl == 0.5 & alpha == 0.03
)
ggplot(tab_samples_rad,
       aes(x = alpha,
           y = proportion_samples,
           group = interaction(spz, spl),
           color = factor(spz),
           linetype = factor(spl))) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_point(size = 1.6) +
  
  # ---- highlight best combination ----
geom_line(
  data = subset(tab_samples_rad, spz == 0.1 & spl == 0.5),
  aes(x = alpha, y = proportion_samples, group = interaction(spz, spl)),
  inherit.aes = FALSE,
  color = "blue",
  linewidth = 1.6
) +
  geom_point(
    data = best_row_samples_rad,
    aes(x = alpha, y = proportion_samples),
    inherit.aes = FALSE,
    color = "blue",
    size = 3
  ) +
  geom_vline(xintercept = 0.03, linewidth = 0.7, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0, 1),
                     name = "Proportion of non-zero samples") +
  scale_x_continuous(breaks = sort(unique(tab_samples_rad$alpha)),
                     name = "alpha") +
  scale_color_brewer(palette = "Blues", name = "spz") +
  scale_linetype_discrete(name = "spl") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    axis.text.x = element_text(size = 14, angle =30, hjust = 1, vjust = 1),
    legend.position = "right",
    #axis.title.y = element_text(margin = margin(r = 10)),
    #axis.title.x = element_text(margin = margin(t = 10)),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  )


library(dplyr)

# ---- set this to your actual feature-count column name ----
feature_col <- "proportion_features"   # change if needed

# ---- view-specific feature proportion from res_rad_alpha_spl_spz$fits ----
make_view_feature_table <- function(res_rad_alpha_spl_spz, view, tol = NULL) {
  if (is.null(tol)) tol <- res_rad_alpha_spl_spz$meta$tol
  if (is.null(res_rad_alpha_spl_spz$fits) || length(res_rad_alpha_spl_spz$fits) == 0) stop("res_rad_alpha_spl_spz$fits missing. You must keep fits in run_fabia_grid().")
  
  prop_view <- sapply(res_rad_alpha_spl_spz$fits, function(x) {
    L <- x$loadings
    if (is.null(L) || nrow(L) == 0 || ncol(L) == 0) return(NA_real_)
    
    idx <- startsWith(rownames(L), paste0(view, "::"))
    Lv <- L[idx, , drop = FALSE]
    if (nrow(Lv) == 0) return(NA_real_)
    
    # row-wise "any factor non-zero"
    mean(apply(abs(Lv) > tol, 1, any))
  })
  
  tab <- data.frame(
    alpha = res_rad_alpha_spl_spz$summary$alpha,
    spz   = res_rad_alpha_spl_spz$summary$spz,
    spl   = res_rad_alpha_spl_spz$summary$spl,
    proportion_features = prop_view
  )
  
  tab <- tab[order(tab$alpha), ]
  tab$proportion_features <- round(tab$proportion_features, 3)
  tab
}

tab_mRNA <- make_view_feature_table(res_rad_alpha_spl_spzX, view = "mRNA")
tab_prot <- make_view_feature_table(res_rad_alpha_spl_spzX, view = "proteomics")

tab_mRNA
tab_prot

best_row_mRNA <- tab_mRNA %>%
  filter(.data[[feature_col]] < 0.49, alpha %in% c(0.02, 0.03)) %>%
  arrange(desc(proportion_features), alpha, .data[[feature_col]]) %>%
  slice(3)
best_row_mRNA <- subset(
  tab_mRNA,
  spz == 0.1 & spl == 0.5 & alpha == 0.03
)
best_row_mRNA <- subset(
  tab_mRNA,
  spz == 0.2 & spl == 0.5 & alpha == 0.03 # 0.348
)
best_row_mRNA <- subset(
  tab_mRNA,
  spz == 0.5 & spl == 0.05 & alpha == 0.03 # 0.348
)

best_row_mRNA
best_alpha_mRNA <- best_row_mRNA$alpha[1]
best_spz_mRNA   <- best_row_mRNA$spz[1]
best_spl_mRNA   <- best_row_mRNA$spl[1]

best_row_prot <- tab_prot %>%
  filter(.data[[feature_col]] < 0.25, alpha %in% c(0.02, 0.03)) %>%
  arrange(desc(proportion_features), alpha, .data[[feature_col]]) %>%
  slice(3)
best_row_prot <- subset(
  tab_prot,
  spz == 0.1 & spl == 0.5 & alpha == 0.03
)
best_row_prot <- subset(
  tab_prot,
  spz == 0.2 & spl == 0.5 & alpha == 0.03
)

best_row_prot <- subset(
  tab_prot,
  spz == 0.5 & spl == 0.05 & alpha == 0.03
)

best_row_prot
best_alpha_prot <- best_row_prot$alpha[1]
best_spz_prot   <- best_row_prot$spz[1]
best_spl_prot   <- best_row_prot$spl[1]

# data for the highlighted curve only
df_best_mRNA <- tab_mRNA %>%
  filter(spz == best_spz_mRNA, spl == best_spl_mRNA)

tab_mRNA$alpha <- as.numeric(as.character(tab_mRNA$alpha))
df_best_mRNA$alpha <- as.numeric(as.character(df_best_mRNA$alpha))
best_alpha_mRNA <- as.numeric(best_alpha_mRNA)

ggplot(tab_mRNA,
       aes(x = alpha,
           y = proportion_features,
           group = interaction(spz, spl),
           color = factor(spz),
           linetype = factor(spl))) +
  geom_line(linewidth = 0.8, alpha = 3) +
  geom_point(size = 1.6) +
  geom_line(data = df_best_mRNA,
            aes(x = alpha, y = proportion_features, group = interaction(spz, spl)),
            inherit.aes = FALSE,
            linewidth = 1.6,
            color = "red") +
  geom_point(data = df_best_mRNA,
             aes(x = alpha, y = proportion_features),
             inherit.aes = FALSE,
             size = 2.2,
             color = "red") +
  geom_vline(xintercept = 0.03, linewidth = 0.7, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.154,
             linewidth = 1,
             linetype = "dashed",
             color = "red") +
  scale_y_continuous(limits = c(0, 1),
                     name = "Proportion of non-zero features") +
  scale_x_continuous(breaks = sort(unique(tab_mRNA$alpha)),
                     name = "alpha") +
  scale_color_brewer(palette = "Blues", name = "spz") +
  scale_linetype_discrete(name = "spl") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    axis.text.x = element_text(size = 14, angle =30, hjust = 1, vjust = 1),
    legend.position = "right",
    #axis.title.y = element_text(margin = margin(r = 10)),
    #axis.title.x = element_text(margin = margin(t = 10)),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  )

# data for the highlighted curve only
df_best_prot <- tab_prot %>%
  filter(spz == best_spz_prot, spl == best_spl_prot)

tab_prot$alpha <- as.numeric(as.character(tab_prot$alpha))
df_best_prot$alpha <- as.numeric(as.character(df_best_prot$alpha))
best_alpha_prot <- as.numeric(best_alpha_prot)

ggplot(tab_prot,
       aes(x = alpha,
           y = proportion_features,
           group = interaction(spz, spl),
           color = factor(spz),
           linetype = factor(spl))) +
  geom_line(linewidth = 0.8, alpha = 3) +
  geom_point(size = 1.6) +
  geom_line(data = df_best_prot,
            aes(x = alpha, y = proportion_features, group = interaction(spz, spl)),
            inherit.aes = FALSE,
            linewidth = 1.6,
            color = "red") +
  geom_point(data = df_best_prot,
             aes(x = alpha, y = proportion_features),
             inherit.aes = FALSE,
             size = 2.2,
             color = "red") +
  geom_vline(xintercept = 0.03, linewidth = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.002,
             linewidth = 1,
             linetype = "dashed",
             color = "red") +
  scale_y_continuous(limits = c(0, 1),
                     name = "Proportion of non-zero features") +
  scale_x_continuous(breaks = sort(unique(tab_prot$alpha)),
                     name = "alpha") +
  scale_color_brewer(palette = "Blues", name = "spz") +
  scale_linetype_discrete(name = "spl") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    axis.text.x = element_text(size = 14, angle =30, hjust = 1, vjust = 1),
    legend.position = "right",
    #axis.title.y = element_text(margin = margin(r = 10)),
    #axis.title.x = element_text(margin = margin(t = 10)),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  )

# ggplot(tab_samples,
#        aes(x = alpha,
#            y = proportion_samples,
#            group = interaction(spz, spl),
#            color = factor(spz),
#            linetype = factor(spl))) +
#   geom_line(linewidth = 0.8, alpha = 0.9) +
#   geom_vline(xintercept = 0.02, linewidth = 0.7, linetype = "dashed") +
#   geom_point(size = 1.6) +
#   scale_y_continuous(limits = c(0, 1)) +
#   scale_color_brewer(palette = "Blues", name = "SPZ") +
#   scale_linetype_discrete(name = "SPL") +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.title = element_blank(),
#     panel.grid.minor = element_blank()
#   )

# ggplot(tab_prot,
#        aes(x = alpha,
#            y = proportion_features,
#            group = interaction(spz, spl),
#            color = factor(spz),
#            linetype = factor(spl))) +
#   geom_line(linewidth = 0.8, alpha = 0.9) +
#   geom_point(size = 1.6) +
#   scale_y_continuous(limits = c(0, 1)) +
#   scale_color_brewer(palette = "Blues", name = "SPZ") +
#   scale_linetype_discrete(name = "SPL") +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.title = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# ggplot(tab_mRNA,
#        aes(x = alpha,
#            y = proportion_features,
#            group = interaction(spz, spl),
#            color = factor(spz),
#            linetype = factor(spl))) +
#   geom_line(linewidth = 0.8, alpha = 0.9) +
#   geom_point(size = 1.6) +
#   scale_y_continuous(limits = c(0, 1)) +
#   scale_color_brewer(palette = "Blues", name = "SPZ") +
#   scale_linetype_discrete(name = "SPL") +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.title = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# ggplot(tab_mRNA,
#        aes(x = alpha,
#            y = proportion_features,
#            group = interaction(spz, spl),
#            color = interaction(spz, spl))) +
#   geom_line(linewidth = 0.8, alpha = 0.9) +
#   geom_point(size = 1.6) +
#   scale_y_continuous(limits = c(0, 1)) +
#   scale_color_viridis_d(name = "SPZ–SPL", option = "C") +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.title = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# ggplot(tab_mRNA,
#        aes(x = alpha,
#            y = proportion_features,
#            group = interaction(spz, spl),
#            color = interaction(spz, spl))) +
#   geom_line(linewidth = 0.8, alpha = 0.9) +
#   geom_point(size = 1.6) +
#   scale_y_continuous(limits = c(0, 1)) +
#   scale_color_viridis_d(name = "SPZ–SPL", option = "C") +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.title = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# tab_mRNA$alpha <- factor(
#   tab_mRNA$alpha,
#   levels = sort(unique(tab_mRNA$alpha))
# )
# ggplot(tab_mRNA,
#        aes(x = alpha,
#            y = proportion_features,
#            group = interaction(spz, spl),
#            color = interaction(spz, spl))) +
#   geom_line(linewidth = 0.8, alpha = 0.9) +
#   geom_point(size = 1.6) +
#   scale_y_continuous(limits = c(0, 1)) +
#   scale_x_discrete(
#     breaks = levels(tab_mRNA$alpha)   # <-- THIS is the key line
#   ) +
#   scale_color_viridis_d(name = "SPZ–SPL", option = "C") +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.title = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# ggplot(tab_mRNA,
#        aes(x = alpha,
#            y = proportion_features,
#            group = interaction(spz, spl),
#            color = interaction(spz, spl))) +
#   geom_line(linewidth = 0.8, alpha = 0.9) +
#   geom_point(size = 1.6) +
#   scale_y_continuous(
#     limits = c(0, 1),
#     name = "Proportion of non-zero features"
#   ) +
#   scale_x_discrete(
#     name = "FABIA (alpha)",
#     breaks = levels(tab_mRNA$alpha),
#     guide = guide_axis(angle = 0)
#   ) +
#   scale_color_viridis_d(
#     name = "SPZ–SPL",
#     option = "C"
#   ) +
#   guides(
#     color = guide_legend(
#       nrow = 3,        # increase to spread wider
#       byrow = TRUE
#     )
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     panel.grid.minor = element_blank(),
#     axis.title.y = element_text(margin = margin(r = 10)),
#     axis.title.x = element_text(margin = margin(t = 10)),
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     legend.box.just = "center",
#     legend.spacing.x = unit(0.8, "cm"),
#     legend.title = element_text(size = 11),
#     legend.text = element_text(size = 7)
#   )
