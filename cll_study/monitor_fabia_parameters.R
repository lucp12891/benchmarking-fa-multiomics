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
custom_grid <- data.frame(
  alpha = c(0.005,0.01, 0.02, 0.025, 0.03, 0.35, 0.04, 0.05, 0.1),
  spz   = c(0.01, 0.01, 0.02,0.025, 0.03, 0.05, 0.10, 0.20, 0.30, 0.5,0.6, 0.8, 0.9, 1.0),
  spl   = c(0.01, 0.01, 0.02,0.025, 0.03, 0.05, 0.10, 0.20, 0.30, 0.5,0.6, 0.8, 0.9, 1.0)
)


res3 <- run_fabia_grid(
  X = list(mRNA = as.matrix(mRNA), proteomics = as.matrix(proteomics)),
  k = 1,
  scale_features = TRUE,
  grid = custom_grid
)

res3$summary

# Table from res3 (samples proportion = row-wise "any factor" non-zero)
tab_samples <- res3$summary[, c("alpha", "spz", "spl", "prop_nz_samples_any")]
names(tab_samples)[4] <- "proportion_samples"

# Order like your screenshot
tab_samples <- tab_samples[order(tab_samples$alpha), ]

# Optional: round for display (keep more digits if you want)
tab_samples$proportion_samples <- round(tab_samples$proportion_samples, 3)

tab_samples

write.csv(tab_samples, "fabia_samples_table.csv", row.names = FALSE)


# Fallback: compute from res3$fits if needed
tol <- res3$meta$tol

tab_samples <- data.frame(
  alpha = res3$summary$alpha,
  spz   = res3$summary$spz,
  spl   = res3$summary$spl,
  proportion_samples = sapply(res3$fits, function(x) {
    Z <- x$scores
    if (is.null(Z) || nrow(Z) == 0 || ncol(Z) == 0) return(NA_real_)
    mean(apply(abs(Z) > tol, 1, any))
  })
)

tab_samples <- tab_samples[order(tab_samples$alpha), ]
tab_samples$proportion_samples <- round(tab_samples$proportion_samples, 3)

tab_samples

library(ggplot2)

df_long <- rbind(
  data.frame(alpha = tab_samples$alpha, metric = "spz", value = tab_samples$spz),
  data.frame(alpha = tab_samples$alpha, metric = "spl", value = tab_samples$spl),
  data.frame(alpha = tab_samples$alpha, metric = "proportion_samples", value = tab_samples$proportion_samples)
)

df_long$metric <- factor(df_long$metric, levels = c("spz", "spl", "proportion_samples"))

ggplot(df_long, aes(x = alpha, y = value, color = metric, shape = metric, linetype = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  scale_x_continuous(breaks = sort(unique(tab_samples$alpha))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "ALPHA", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---- view-specific feature proportion from res3$fits ----
make_view_feature_table <- function(res3, view, tol = NULL) {
  if (is.null(tol)) tol <- res3$meta$tol
  if (is.null(res3$fits) || length(res3$fits) == 0) stop("res3$fits missing. You must keep fits in run_fabia_grid().")
  
  prop_view <- sapply(res3$fits, function(x) {
    L <- x$loadings
    if (is.null(L) || nrow(L) == 0 || ncol(L) == 0) return(NA_real_)
    
    idx <- startsWith(rownames(L), paste0(view, "::"))
    Lv <- L[idx, , drop = FALSE]
    if (nrow(Lv) == 0) return(NA_real_)
    
    # row-wise "any factor non-zero"
    mean(apply(abs(Lv) > tol, 1, any))
  })
  
  tab <- data.frame(
    alpha = res3$summary$alpha,
    spz   = res3$summary$spz,
    spl   = res3$summary$spl,
    proportion_features = prop_view
  )
  
  tab <- tab[order(tab$alpha), ]
  tab$proportion_features <- round(tab$proportion_features, 3)
  tab
}

tab_mRNA <- make_view_feature_table(res3, view = "mRNA")
tab_prot <- make_view_feature_table(res3, view = "proteomics")

tab_mRNA
tab_prot
