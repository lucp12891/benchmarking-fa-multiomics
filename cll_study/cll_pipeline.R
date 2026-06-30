#===========================================
# Libraries
#===========================================
library(fabia); library(MOFA2); library(MOFAdata);library(MOFAdata);library(data.table); 
library(gridExtra); library(ggplot2); library(tidyverse); library(dplyr); library(GFA);
library(dplyr); library(reshape2); library(patchwork); library(grid)

# ---------------------------
# Factorization mega-pipeline
# ---------------------------
# Input data:
# - X can be either:
#     * a single numeric matrix (features x samples)
#     * a named list of numeric matrices (each features x samples). All must share samples.
# - 'methods' any subset of c("FABIA","MOFA","GFA","MFA")
#
# Returns a named list with one entry per method:
#   $<METHOD>$scores   : matrix (samples x k)
#   $<METHOD>$loadings : matrix (features(±view) x k)  (row names are features; if multiple views, rows are stacked with "view::feature")
#   $<METHOD>$meta     : list with method-specific bits
#
# You can set k (target #factors), scaling, and MFA group sizes.
#
# Required CRAN packages:
#   fabia, MOFA2, GFA, FactoMineR, Matrix
# --------------------------------------

# ---------------------------
# Helper: build comparison object (loadings per method, same feature order)
# ---------------------------
build_loading_comparison <- function(results) {
  stopifnot(length(results) >= 1)
  # union of feature names across methods
  feat_all <- unique(unlist(lapply(results, function(r) rownames(r$loadings))))
  out <- lapply(results, function(r) {
    L <- r$loadings
    # align to union (pad missing with 0)
    M <- matrix(0, nrow = length(feat_all), ncol = ncol(L),
                dimnames = list(feat_all, colnames(L)))
    M[rownames(L), ] <- L
    M
  })
  return(out)  # list: method -> (features x k)
}

run_factor_pipeline <- function(
    X,
    methods = c("FABIA","MOFA","GFA","MFA"),
    k = 10,
    scale_features = TRUE,
    mfa_group_sizes = NULL,
    verbose = TRUE
) {
  .stop_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is required.")
  .stack_views <- function(lst) {
    do.call(rbind, lapply(names(lst), function(v) {
      M <- lst[[v]]
      rn <- rownames(M); if (is.null(rn)) rn <- paste0("f", seq_len(nrow(M)))
      rownames(M) <- paste0(v, "::", rn)
      M
    }))
  }
  .as_list_views <- function(X) if (is.list(X)) X else list(omic1 = X)
  .nzcol_mask <- function(M) {
    # keep columns with non-zero variance (ignoring NAs)
    if (!is.matrix(M)) M <- as.matrix(M)
    if (ncol(M) == 0) return(rep(FALSE, 0))
    apply(M, 2, function(v) {
      vv <- v[is.finite(v)]
      if (!length(vv)) return(FALSE)
      sd(vv) > 0
    })
  }
  
  methods <- toupper(methods)
  
  # ---- input checks & harmonize ----
  if (!(is.matrix(X) || is.list(X))) stop("X must be a matrix or a list of matrices.")
  V <- .as_list_views(X)
  sample_names <- colnames(V[[1]])
  if (is.null(sample_names)) stop("All matrices must have column (sample) names.")
  for (nm in names(V)) {
    v <- V[[nm]]
    if (is.null(colnames(v))) stop("All matrices must have column (sample) names.")
    if (!setequal(colnames(v), sample_names))
      stop("All views must contain the same sample names.")
  }
  # reorder all views to same sample order
  V <- lapply(V, function(m) m[, sample_names, drop = FALSE])
  
  # optional scaling (feature-wise)
  if (scale_features) {
    V <- lapply(V, function(m) {
      m[] <- scale(t(m)) |> t()
      m[is.na(m)] <- 0
      m
    })
  }
  
  results <- list()
  
  # ----------------
  # FABIA (concat)
  # ----------------
  if ("FABIA" %in% methods) {
    .stop_pkg("fabia")
    if (verbose) message("Running FABIA ...")
    M_fabia <- .stack_views(V)   # features x samples
    p <- min(nrow(M_fabia), k)
    set.seed(21)#991489
    fit <- fabia(X = M_fabia, p = p, alpha = 0.01, cyc = 1000, spz = 0.5, lap = 1.0, nL = 2)
    # In fabia objects, loadings ~ L (features x k), factors/scores ~ Z (samples x k)
    
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
    results$FABIA <- list(scores = Z, loadings = L, meta = list(n_factors = ncol(Z)))
  }
  
  # -------------
  # MOFA2 (train with all features, let MOFA handle scaling)
  # -------------
  # if ("MOFA" %in% methods) {
  #   .stop_pkg("MOFA2")
  #   if (verbose) message("Running MOFA2 ...")
  #   library(MOFA2)
  #   
  #   # Ensure input matrices are numeric with rownames
  #   V_full <- lapply(V, function(m) {
  #     M <- as.matrix(m)
  #     storage.mode(M) <- "double"
  #     if (is.null(rownames(M))) rownames(M) <- paste0("f", seq_len(nrow(M)))
  #     M
  #   })
  #   
  #   # Build MOFA object
  #   mofa_obj <- MOFA2::create_mofa(V_full)
  #   
  #   # --- Data options: let MOFA do the scaling
  #   data_opts <- MOFA2::get_default_data_options(mofa_obj)
  #   data_opts$scale_views  <- TRUE
  #   data_opts$scale_groups <- TRUE
  #   
  #   # --- Model options
  #   model_opts <- MOFA2::get_default_model_options(mofa_obj)
  #   model_opts$num_factors <- k
  #   
  #   # --- Training options
  #   train_opts <- MOFA2::get_default_training_options(mofa_obj)
  #   #train_opts$verbose <- isTRUE(verbose)
  #   train_opts$seed    <- 21#12367
  #   train_opts$convergence_mode <- "slow"
  #   train_opts$maxiter <- 2000
  #   
  #   # --- Prepare and run MOFA
  #   mofa_obj <- MOFA2::prepare_mofa(
  #     object = mofa_obj,
  #     data_options     = data_opts,
  #     model_options    = model_opts,
  #     training_options = train_opts
  #   )
  # 
  #   outfile <- paste0(getwd(),"model_raw_pipeline.hdf5")
  #   #file.path(
  #   #   tempdir(),
  #   #   paste0("mofa_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".hdf5")
  #   # )
  #   
  #   fit_trained <- MOFA2::run_mofa(mofa_obj, outfile, use_basilisk = TRUE)
  #   
  #   # --- Load model (keep inactive factors)
  #   fit <- tryCatch({
  #     MOFA2::load_model(outfile, remove_inactive_factors = FALSE)
  #   }, error = function(e) {
  #     if (methods::is(fit_trained, "MOFA")) {
  #       warning("MOFA2: load_model() failed; using in-memory model. Inactive factors may be pruned.")
  #       fit_trained
  #     } else {
  #       stop("MOFA2: load_model() failed and run_mofa() did not return a MOFA object: ",
  #            conditionMessage(e))
  #     }
  #   })
  #   
  #   # --- Extract scores
  #   Z_list <- MOFA2::get_factors(fit, factors = "all", as.data.frame = FALSE)
  #   Z <- as.matrix(Z_list[[1]])
  #   if (!is.null(Z) && ncol(Z) > 0) {
  #     colnames(Z) <- paste0("F", seq_len(ncol(Z)))
  #   } else {
  #     Z <- matrix(0, nrow = length(sample_names), ncol = 0,
  #                 dimnames = list(sample_names, NULL))
  #   }
  #   
  #   # --- Extract loadings
  #   W_list <- MOFA2::get_weights(fit, views = "all", factors = "all", as.data.frame = FALSE)
  #   L <- .stack_views(W_list)
  #   if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
  #   
  #   results$MOFA <- list(
  #     scores   = Z,
  #     loadings = L,
  #     meta     = list(
  #       n_factors = ncol(Z),
  #       trained_hdf5 = outfile,
  #       scaling_in_mofa = TRUE
  #     )
  #   )
  # }
  if ("MOFA" %in% methods) {
    .stop_pkg("MOFA2")
    if (verbose) message("Running MOFA2 ...")
    library(MOFA2)
    
    # Ensure input matrices are numeric with rownames
    V_full <- lapply(V, function(m) {
      M <- as.matrix(m)
      storage.mode(M) <- "double"
      if (is.null(rownames(M))) rownames(M) <- paste0("f", seq_len(nrow(M)))
      M
    })
    
    # Build MOFA object
    mofa_obj <- MOFA2::create_mofa(V_full)
    
    # --- Data options: let MOFA do the scaling
    data_opts <- MOFA2::get_default_data_options(mofa_obj)
    data_opts$scale_views  <- TRUE
    data_opts$scale_groups <- TRUE
    
    # --- Model options
    model_opts <- MOFA2::get_default_model_options(mofa_obj)
    model_opts$num_factors <- k
    
    # --- Training options
    train_opts <- MOFA2::get_default_training_options(mofa_obj)
    train_opts$seed    <- 21
    train_opts$convergence_mode <- "slow"
    train_opts$maxiter <- 2000
    
    # --- Prepare and run MOFA
    mofa_obj <- MOFA2::prepare_mofa(
      object = mofa_obj,
      data_options     = data_opts,
      model_options    = model_opts,
      training_options = train_opts
    )
    
    # Safer path construction
    outfile <- file.path(getwd(), "model_raw_pipeline.hdf5")
    
    fit_trained <- MOFA2::run_mofa(mofa_obj, outfile, use_basilisk = TRUE)
    
    # --- Load model (keep inactive factors)
    fit <- tryCatch({
      MOFA2::load_model(outfile, remove_inactive_factors = FALSE)
    }, error = function(e) {
      if (methods::is(fit_trained, "MOFA")) {
        warning("MOFA2: load_model() failed; using in-memory model. Inactive factors may be pruned.")
        fit_trained
      } else {
        stop("MOFA2: load_model() failed and run_mofa() did not return a MOFA object: ",
             conditionMessage(e))
      }
    })
    
    # --- Extract scores
    Z_list <- MOFA2::get_factors(fit, factors = "all", as.data.frame = FALSE)
    Z <- as.matrix(Z_list[[1]])
    if (!is.null(Z) && ncol(Z) > 0) {
      colnames(Z) <- paste0("F", seq_len(ncol(Z)))
    } else {
      Z <- matrix(0, nrow = length(sample_names), ncol = 0,
                  dimnames = list(sample_names, NULL))
    }
    
    # --- Extract loadings
    W_list <- MOFA2::get_weights(fit, views = "all", factors = "all", as.data.frame = FALSE)
    L <- .stack_views(W_list)
    if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
    
    # --- New: variance explained (plot + table)
    p_varExp <- tryCatch(
      MOFA2::plot_variance_explained(fit, max_r2 = 15),
      error = function(e) { warning("plot_variance_explained failed: ", conditionMessage(e)); NULL }
    )
    varExp_tbl <- tryCatch(
      MOFA2::calculate_variance_explained(fit),
      error = function(e) NULL
    )
    
    # --- New: factor-by-group boxplots
    p_factorByGroup <- tryCatch({
      if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 not installed")
      # Get sample metadata from the fitted object if available
      smd <- tryCatch(MOFA2::samples_metadata(fit), error = function(e) NULL)
      # Build a minimal metadata if missing
      if (is.null(smd) || nrow(smd) == 0) {
        smd <- data.frame(sample = rownames(Z), Group = "All", row.names = rownames(Z))
      } else {
        # ensure rownames and a Group column
        if (is.null(rownames(smd)) && "sample" %in% names(smd)) {
          rownames(smd) <- smd$sample
        }
        smd <- smd[rownames(Z), , drop = FALSE]
        if (!("Group" %in% colnames(smd))) {
          # Try to pick a categorical column as Group; else fallback to "All"
          cat_cols <- names(which(vapply(smd, function(x) is.factor(x) || is.character(x), logical(1))))
          if (length(cat_cols) > 0) {
            smd$Group <- smd[[cat_cols[1]]]
          } else {
            smd$Group <- "All"
          }
        }
      }
      
      # Long table of factor scores
      df_long <- utils::stack(as.data.frame(Z))
      df_long$sample <- rep(rownames(Z), times = ncol(Z))
      names(df_long) <- c("score", "Factor", "sample")
      df_long$Factor <- as.character(df_long$Factor)
      df_long$Group  <- smd$Group[match(df_long$sample, rownames(smd))]
      
      ggplot2::ggplot(df_long, ggplot2::aes(x = Factor, y = score, fill = Group)) +
        ggplot2::geom_boxplot(outlier.size = 0.5) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "MOFA factor scores by group",
          x = "Factor",
          y = "Score",
          fill = "Group"
        )
    }, error = function(e) { warning("factor-by-group plot failed: ", conditionMessage(e)); NULL })
    
    results$MOFA <- list(
      scores   = Z,
      loadings = L,
      meta     = list(
        n_factors = ncol(Z),
        trained_hdf5 = outfile,
        scaling_in_mofa = TRUE,
        variance_explained = varExp_tbl
      ),
      plots = list(
        variance_by_factor = p_varExp,      # ggplot (may be NULL on error)
        factors_by_group   = p_factorByGroup # ggplot (may be NULL on error)
      )
    )
  }
  
  # -------------
  # GFA (multiview) — single run, fixed labeling & single-matrix W support
  # -------------
  if ("GFA" %in% methods) {
    .stop_pkg("GFA")
    if (verbose) message("Running GFA ...")
    library(GFA)
    
    ## Ensure views have stable names: 'rna' and 'prot' for 2 views, else view1, view2, ...
    if (is.null(names(V)) || any(!nzchar(names(V)))) {
      names(V) <- if (length(V) == 2) c("rna", "prot") else paste0("view", seq_along(V))
    }
    
    ## GFA expects samples x features; keep feature names
    V_gfa <- lapply(V, t)                 # each: samples x features
    names(V_gfa) <- names(V)              # keep labels ('rna','prot')
    
    ## Options (yours)
    model_option <- GFA::getDefaultOpts()
    model_option$iter.burnin <- 500
    model_option$iter.max    <- 1000
    
    ## Normalization (as you had)
    norm <- GFA::normalizeData(V_gfa, type = "center")
    
    ## Fit once
    fit <- GFA::gfa(norm$train, K = k, opts = model_option)
    
    ## ---------- SCORES (samples x k) ----------
    Z <- {
      Zcand <- if (is.list(fit$X) && length(fit$X) >= 1) fit$X[[1]] else fit$X
      if (is.null(Zcand)) {
        matrix(0, nrow = nrow(V_gfa[[1]]), ncol = 0,
               dimnames = list(rownames(V_gfa[[1]]), NULL))
      } else {
        Zcand <- as.matrix(Zcand)
        if (is.null(rownames(Zcand))) rownames(Zcand) <- rownames(V_gfa[[1]])
        if (NCOL(Zcand) > 0 && is.null(colnames(Zcand))) colnames(Zcand) <- paste0("F", seq_len(NCOL(Zcand)))
        Zcand
      }
    }
    
    ## ---------- LOADINGS (features x k), with 'rna::' / 'prot::' prefixes ----------
    make_empty_block <- function(i) {
      nf <- ncol(V_gfa[[i]])
      feats <- colnames(V_gfa[[i]]); if (is.null(feats)) feats <- paste0("feature_", seq_len(nf))
      M <- matrix(0, nrow = nf, ncol = 0)
      rownames(M) <- paste0(names(V_gfa)[i], "::", feats)
      M
    }
    
    W_list <- list()
    
    if (is.list(fit$W) && length(fit$W) == length(V_gfa)) {
      ## Standard case: one W per view
      W_list <- lapply(seq_along(fit$W), function(i) {
        Wi <- as.matrix(fit$W[[i]])
        if (is.null(dim(Wi))) Wi <- matrix(Wi, nrow = 1)
        
        n_feat <- ncol(V_gfa[[i]])               # features in this view
        feats  <- colnames(V_gfa[[i]]); if (is.null(feats)) feats <- paste0("feature_", seq_len(n_feat))
        
        ## Orient to (features x k)
        if (nrow(Wi) == n_feat) {
          # ok
        } else if (ncol(Wi) == n_feat) {
          Wi <- t(Wi)
        } else if (n_feat == 1 && (nrow(Wi) == 1 || ncol(Wi) == 1)) {
          if (nrow(Wi) != 1) Wi <- t(Wi)
        } else {
          # unreconcilable for this view → empty, correctly labeled
          return(make_empty_block(i))
        }
        
        rownames(Wi) <- paste0(names(V_gfa)[i], "::", feats)
        if (NCOL(Wi) > 0 && is.null(colnames(Wi))) colnames(Wi) <- paste0("F", seq_len(NCOL(Wi)))
        Wi
      })
      
    } else if (is.matrix(fit$W)) {
      ## All views merged into a single W matrix
      Wi_all <- as.matrix(fit$W)
      
      view_feat   <- vapply(V_gfa, ncol, integer(1))   # features per view
      total_feat  <- sum(view_feat)
      
      ## Orient to (total_features x k)
      if (nrow(Wi_all) == total_feat) {
        # ok
      } else if (ncol(Wi_all) == total_feat) {
        Wi_all <- t(Wi_all)
      } else if (nrow(Wi_all) == ncol(Wi_all) && nrow(Wi_all) < total_feat) {
        warning(sprintf("GFA returned square W (%d x %d); cannot map to %d features. Returning empty blocks.",
                        nrow(Wi_all), ncol(Wi_all), total_feat))
        W_list <- lapply(seq_along(V_gfa), make_empty_block)
      } else {
        warning(sprintf("Unexpected W shape: %d x %d (expected total features %d). Returning empty blocks.",
                        nrow(Wi_all), ncol(Wi_all), total_feat))
        W_list <- lapply(seq_along(V_gfa), make_empty_block)
      }
      
      ## If Wi_all is oriented correctly, split rows into per-view blocks
      if (length(W_list) == 0) {
        cuts <- c(0, cumsum(view_feat))
        W_list <- vector("list", length(V_gfa))
        for (i in seq_along(V_gfa)) {
          rows  <- (cuts[i] + 1):cuts[i + 1]
          Wi    <- Wi_all[rows, , drop = FALSE]            # features_i x k
          feats <- colnames(V_gfa[[i]]); if (is.null(feats)) feats <- paste0("feature_", seq_len(nrow(Wi)))
          rownames(Wi) <- paste0(names(V_gfa)[i], "::", feats)
          if (NCOL(Wi) > 0 && is.null(colnames(Wi))) colnames(Wi) <- paste0("F", seq_len(NCOL(Wi)))
          W_list[[i]] <- Wi
        }
      }
      
    } else {
      ## Unknown W type → labeled empty blocks
      W_list <- lapply(seq_along(V_gfa), make_empty_block)
    }
    
    L <- if (length(W_list)) do.call(rbind, W_list) else NULL
    
    results$GFA <- list(
      scores   = Z,
      loadings = L,
      meta     = list(
        n_factors   = NCOL(Z),
        iter.max    = model_option$iter.max,
        iter.burnin = model_option$iter.burnin,
        K           = k
      )
    )
  }
  
  # -------------
  # MFA (FactoMineR) — filter zero-variance, NULL/empty-safe; robust loadings from quanti.var$coord
  # -------------
  # if ("MFA" %in% methods) {
  #   .stop_pkg("FactoMineR")
  #   if (verbose) message("Running MFA ...")
  #   library(FactoMineR)
  #   # samples x features per view (optionally filter zero-variance columns)
  #   Vt <- lapply(V, t)
  #   nz_mask <- function(M) if (ncol(M)) apply(M, 2, function(v) sd(v[is.finite(v)]) > 0) else logical(0)
  #   Vt <- lapply(Vt, function(m) { keep <- nz_mask(m); m[, keep, drop = FALSE] })
  #   
  #   # bind blocks; compute (or accept) group sizes
  #   X_mfa <- do.call(cbind, Vt)                  # samples x all_features
  #   if (is.null(mfa_group_sizes)) mfa_group_sizes <- vapply(Vt, ncol, 1L)
  #   stopifnot(sum(mfa_group_sizes) == ncol(X_mfa))
  #   
  #   # data.frame is friendlier to FactoMineR
  #   X_df <- as.data.frame(X_mfa)
  #   if (!is.null(rownames(X_mfa))) rownames(X_df) <- rownames(X_mfa)
  #   
  #   fit <- FactoMineR::MFA(
  #     X_df,
  #     group = mfa_group_sizes,
  #     type  = rep("c", length(mfa_group_sizes)),
  #     ncp   = k,
  #     graph = FALSE
  #   )
  #   
  #   ## Scores (samples x q)
  #   Z <- fit$ind$coord
  #   if (is.null(Z)) {
  #     Z <- matrix(0, nrow = nrow(X_df), ncol = 0,
  #                 dimnames = list(rownames(X_df), NULL))
  #     q <- 0L
  #   } else {
  #     Z <- as.matrix(Z); q <- ncol(Z)
  #     if (q > 0) colnames(Z) <- paste0("F", seq_len(q))
  #   }
  #   
  #   ## Loadings (features x q)
  #   # 1) Try global var coords; 2) else use quanti.var$coord; 3) else empty
  #   L <- fit$var$coord
  #   if (is.null(L)) {
  #     Lq <- fit$quanti.var$coord
  #     if (!is.null(Lq)) {
  #       # keep only dimensions (Dim.1, Dim.2, ...) and coerce to matrix
  #       dim_cols <- grep("^Dim\\.", colnames(Lq))
  #       L <- as.matrix(Lq[, dim_cols, drop = FALSE])
  #       # rename columns to F1..Fq (q may be <= k)
  #       if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
  #       
  #       # Reorder rows to match our concatenated block order
  #       # Build desired order of base variable names (per group, in the same order as X_mfa)
  #       want_base <- unlist(lapply(seq_along(Vt), function(i) colnames(Vt[[i]])), use.names = FALSE)
  #       
  #       # Row names in L often look like "group.var" or just "var"
  #       row_base <- sub("^.*?\\.", "", rownames(L))  # strip optional "group."
  #       # Try to match per-group to avoid cross-group duplicates
  #       offs <- c(0, cumsum(mfa_group_sizes)); offs <- offs[-length(offs)]
  #       idx <- integer(0)
  #       for (i in seq_along(Vt)) {
  #         g_vars <- colnames(Vt[[i]])
  #         # prefer matches where original rowname already equals "group.var"
  #         # but fall back to base-name match
  #         rnames <- rownames(L)
  #         cand_group <- paste0(names(V)[i], ".", g_vars)
  #         m <- match(cand_group, rnames)
  #         m[is.na(m)] <- match(g_vars[is.na(m)], row_base)
  #         idx <- c(idx, m)
  #       }
  #       # guard: drop any NAs (shouldn’t happen unless names collided)
  #       keep_ok <- !is.na(idx)
  #       L <- L[idx[keep_ok], , drop = FALSE]
  #       # final nice rownames with view prefixes
  #       nice_rn <- unlist(lapply(seq_along(Vt), function(i) paste0(names(V)[i], "::", colnames(Vt[[i]]))),
  #                         use.names = FALSE)
  #       rownames(L) <- nice_rn[keep_ok]
  #     }
  #   } else {
  #     L <- as.matrix(L)
  #     if (ncol(L) > 0) {
  #       # trim to q if needed
  #       if (!is.null(q) && q > 0 && ncol(L) > q) L <- L[, seq_len(q), drop = FALSE]
  #       colnames(L) <- paste0("F", seq_len(ncol(L)))
  #     }
  #     # add view prefixes in our concatenated order
  #     if (is.null(colnames(X_mfa))) colnames(X_mfa) <- paste0("v", seq_len(ncol(X_mfa)))
  #     block_names <- rep(names(V), times = mfa_group_sizes)
  #     rownames(L) <- paste0(block_names, "::", colnames(X_mfa))
  #   }
  #   
  #   # If still NULL/empty, return a dimensionally correct empty matrix
  #   if (is.null(L)) {
  #     rn <- unlist(lapply(seq_along(Vt), function(i)
  #       paste0(names(V)[i], "::", colnames(Vt[[i]]))), use.names = FALSE)
  #     L <- matrix(0, nrow = length(rn), ncol = 0, dimnames = list(rn, NULL))
  #   }
  #   
  #   results$MFA <- list(scores = Z, loadings = L, meta = list(n_factors = ncol(Z)))
  # }
  if ("MFA" %in% methods) {
    .stop_pkg("FactoMineR")
    if (verbose) message("Running MFA ...")
    library(FactoMineR)
    
    # transpose to samples x features per view
    Vt <- lapply(V, t)
    
    # add tiny jitter to any zero-variance columns so nothing is dropped
    jitter_cols <- function(M, eps = 1e-8) {
      if (ncol(M) == 0) return(M)
      sds <- apply(M, 2, sd, na.rm = TRUE)
      bad <- which(!is.finite(sds) | sds == 0)
      if (length(bad)) {
        M[, bad] <- M[, bad, drop = FALSE] +
          matrix(rnorm(nrow(M) * length(bad), 0, eps), nrow(M))
      }
      M
    }
    Vt <- lapply(Vt, jitter_cols)
    
    # bind blocks; compute (or accept) group sizes
    X_mfa <- do.call(cbind, Vt)
    if (is.null(mfa_group_sizes)) mfa_group_sizes <- vapply(Vt, ncol, 1L)
    stopifnot(sum(mfa_group_sizes) == ncol(X_mfa))
    
    # data.frame is friendlier to FactoMineR
    X_df <- as.data.frame(X_mfa)
    if (!is.null(rownames(X_mfa))) rownames(X_df) <- rownames(X_mfa)
    
    fit <- FactoMineR::MFA(
      X_df,
      group = mfa_group_sizes,
      type  = rep("c", length(mfa_group_sizes)),
      ncp   = k,
      graph = FALSE
    )
    
    ## Scores
    Z <- fit$ind$coord
    if (is.null(Z)) {
      Z <- matrix(0, nrow = nrow(X_df), ncol = 0,
                  dimnames = list(rownames(X_df), NULL))
      q <- 0L
    } else {
      Z <- as.matrix(Z); q <- ncol(Z)
      if (q > 0) colnames(Z) <- paste0("F", seq_len(q))
    }
    
    ## Loadings
    L <- fit$var$coord
    if (is.null(L)) {
      Lq <- fit$quanti.var$coord
      if (!is.null(Lq)) {
        dim_cols <- grep("^Dim\\.", colnames(Lq))
        L <- as.matrix(Lq[, dim_cols, drop = FALSE])
        if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
        
        # reorder and prefix
        nice_rn <- unlist(lapply(seq_along(Vt),
                                 function(i) paste0(names(V)[i], "::", colnames(Vt[[i]]))),
                          use.names = FALSE)
        rownames(L) <- nice_rn
      }
    } else {
      L <- as.matrix(L)
      if (ncol(L) > 0 && !is.null(q) && q > 0 && ncol(L) > q)
        L <- L[, seq_len(q), drop = FALSE]
      colnames(L) <- paste0("F", seq_len(ncol(L)))
      block_names <- rep(names(V), times = mfa_group_sizes)
      rownames(L) <- paste0(block_names, "::", colnames(X_mfa))
    }
    
    if (is.null(L)) {
      rn <- unlist(lapply(seq_along(Vt),
                          function(i) paste0(names(V)[i], "::", colnames(Vt[[i]]))),
                   use.names = FALSE)
      L <- matrix(0, nrow = length(rn), ncol = 0, dimnames = list(rn, NULL))
    }
    
    results$MFA <- list(scores = Z, loadings = L,
                        meta = list(n_factors = ncol(Z)))
  }
  
  
  if (verbose) message("Done. Methods fit: ", paste(names(results), collapse = ", "))
  results
}

#===========================================
# Create Benchmark Datasets
#===========================================
set.seed(123)
benchmark_res <- run_factor_pipeline(
  list(drugs = as.matrix(drugs), methyl = as.matrix(methylation)),
  methods = c("fabia","mofa","gfa","mfa"),
  k = 2,
  scale_features = TRUE
)

#============================================
# Merge 'scores' with metadata by sample for each method in a benchmark object
#============================================
merge_scores_with_metadata <- function(bench, metadata, sample_col = "sample",
                                       methods = intersect(names(bench), c("FABIA","MOFA","GFA","MFA","GT")),
                                       join = c("left","inner")) {
  join <- match.arg(join)
  
  # ensure metadata has the join key and unique rows
  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!(sample_col %in% names(meta))) {
    if (!is.null(rownames(meta))) {
      meta[[sample_col]] <- rownames(meta)
    } else {
      stop("`metadata` must have a sample column named ", sample_col, " or rownames to use as sample IDs.")
    }
  }
  meta[[sample_col]] <- as.character(meta[[sample_col]])
  meta <- distinct(meta, .data[[sample_col]], .keep_all = TRUE)
  
  by_method <- list()
  for (m in methods) {
    if (is.null(bench[[m]]$scores)) next
    
    S <- bench[[m]]$scores
    # coerce to data.frame and add 'sample' column from rownames if needed
    Sdf <- as.data.frame(S, check.names = FALSE, stringsAsFactors = FALSE)
    if (!("sample" %in% names(Sdf))) {
      rn <- rownames(Sdf)
      if (is.null(rn)) stop("Scores for ", m, " have no rownames and no 'sample' column.")
      Sdf$sample <- rn
    }
    Sdf$sample <- as.character(Sdf$sample)
    
    # join: scores$sample -> metadata[[sample_col]]
    by_map <- setNames(sample_col, "sample")
    out <- if (join == "left") {
      left_join(Sdf, meta, by = by_map)
    } else {
      inner_join(Sdf, meta, by = by_map)
    }
    
    # useful ordering: sample, factors..., metadata...
    fac_cols <- grep("^F\\d+(\\s*\\(|$)", names(Sdf), value = TRUE)
    front <- c("sample", fac_cols)
    out <- out[, c(intersect(front, names(out)), setdiff(names(out), front)), drop = FALSE]
    out$method <- m
    by_method[[m]] <- out
  }
  
  all_long <- bind_rows(by_method)
  list(by_method = by_method, all_long = all_long)
}

library(dplyr)

# Attach metadata to scores for every method in a benchmark object
attach_scores_metadata <- function(bench,
                                   metadata,
                                   sample_col = "sample",
                                   methods = intersect(names(bench), c("FABIA","MOFA","GFA","MFA","GT")),
                                   store_name = "scores_with_meta",
                                   replace_scores = FALSE) {
  # --- prep metadata ---
  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!(sample_col %in% names(meta))) {
    if (!is.null(rownames(meta))) {
      meta[[sample_col]] <- rownames(meta)
    } else {
      stop("`metadata` must have a '", sample_col, "' column or rownames.")
    }
  }
  meta[[sample_col]] <- as.character(meta[[sample_col]])
  meta <- distinct(meta, .data[[sample_col]], .keep_all = TRUE)
  
  for (m in methods) {
    if (is.null(bench[[m]]$scores)) next
    
    # coerce scores to data frame; pull sample IDs from rownames if needed
    Sdf <- as.data.frame(bench[[m]]$scores, check.names = FALSE, stringsAsFactors = FALSE)
    if (!("sample" %in% names(Sdf))) {
      rn <- rownames(Sdf)
      if (is.null(rn)) stop("Scores for ", m, " have no rownames and no 'sample' column.")
      Sdf$sample <- rn
    }
    Sdf$sample <- as.character(Sdf$sample)
    
    # left-join by sample
    merged <- left_join(Sdf, meta, by = setNames(sample_col, "sample"))
    merged$method <- m
    
    if (isTRUE(replace_scores)) {
      # keep original numeric matrix under a backup name, replace scores with merged df
      bench[[m]]$scores_matrix <- bench[[m]]$scores
      bench[[m]]$scores <- merged
    } else {
      # store alongside the original scores matrix
      bench[[m]][[store_name]] <- merged
    }
  }
  bench
}

mergedx <- attach_scores_metadata(benchmark_res, CLL_metadata, sample_col = "sample")

#============================================
# Plot based on metadata
#============================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)   # unit()

# ---------------- utilities ----------------
resolve_factor_col <- function(df, method, factor_idx = 1) {
  cand  <- paste0("F", factor_idx, " (", method, ")")
  short <- paste0("F", factor_idx)
  if (cand  %in% names(df)) return(cand)
  if (short %in% names(df)) return(short)
  hit <- grep(paste0("^F", factor_idx, "(\\b|\\s*\\()"), names(df), value = TRUE)
  if (length(hit)) hit[1] else stop("Could not find F", factor_idx, " for ", method)
}

# robust palette builder: user overrides defaults; input can be NULL, vector, or list
build_palette <- function(pal = NULL) {
  default <- c("0"="#E37449", "1"="#00366C", "NA"="#999999")#
  if (is.null(pal)) return(default)
  if (is.list(pal)) pal <- unlist(pal, use.names = TRUE)
  pal <- as.character(pal)
  if (is.null(names(pal))) stop("Custom palette must be a *named* vector with names '0','1','NA'.")
  default[names(pal)] <- pal
  default
}

# dataframe for one method
scores_df_for_method <- function(bench, method, factor_idx = 1, group_col = "IGHV") {
  x <- bench[[method]][["scores_with_meta"]]
  if (is.null(x)) stop("scores_with_meta missing for ", method, ". Run attach_scores_metadata() first.")
  score_col <- resolve_factor_col(x, method, factor_idx)
  
  grp_raw <- if (group_col %in% names(x)) x[[group_col]] else NA
  grp_chr <- dplyr::case_when(
    is.na(grp_raw) ~ "NA",
    as.character(grp_raw) %in% c("0","1") ~ as.character(grp_raw),
    TRUE ~ "NA"
  )
  
  tibble(
    sample      = if ("sample" %in% names(x)) x$sample else seq_len(nrow(x)),
    .sample_idx = seq_len(nrow(x)),
    .score      = as.numeric(x[[score_col]]),
    .group      = grp_chr,
    method      = method
  )
}

# ---------------- panel ----------------
plot_scores_panel <- function(
    df, y_lab, panel_tag = NULL,
    tag_position = c("tl_in","none"),
    palette = NULL
) {
  tag_position <- match.arg(tag_position)
  pal <- build_palette(palette)
  df$.group <- factor(df$.group, levels = c("0","1","NA"))
  
  p <- ggplot(df, aes(x = .sample_idx, y = .score, fill = .group, color = .group)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(shape = 21, size = 2.7, stroke = 0.8, alpha = 0.95) +
    # ⬇️ one legend only (from fill); hide color legend
    scale_fill_manual(values = pal, drop = FALSE, name = ".group") +
    scale_color_manual(values = pal, drop = FALSE, guide = "none") +
    labs(x = "Samples", y = y_lab) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.x = element_text(size = 14, face = "plain", colour = "black"),#1f3765 #0072B2
      axis.title.y = element_text(size = 14, face = "plain", colour = "black"),
      plot.margin = margin(6, 6, 6, 6)
    )
  
  if (tag_position == "tl_in" && !is.null(panel_tag) && nzchar(panel_tag)) {
    p <- p + annotate("label", x = -Inf, y = Inf, label = panel_tag,
                      hjust = -0.1, vjust = 1.1, size = 4.5,
                      fill = "grey90", label.r = unit(0.15, "lines"))
  }
  p
}

# small plot used for OUTSIDE tag row
make_tag_plot <- function(tag) {
  if (is.null(tag) || !nzchar(tag)) return(patchwork::plot_spacer())
  ggplot() +
    annotate("label", x = 0, y = 0, label = tag,
             size = 4.5, fill = "grey90", label.r = unit(0.15, "lines")) +
    xlim(-1, 1) + ylim(-1, 1) + theme_void() +
    theme(plot.margin = margin(0, 6, 0, 6))
}
# # tag_position: "tl_in" (inside, top-left), "none" (no tag)
# plot_scores_panel <- function(
#     df, y_lab, panel_tag = NULL,
#     tag_position = c("tl_in","none"),
#     palette = NULL
# ) {
#   tag_position <- match.arg(tag_position)
#   pal <- build_palette(palette)
#   df$.group <- factor(df$.group, levels = c("0","1","NA"))
#   
#   p <- ggplot(df, aes(x = .sample_idx, y = .score, fill = .group, color = .group)) +
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     geom_point(shape = 21, size = 2.7, stroke = 0.8, alpha = 0.95) +
#     scale_fill_manual(values = pal, drop = FALSE) +
#     scale_color_manual(values = pal, drop = FALSE) +
#     labs(x = "Samples", y = y_lab) +
#     coord_cartesian(clip = "off") +
#     theme_minimal(base_size = 13) +
#     theme(
#       panel.grid.minor = element_blank(),
#       panel.grid.major.x = element_blank(),
#       axis.title.x = element_text(size = 14, face = "bold", colour = "#1f3765"),
#       axis.title.y = element_text(size = 14, face = "bold", colour = "#1f3765"),
#       plot.margin = margin(6, 6, 6, 6)
#     )
#   
#   if (tag_position == "tl_in" && !is.null(panel_tag) && nzchar(panel_tag)) {
#     p <- p + annotate(
#       "label", x = -Inf, y = Inf, label = panel_tag,
#       hjust = -0.1, vjust = 1.1, size = 4.5,
#       fill = "grey90", label.r = unit(0.15, "lines")
#     )
#   }
#   p
# }
# 
# # small plot used for OUTSIDE tag row
# make_tag_plot <- function(tag) {
#   if (is.null(tag) || !nzchar(tag)) return(patchwork::plot_spacer())
#   ggplot() +
#     annotate("label", x = 0, y = 0, label = tag,
#              size = 4.5, fill = "grey90", label.r = unit(0.15, "lines")) +
#     xlim(-1, 1) + ylim(-1, 1) + theme_void() +
#     theme(plot.margin = margin(0, 6, 0, 6))
# }

plot_scores_scatter_grid <- function(
    bench,
    factor_idx      = 1,
    group_col       = "IGHV",
    methods         = c("MOFA","FABIA","MFA","GFA"),
    panel_tags      = NULL,                 # e.g., c("A","B","C","D") or NULL
    tag_position    = c("tl_in","tl_out","none"),
    title_text      = "Factor scores (scatter plots)",
    palette         = NULL,                 # passed to plot_scores_panel()
    legend_position = "right",              # "right" or "bottom"
    title_height    = 0.15,                 # relative space for title row
    tag_row_height  = 0.10                  # relative space for outside tag row
) {
  stopifnot(length(methods) == 4)
  tag_position <- match.arg(tag_position)
  
  # local tag-plot helper (only used when tag_position == "tl_out")
  .make_tag_plot <- function(tag) {
    if (is.null(tag) || !nzchar(tag)) return(patchwork::plot_spacer())
    ggplot() +
      annotate("label", x = 0, y = 0, label = tag,
               size = 4.5, fill = "grey90", label.r = grid::unit(0.15, "lines")) +
      xlim(-1, 1) + ylim(-1, 1) + theme_void() +
      theme(plot.margin = margin(0, 6, 0, 6))
  }
  
  # data for the four panels (relies on your scores_df_for_method)
  dfs <- lapply(methods, function(m) scores_df_for_method(bench, m, factor_idx, group_col))
  names(dfs) <- methods
  
  # normalize panel tags
  if (is.null(panel_tags)) tags <- rep("", 4) else tags <- rep_len(panel_tags, 4)
  all_empty <- all(!nzchar(tags))
  
  # build panels (relies on your plot_scores_panel)
  pA <- plot_scores_panel(dfs[[methods[1]]], paste0("Score (", methods[1], ")"),
                          panel_tag = if (tag_position == "tl_in") tags[1] else NULL,
                          tag_position = if (tag_position == "tl_in") "tl_in" else "none",
                          palette = palette)
  pB <- plot_scores_panel(dfs[[methods[2]]], paste0("Score (", methods[2], ")"),
                          panel_tag = if (tag_position == "tl_in") tags[2] else NULL,
                          tag_position = if (tag_position == "tl_in") "tl_in" else "none",
                          palette = palette)
  pC <- plot_scores_panel(dfs[[methods[3]]], paste0("Score (", methods[3], ")"),
                          panel_tag = if (tag_position == "tl_in") tags[3] else NULL,
                          tag_position = if (tag_position == "tl_in") "tl_in" else "none",
                          palette = palette)
  pD <- plot_scores_panel(dfs[[methods[4]]], paste0("Score (", methods[4], ")"),
                          panel_tag = if (tag_position == "tl_in") tags[4] else NULL,
                          tag_position = if (tag_position == "tl_in") "tl_in" else "none",
                          palette = palette)
  
  # put the legend settings onto each panel (avoids the '&' operator)
  leg_theme <- theme(legend.position = legend_position,
                     legend.title = element_text(face = "bold"))
  pA <- pA + leg_theme; pB <- pB + leg_theme; pC <- pC + leg_theme; pD <- pD + leg_theme
  
  # title
  title_grob <- ggplot() +
    annotate("text", x = 0, y = 0, label = title_text,
             fontface = 2, size = 5.5, colour = "#6a1b9a") +
    theme_void() + theme(plot.margin = margin(0,0,6,0))
  
  # 2×2 grid of panels
  panel_grid <- (pA | pB) / (pC | pD)
  
  # assemble with optional outside tag row and collect guides → one legend
  if (tag_position == "tl_out" && !all_empty) {
    tag_row <- .make_tag_plot(tags[1]) | .make_tag_plot(tags[2]) |
      .make_tag_plot(tags[3]) | .make_tag_plot(tags[4])
    title_grob / tag_row / panel_grid +
      plot_layout(heights = c(title_height, tag_row_height, 1), guides = "collect")
  } else {
    title_grob / panel_grid +
      plot_layout(heights = c(title_height, 1), guides = "collect")
  }
}

p_nonex <- plot_scores_scatter_grid(
  mergedx, factor_idx = 1, group_col = "IGHV",
  methods = c("MOFA","FABIA","MFA","GFA"),
  panel_tags = NULL,
  tag_position = "none"
)
print(p_nonex)

#============================================
# Create scores dataframes
#============================================
benchmark_res = mergedx
# Scores
fabia_score_df <- data.frame(benchmark_res$FABIA$scores); fabia_score_df$sample <- rownames(fabia_score_df)
# rename columns that start with "F"
colnames(fabia_score_df) <- sub("^F(\\d+)$", "F\\1 (FABIA)", colnames(fabia_score_df))

mofa_score_df <- data.frame(benchmark_res$MOFA$scores); mofa_score_df$sample <- rownames(mofa_score_df)
# rename columns that start with "F"
colnames(mofa_score_df) <- sub("^F(\\d+)$", "F\\1 (MOFA)", colnames(mofa_score_df))

mfa_score_df <- data.frame(benchmark_res$MFA$scores); mfa_score_df$sample <- rownames(mfa_score_df)
# rename columns that start with "F"
colnames(mfa_score_df) <- sub("^F(\\d+)$", "F\\1 (MFA)", colnames(mfa_score_df))

gfa_score_df <- data.frame(benchmark_res$GFA$scores); gfa_score_df$sample <- rownames(gfa_score_df)
# rename columns that start with "F"
colnames(gfa_score_df) <- sub("^F(\\d+)$", "F\\1 (GFA)", colnames(gfa_score_df))

names(fabia_score_df)
names(mofa_score_df)
names(mfa_score_df)
names(gfa_score_df)

library(dplyr)
library(ggplot2)
library(reshape2)

# LOADINGS
fabia_loading_df <- data.frame(benchmark_res$FABIA$loading); fabia_loading_df$feature <- rownames(fabia_loading_df)
# rename columns that start with "F"
colnames(fabia_loading_df) <- sub("^F(\\d+)$", "F\\1 (FABIA)", colnames(fabia_loading_df))

mofa_loading_df <- data.frame(benchmark_res$MOFA$loading); mofa_loading_df$feature <- rownames(mofa_loading_df)
# rename columns that start with "F"
colnames(mofa_loading_df) <- sub("^F(\\d+)$", "F\\1 (MOFA)", colnames(mofa_loading_df))

mfa_loading_df <- data.frame(benchmark_res$MFA$loading); mfa_loading_df$feature <- rownames(mfa_loading_df)
# rename columns that start with "F"
colnames(mfa_loading_df) <- sub("^F(\\d+)$", "F\\1 (MFA)", colnames(mfa_loading_df))

gfa_loading_df <- data.frame(benchmark_res$GFA$loading); gfa_loading_df$feature <- rownames(gfa_loading_df)
# rename columns that start with "F"
colnames(gfa_loading_df) <- sub("^F(\\d+)$", "F\\1 (GFA)", colnames(gfa_loading_df))

names(fabia_loading_df)
names(mofa_loading_df)
names(mfa_loading_df)
names(gfa_loading_df)


# Helper function: split loading df by prefix
# split_by_prefix <- function(df, prefix_s16 = "^drugs", prefix_met = "^met", prefix_lip = "^lip") {
#   df_s16 <- df[grepl(prefix_s16, df$feature), , drop = FALSE]
#   df_met <- df[grepl(prefix_met, df$feature), , drop = FALSE]
#   df_lip <- df[grepl(prefix_lip, df$feature), , drop = FALSE]
#   list(s16 = df_s16, met = df_met, lip = df_lip)
# }
split_by_pattern <- function(df, column, patterns) {
  if (!column %in% names(df)) {
    stop("Specified column not found in data frame")
  }
  if (is.null(names(patterns))) {
    stop("patterns must be a named character vector or list")
  }
  
  lapply(patterns, function(pat) {
    df[grepl(pat, df[[column]]), , drop = FALSE]
  })
}

patterns <- c(
  drugs = "^drugs",
  methyl = "^methyl"
)

# Apply to each metho
fabia_split <- split_by_pattern(
  df = fabia_loading_df,
  column = "feature",
  patterns = patterns
)

mofa_split <- split_by_pattern(
  df = mofa_loading_df,
  column = "feature",
  patterns = patterns
)

mfa_split <- split_by_pattern(
  df = mfa_loading_df,
  column = "feature",
  patterns = patterns
)

gfa_split <- split_by_pattern(
  df = gfa_loading_df,
  column = "feature",
  patterns = patterns
)

library(dplyr)
library(ggplot2)
library(reshape2)

# SCORES VISUALIZATION
names(fabia_score_df)
names(mofa_score_df)
names(mfa_score_df)
names(gfa_score_df)

#### ============================================================
#### 1) Ground truth extraction (scores + per-omic loadings)
####    - Scores rows:  sample_1, sample_2, ...
####    - Loadings rows: "rna::omic1_feature_i", "prot::omic2_feature_i", ...
####    - Per-omic data.frames also carry a `feature` column for merges.
#### ============================================================

# ============================================================
# Robust benchmark for factor methods with differing #columns
# ============================================================
# ---------- Pairs-plot panels (base graphics) -------------  # NEW
.panel_hist <- function(x) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE)
  y <- if (max(h$counts) == 0) h$counts else h$counts / max(h$counts)
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], y, col = "grey85", border = "white")
  dx <- try(density(x, na.rm = TRUE), silent = TRUE)
  if (!inherits(dx, "try-error")) lines(dx$x, dx$y / max(dx$y), lwd = 1)
}

.panel_cor <- function(x, y) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r  <- suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs"))
  txt <- if (is.finite(r)) sprintf("%.2f", abs(r)) else "NA"
  text(0.5, 0.5, txt, cex = 2.4, font = 2)
}

.panel_scatter <- function(x, y) {
  points(x, y, pch = 8, cex = 0.5)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) >= 3) abline(lm(y ~ x), col = "red", lwd = 1)
}

.make_pairs_plot <- function(mat, main = "", oma_top = 2) {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(oma = c(0, 0, oma_top, 0))
  pairs(
    mat,
    diag.panel = .panel_hist,
    upper.panel = .panel_cor,
    lower.panel = .panel_scatter
  )
  if (nzchar(main)) mtext(main, outer = TRUE, line = 0.5, cex = 1.2, font = 2)
  recordPlot()
}

# Flip 'v' so that cor(ref, v) >= 0 (pairwise NA-safe). If cor is NA, do nothing.
.align_sign_vec <- function(ref, v) {
  r <- suppressWarnings(stats::cor(ref, v, use = "pairwise.complete.obs"))
  if (is.finite(r) && r < 0) -v else v
}

# ---------- Reference-based factor alignment ----------------  # NEW
# Map each method's factors to a reference method using Hungarian on |cor|.
.align_to_reference_scores <- function(scores_list, ref, match_factors_fun) {
  ref_mat <- scores_list[[ref]]
  ref_k   <- ncol(ref_mat)
  maps <- list()
  for (m in names(scores_list)) if (m != ref) {
    common <- intersect(rownames(ref_mat), rownames(scores_list[[m]]))
    A <- ref_mat[common, , drop = FALSE]
    B <- scores_list[[m]][common, , drop = FALSE]
    C <- stats::cor(A, B, use = "pairwise.complete.obs")
    ms <- match_factors_fun(C)
    map <- rep(NA_integer_, ref_k)
    if (!is.null(ms) && nrow(ms) > 0) map[ms$A] <- ms$B
    maps[[m]] <- map
  }
  maps
}

.align_to_reference_loadings <- function(loadings_list, ref, get_prefix, get_suffix, match_factors_fun) {
  ref_mat <- loadings_list[[ref]]
  ref_k   <- ncol(ref_mat)
  pref_ref <- get_prefix(rownames(ref_mat))
  ds <- unique(pref_ref)
  out <- setNames(vector("list", length(ds)), ds)
  for (om in ds) {
    idx_ref <- which(pref_ref == om)
    ids_ref <- get_suffix(rownames(ref_mat)[idx_ref])
    maps_om <- list()
    for (m in names(loadings_list)) if (m != ref) {
      pref_m <- get_prefix(rownames(loadings_list[[m]]))
      idx_m  <- which(pref_m == om)
      if (!length(idx_m)) { maps_om[[m]] <- rep(NA_integer_, ref_k); next }
      ids_m  <- get_suffix(rownames(loadings_list[[m]])[idx_m])
      common <- intersect(ids_ref, ids_m)
      if (!length(common)) { maps_om[[m]] <- rep(NA_integer_, ref_k); next }
      A <- ref_mat[idx_ref[match(common, ids_ref)], , drop = FALSE]
      B <- loadings_list[[m]][idx_m[match(common, ids_m)], , drop = FALSE]
      C <- stats::cor(A, B, use = "pairwise.complete.obs")
      ms <- match_factors_fun(C)
      map <- rep(NA_integer_, ref_k)
      if (!is.null(ms) && nrow(ms) > 0) map[ms$A] <- ms$B
      maps_om[[m]] <- map
    }
    out[[om]] <- maps_om
  }
  out
}

benchmark_factor_methods <- function(
    benchmark_res,
    methods = c("FABIA","MOFA","GFA","MFA","GT"),
    rename_cols = TRUE,          # label factor cols "F1 (METHOD)" etc.
    center_scale = TRUE,         # z-score each factor column
    feature_id_split = "::",     # used to align features (suffix after this)
    ground_truth_feature = NULL, # optional named 0/1 vector
    ground_truth_sample  = NULL, # optional named 0/1 vector
    plot_theme_base_size = 11){
  # ---- deps ----
  req <- c("ggplot2","reshape2")
  miss <- req[!vapply(req, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss)) stop("Please install required packages: ", paste(miss, collapse=", "))
  if (!requireNamespace("clue", quietly = TRUE)) stop("Please install package 'clue'")
  
  `%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
  
  # ---------- robust numeric ingestion ----------
  coerce_numeric_cols <- function(df) {
    df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
    f_like <- grepl("^F[0-9]+\\s*\\(", names(df))
    for (j in which(f_like)) if (!is.numeric(df[[j]])) df[[j]] <- suppressWarnings(as.numeric(df[[j]]))
    non_f <- which(!f_like)
    for (j in non_f) if (!is.numeric(df[[j]])) {
      tc <- suppressWarnings(type.convert(df[[j]], as.is = TRUE))
      if (is.numeric(tc)) df[[j]] <- tc
    }
    df
  }
  
  to_numeric_matrix <- function(obj, method_name, id_pref = c("sample","feature")) {
    df <- coerce_numeric_cols(obj)
    rn <- rownames(df)
    if (is.null(rn) || anyNA(rn) || any(rn == "")) {
      hit <- intersect(id_pref, names(df))
      if (length(hit)) rn <- as.character(df[[hit[1]]])
    }
    f_like <- grepl("^F[0-9]+\\s*\\(", names(df))
    keep <- f_like & vapply(df, is.numeric, logical(1))
    if (!any(keep)) {
      drop_names <- c("sample","feature", grep("^signal_", names(df), value = TRUE))
      keep <- vapply(df, is.numeric, logical(1)) & !names(df) %in% drop_names
    }
    if (!any(keep)) {
      stop("No numeric columns found for method ", method_name,
           ". Columns were: ", paste(names(df), collapse = ", "))
    }
    mat <- as.matrix(df[, keep, drop = FALSE])
    storage.mode(mat) <- "double"
    if (is.null(rn)) {
      rn <- if (identical(id_pref[1], "sample")) paste0("sample_", seq_len(nrow(mat)))
      else paste0("feature_", seq_len(nrow(mat)))
    }
    rownames(mat) <- rn
    mat
  }
  
  # robust scaling that never drops dimensions
  safe_scale <- function(X) {
    if (!center_scale) return(as.matrix(X))
    X <- as.matrix(X)
    if (ncol(X) == 0L) return(X)
    ok <- colSums(!is.na(X)) > 0
    if (!any(ok)) return(X)
    tmp <- scale(X[, ok, drop = FALSE])
    if (is.null(dim(tmp))) {
      tmp <- matrix(tmp, nrow = nrow(X), ncol = sum(ok),
                    dimnames = list(rownames(X), colnames(X)[ok]))
    }
    X[, ok] <- as.matrix(tmp)
    X
  }
  
  rename_factor_cols <- function(X, method_name) {
    k <- ncol(X)
    if (rename_cols) {
      colnames(X) <- paste0("F", seq_len(k), " (", method_name, ")")
    } else if (is.null(colnames(X))) {
      colnames(X) <- paste0("F", seq_len(k))
    }
    X
  }
  
  # ---- 0) Which methods are present ----
  available <- intersect(methods, names(benchmark_res))
  if (length(available) < 2)
    stop("Need at least two methods in 'benchmark_res'. Found: ", paste(available, collapse=", "))
  
  # ---- 1) Collect & standardize matrices ----
  get_scores <- function(m) {
    X <- to_numeric_matrix(benchmark_res[[m]]$scores, m, id_pref = "sample")
    X <- rename_factor_cols(X, m)
    X <- safe_scale(X)
    X
  }
  get_loadings <- function(m) {
    X <- to_numeric_matrix(benchmark_res[[m]]$loadings, m, id_pref = "feature")
    X <- rename_factor_cols(X, m)
    X <- safe_scale(X)
    X
  }
  
  scores_list   <- lapply(available, get_scores);   names(scores_list)   <- available
  loadings_list <- lapply(available, get_loadings); names(loadings_list) <- available
  
  # ---- 2) Helpers for per-dataset loadings --------------------------------
  get_prefix <- function(rn) sub("::.*$", "", rn)       # "rna", "prot", "omic3", ...
  get_suffix <- function(rn) sub("^.*::", "", rn)       # "omic1_feature_123", ...
  normalize_ids <- function(xnames) {
    if (is.null(xnames)) return(NULL)
    if (is.null(feature_id_split)) return(xnames)
    parts <- strsplit(xnames, split = feature_id_split, fixed = TRUE)
    vapply(parts, function(p) if (length(p)) tail(p, 1) else "", character(1))
  }
  
  # Per-dataset loading correlations
  cor_mat_loadings_peromic <- function(A, B) {
    pA <- get_prefix(rownames(A)); pB <- get_prefix(rownames(B))
    omics <- intersect(unique(pA), unique(pB))
    out <- vector("list", length(omics)); names(out) <- omics
    for (om in omics) {
      idxA <- which(pA == om); idxB <- which(pB == om)
      if (!length(idxA) || !length(idxB)) next
      idsA <- get_suffix(rownames(A)[idxA])
      idsB <- get_suffix(rownames(B)[idxB])
      common <- intersect(idsA, idsB)
      if (!length(common)) next
      AA <- A[idxA[match(common, idsA)], , drop = FALSE]
      BB <- B[idxB[match(common, idsB)], , drop = FALSE]
      out[[om]] <- stats::cor(AA, BB, use = "pairwise.complete.obs")
    }
    Filter(Negate(is.null), out)
  }
  
  # --- aligned vectors for scattering loadings within a dataset ---
  get_aligned_loading_vectors <- function(A, B, dataset, i, j) {
    pA <- get_prefix(rownames(A)); pB <- get_prefix(rownames(B))
    idxA <- which(pA == dataset); idxB <- which(pB == dataset)
    if (!length(idxA) || !length(idxB)) return(NULL)
    idsA <- get_suffix(rownames(A)[idxA]); idsB <- get_suffix(rownames(B)[idxB])
    common <- intersect(idsA, idsB)
    if (!length(common)) return(NULL)
    x <- A[idxA[match(common, idsA)], i]
    y <- B[idxB[match(common, idsB)], j]
    list(x = x, y = y)
  }
  
  # ---- plotting helpers (scatter with r) -----------------------------------
  make_scatter <- function(x, y, title, xlab, ylab) {
    df <- data.frame(x = as.numeric(x), y = as.numeric(y))
    r  <- suppressWarnings(stats::cor(df$x, df$y, use = "pairwise.complete.obs"))
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(alpha = 0.6, size = 1) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, linetype = "dashed") +
      ggplot2::annotate("text", x = Inf, y = Inf, label = sprintf("r = %.3f", r),
                        hjust = 1.05, vjust = 1.5, size = 4) +
      ggplot2::labs(title = title, x = xlab, y = ylab) +
      ggplot2::theme_minimal(base_size = plot_theme_base_size) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }
  
  # ---- 3) Factor matching (Hungarian on |corr|), NA-safe -------------------
  match_factors <- function(C) {
    if (is.null(C) || any(dim(C) == 0)) return(NULL)
    C[!is.finite(C)] <- NA_real_
    absC <- abs(C); absC[is.na(absC)] <- 0
    ka <- nrow(absC); kb <- ncol(absC)
    maxv <- if (all(absC == 0)) 0 else max(absC, na.rm = TRUE)
    cost <- maxv - absC
    pad_val <- if (is.finite(suppressWarnings(max(cost, na.rm = TRUE))))
      suppressWarnings(max(cost, na.rm = TRUE)) + 1 else 1
    if (ka > kb) cost <- cbind(cost, matrix(pad_val, nrow = ka, ncol = ka - kb))
    if (kb > ka) cost <- rbind(cost, matrix(pad_val, nrow = kb - ka, ncol = kb))
    storage.mode(cost) <- "double"
    ass <- clue::solve_LSAP(cost)
    real_rows <- seq_len(ka)
    j_all     <- as.integer(ass[real_rows])
    keep      <- which(j_all <= kb)
    if (!length(keep)) return(NULL)
    i_keep <- real_rows[keep]; j_keep <- j_all[keep]
    data.frame(A = i_keep, B = j_keep, corr = mapply(function(i,j) C[i,j], i_keep, j_keep),
               stringsAsFactors = FALSE)
  }
  
  # ---- 4) Pairwise similarities (+ scatters) -------------------------------
  method_pairs <- t(combn(available, 2))
  pair_summaries_scores   <- list()
  pair_summaries_loadings <- list()
  pair_heatmaps_scores    <- list()
  pair_heatmaps_loadings  <- list()
  pair_scatter_scores     <- list()   # NEW
  pair_scatter_loadings   <- list()   # NEW
  
  for (r in seq_len(nrow(method_pairs))) {
    m1 <- method_pairs[r,1]; m2 <- method_pairs[r,2]
    
    # ----- SCORES -----
    common_samp <- intersect(rownames(scores_list[[m1]]), rownames(scores_list[[m2]]))
    A_s <- scores_list[[m1]][common_samp, , drop = FALSE]
    B_s <- scores_list[[m2]][common_samp, , drop = FALSE]
    Cs  <- stats::cor(A_s, B_s, use = "pairwise.complete.obs")
    ms  <- match_factors(Cs)
    
    # summary + heatmap
    pair_summaries_scores[[paste(m1,m2,sep="|")]] <- list(
      summary = data.frame(
        pair = paste(m1, m2, sep = " vs "),
        k_m1_scores = ncol(A_s),
        k_m2_scores = ncol(B_s),
        k_matched_scores = if (!is.null(ms)) nrow(ms) else 0L,
        mean_abs_corr_scores = if (!is.null(ms)) mean(abs(ms$corr), na.rm = TRUE) else NA_real_,
        median_abs_corr_scores = if (!is.null(ms)) stats::median(abs(ms$corr), na.rm = TRUE) else NA_real_,
        max_abs_corr_scores = if (!is.null(ms)) max(abs(ms$corr), na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      ),
      scores_cor = Cs, scores_match = ms
    )
    
    dfh <- reshape2::melt(Cs,
                          varnames = c(paste0("Factor (", m1,")"),
                                       paste0("Factor (", m2,")")),
                          value.name = "Correlation")
    xcol <- names(dfh)[2]; ycol <- names(dfh)[1]
    pair_heatmaps_scores[[paste(m1,m2,sep="|")]] <-
      ggplot2::ggplot(dfh, ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]], fill = .data[["Correlation"]])) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(limits = c(-1,1), midpoint = 0, name = "Corr") +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = paste0("Scores correlation: ", m1, " vs ", m2), x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = plot_theme_base_size) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold"))
    
    # score scatters for matched pairs
    ss_plots <- list()
    if (!is.null(ms) && nrow(ms) > 0) {
      for (k in seq_len(nrow(ms))) {
        i <- ms$A[k]; j <- ms$B[k]
        x <- A_s[, i]; y <- B_s[, j]
        title <- sprintf("Scores scatter: %s[%s] vs %s[%s]",
                         m1, colnames(A_s)[i], m2, colnames(B_s)[j])
        ss_plots[[k]] <- make_scatter(x, y, title,
                                      xlab = paste0(colnames(A_s)[i]),
                                      ylab = paste0(colnames(B_s)[j]))
      }
    }
    pair_scatter_scores[[paste(m1,m2,sep="|")]] <- ss_plots
    
    # ----- LOADINGS (per dataset) -----
    Cl_list <- cor_mat_loadings_peromic(loadings_list[[m1]], loadings_list[[m2]])
    ml_list <- lapply(Cl_list, match_factors)
    
    # per-dataset summary + heatmap + scatters
    rows <- list()
    scat_by_ds <- list()
    for (om in names(Cl_list)) {
      ml <- ml_list[[om]]
      rows[[length(rows)+1]] <- data.frame(
        pair = paste(m1, m2, sep = " vs "),
        dataset = om,
        k_m1_loadings = ncol(loadings_list[[m1]]),
        k_m2_loadings = ncol(loadings_list[[m2]]),
        k_matched_loadings = if (!is.null(ml)) nrow(ml) else 0L,
        mean_abs_corr_loadings = if (!is.null(ml)) mean(abs(ml$corr), na.rm = TRUE) else NA_real_,
        median_abs_corr_loadings = if (!is.null(ml)) stats::median(abs(ml$corr), na.rm = TRUE) else NA_real_,
        max_abs_corr_loadings = if (!is.null(ml)) max(abs(ml$corr), na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
      
      # heatmap for this dataset
      dfhL <- reshape2::melt(Cl_list[[om]],
                             varnames = c(paste0("Factor (", m1,")"),
                                          paste0("Factor (", m2,")")),
                             value.name = "Correlation")
      xL <- names(dfhL)[2]; yL <- names(dfhL)[1]
      pair_heatmaps_loadings[[paste(m1,m2,om,sep="|")]] <-
        ggplot2::ggplot(dfhL, ggplot2::aes(x = .data[[xL]], y = .data[[yL]], fill = .data[["Correlation"]])) +
        ggplot2::geom_tile(color = "white", linewidth = 0.2) +
        ggplot2::scale_fill_gradient2(limits = c(-1,1), midpoint = 0, name = "Corr") +
        ggplot2::coord_fixed() +
        ggplot2::labs(title = paste0("Loadings correlation (", om, "): ", m1, " vs ", m2),
                      x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = plot_theme_base_size) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       panel.grid = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(face = "bold"))
      
      # loading scatters for matched pairs within this dataset
      ds_plots <- list()
      if (!is.null(ml) && nrow(ml) > 0) {
        for (k in seq_len(nrow(ml))) {
          i <- ml$A[k]; j <- ml$B[k]
          vecs <- get_aligned_loading_vectors(loadings_list[[m1]], loadings_list[[m2]], om, i, j)
          if (is.null(vecs)) next
          title <- sprintf("Loadings scatter (%s): %s[%s] vs %s[%s]",
                           om, m1, colnames(loadings_list[[m1]])[i], m2, colnames(loadings_list[[m2]])[j])
          ds_plots[[k]] <- make_scatter(vecs$x, vecs$y, title,
                                        xlab = paste0(colnames(loadings_list[[m1]])[i]),
                                        ylab = paste0(colnames(loadings_list[[m2]])[j]))
        }
      }
      scat_by_ds[[om]] <- ds_plots
    }
    
    pair_summaries_loadings[[paste(m1,m2,sep="|")]] <- list(
      cor   = Cl_list,
      match = ml_list,
      summary = if (length(rows)) do.call(rbind, rows) else NULL
    )
    pair_scatter_loadings[[paste(m1,m2,sep="|")]] <- scat_by_ds
  }
  
  # ======================== 5b) Pairs-style plots (NEW) =========================
  # Reference method for alignment
  ref_method <- if ("FABIA" %in% available) "FABIA" else available[1]
  
  # ---- Scores: build a pairs plot per reference factor -------------------------
  maps_scores <- .align_to_reference_scores(scores_list, ref_method, match_factors)
  
  # samples common to all selected methods
  common_samples_all <- Reduce(intersect, lapply(scores_list, rownames))
  pairs_scores_by_factor <- list()
  ref_k <- ncol(scores_list[[ref_method]])
  
  # for (i in seq_len(ref_k)) {
  #   cols <- list()
  #   for (m in available) {
  #     if (m == ref_method) {
  #       cols[[paste0("Score (", m, ")")]] <- scores_list[[m]][common_samples_all, i]
  #     } else {
  #       j <- maps_scores[[m]][i]
  #       if (is.na(j)) next
  #       cols[[paste0("Score (", m, ")")]] <- scores_list[[m]][common_samples_all, j]
  #     }
  #   }
  #   if (length(cols) >= 2) {
  #     mat <- as.data.frame(cols, check.names = FALSE)
  #     pairs_scores_by_factor[[paste0("F", i)]] <-
  #       .make_pairs_plot(mat, main = sprintf("Factor scores (cor.) — aligned to %s F%d", ref_method, i))
  #   }
  # }
  for (i in seq_len(ref_k)) {
    cols <- list()
    ref_vec <- scores_list[[ref_method]][common_samples_all, i]
    
    for (m in available) {
      if (m == ref_method) {
        cols[[paste0("Score (", m, ")")]] <- ref_vec
      } else {
        j <- maps_scores[[m]][i]
        if (is.na(j)) next
        v <- scores_list[[m]][common_samples_all, j]
        v <- .align_sign_vec(ref_vec, v)        # << align sign to reference
        cols[[paste0("Score (", m, ")")]] <- v
      }
    }
    if (length(cols) >= 2) {
      mat <- as.data.frame(cols, check.names = FALSE)
      pairs_scores_by_factor[[paste0("F", i)]] <-
        .make_pairs_plot(mat, main = sprintf("Factor scores (Pearson cor.)"))
    }
  }
  
  # ---- Per-omic loadings: pairs plot per dataset & factor ----------------------
  maps_loadings <- .align_to_reference_loadings(loadings_list, ref_method, get_prefix, get_suffix, match_factors)
  pairs_loadings_by_dataset_and_factor <- list()
  
  # datasets seen in reference loadings
  datasets_ref <- unique(sub("::.*$", "", rownames(loadings_list[[ref_method]])))
  
  for (om in datasets_ref) {
    pref_ref <- get_prefix(rownames(loadings_list[[ref_method]]))
    idx_ref  <- which(pref_ref == om)
    ids_ref  <- get_suffix(rownames(loadings_list[[ref_method]])[idx_ref])
    
    # universal common features across all methods for this dataset
    common_ids <- ids_ref
    for (m in available) if (m != ref_method) {
      pref_m <- get_prefix(rownames(loadings_list[[m]]))
      ids_m  <- get_suffix(rownames(loadings_list[[m]])[pref_m == om])
      common_ids <- intersect(common_ids, ids_m)
    }
    if (!length(common_ids)) next
    
    # row indices by method for the same features in the same order
    row_idx <- list()
    row_idx[[ref_method]] <- idx_ref[match(common_ids, ids_ref)]
    for (m in available) if (m != ref_method) {
      pref_m <- get_prefix(rownames(loadings_list[[m]]))
      idx_m  <- which(pref_m == om)
      ids_m  <- get_suffix(rownames(loadings_list[[m]])[idx_m])
      row_idx[[m]] <- idx_m[match(common_ids, ids_m)]
    }
    
    # per-factor pairs
    ref_k <- ncol(loadings_list[[ref_method]])
  #   for (i in seq_len(ref_k)) {
  #     cols <- list()
  #     for (m in available) {
  #       if (m == ref_method) {
  #         cols[[paste0(m)]] <- loadings_list[[m]][row_idx[[m]], i]
  #       } else {
  #         j <- maps_loadings[[om]][[m]][i]
  #         if (is.na(j)) next
  #         cols[[paste0(m)]] <- loadings_list[[m]][row_idx[[m]], j]
  #       }
  #     }
  #     if (length(cols) >= 2) {
  #       mat <- as.data.frame(cols, check.names = FALSE)
  #       lab <- paste0("Loadings (", om, ") cor. — aligned to ", ref_method, " F", i)
  #       pairs_loadings_by_dataset_and_factor[[om]][[paste0("F", i)]] <- .make_pairs_plot(mat, main = lab)
  #     }
  #   }

    for (i in seq_len(ref_k)) {
      cols <- list()
      ref_vec <- loadings_list[[ref_method]][row_idx[[ref_method]], i]  # same features/order
      
      for (m in available) {
        if (m == ref_method) {
          cols[[paste0(m)]] <- ref_vec
        } else {
          j <- maps_loadings[[om]][[m]][i]
          if (is.na(j)) next
          v <- loadings_list[[m]][row_idx[[m]], j]
          v <- .align_sign_vec(ref_vec, v)        # << align sign to reference
          cols[[paste0(m)]] <- v
        }
      }
      if (length(cols) >= 2) {
        mat <- as.data.frame(cols, check.names = FALSE)
        lab <- paste0("Loadings (", om, ") — Pearson cor. ")
        pairs_loadings_by_dataset_and_factor[[om]][[paste0("F", i)]] <- .make_pairs_plot(mat, main = lab)
      }
    }
  }
    
  # ---- 6) Global all-in-one correlation heatmap (scores) ----
  merge_scores <- function(lst){
    X <- lst[[1]]; colnames(X) <- make.unique(colnames(X))
    for (i in seq(2, length(lst))) {
      Y <- lst[[i]]; colnames(Y) <- make.unique(colnames(Y))
      common <- intersect(rownames(X), rownames(Y))
      X <- cbind(X[common,,drop=FALSE], Y[common,,drop=FALSE])
    }
    X
  }
  all_scores_mat <- merge_scores(scores_list)
  all_cor <- stats::cor(all_scores_mat, use = "pairwise.complete.obs")
  
  get_method <- function(nm) sub(".*\\(([^)]+)\\)\\s*$", "\\1", nm)
  get_fnum   <- function(nm) suppressWarnings(as.integer(sub("^\\s*F(\\d+).*", "\\1", nm)))
  mth <- vapply(colnames(all_cor), get_method, character(1))
  fno <- vapply(colnames(all_cor), get_fnum, numeric(1))
  lev <- unique(mth)
  ord <- order(match(mth, lev), fno, colnames(all_cor))
  all_cor <- all_cor[ord, ord, drop = FALSE]
  
  dfall <- reshape2::melt(all_cor, varnames = c("Factor1","Factor2"), value.name = "Correlation")
  dfall$Factor1 <- factor(dfall$Factor1, levels = rownames(all_cor))
  dfall$Factor2 <- factor(dfall$Factor2, levels = colnames(all_cor))
  
  tab <- as.data.frame(table(mth[ord]), stringsAsFactors = FALSE)
  tab$cum <- cumsum(tab$Freq); bounds <- tab$cum + 0.5
  
  p_all <- ggplot2::ggplot(dfall, ggplot2::aes(x = Factor2, y = Factor1, fill = Correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(limits = c(-1,1), midpoint = 0, name = "Corr") +
    ggplot2::geom_vline(xintercept = bounds, linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = bounds, linewidth = 0.3) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "All methods — factor correlation (scores)", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = plot_theme_base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))
  
  # ---- 7) Optional AUCs vs truth (unchanged) ----
  normalize_ids_suffix <- function(xnames) normalize_ids(xnames)
  auc_results <- NULL
  if (!is.null(ground_truth_feature) || !is.null(ground_truth_sample)) {
    if (!requireNamespace("pROC", quietly = TRUE)) stop("To compute AUC, please install 'pROC'.")
    calc_auc <- function(scores, truth) {
      truth <- as.numeric(truth)
      if (!length(scores) || !length(truth) || length(scores) != length(truth)) return(NA_real_)
      roc <- pROC::roc(response = truth, predictor = scores, quiet = TRUE, direction = "<")
      as.numeric(pROC::auc(roc))
    }
    rows_auc <- list()
    if (!is.null(ground_truth_feature)) {
      for (m in available) {
        L <- loadings_list[[m]]
        s <- apply(abs(L), 1, max, na.rm = TRUE)
        idsL <- normalize_ids_suffix(names(s))
        idsT <- names(ground_truth_feature)
        if (!is.null(idsT)) idsT <- normalize_ids_suffix(idsT)
        common <- intersect(idsL, idsT)
        if (length(common)) {
          rows_auc[[length(rows_auc)+1]] <- data.frame(
            target="feature", method=m,
            AUC=calc_auc(scores = s[match(common, idsL)],
                         truth  = ground_truth_feature[match(common, idsT)]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (!is.null(ground_truth_sample)) {
      for (m in available) {
        S <- scores_list[[m]]
        ss <- apply(abs(S), 1, max, na.rm = TRUE)
        common <- intersect(names(ss), names(ground_truth_sample))
        if (length(common)) {
          rows_auc[[length(rows_auc)+1]] <- data.frame(
            target="sample", method=m,
            AUC=calc_auc(scores = ss[common], truth = ground_truth_sample[common]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    auc_results <- if (length(rows_auc)) do.call(rbind, rows_auc) else NULL
  }
  
  summary_table_scores   <- do.call(rbind, lapply(pair_summaries_scores, `[[`, "summary"))
  rownames(summary_table_scores) <- NULL
  summary_table_loadings <- do.call(rbind, lapply(pair_summaries_loadings, `[[`, "summary"))
  rownames(summary_table_loadings) <- NULL
  
  list(
    available_methods        = available,
    scores                   = scores_list,
    loadings                 = loadings_list,
    pairwise_scores          = pair_summaries_scores,
    pairwise_loadings        = pair_summaries_loadings,  # per-dataset
    pair_heatmaps_scores     = pair_heatmaps_scores,
    pair_heatmaps_loadings   = pair_heatmaps_loadings,   # per dataset
    pair_scatter_scores      = pair_scatter_scores,      # << NEW
    pair_scatter_loadings    = pair_scatter_loadings,    # << NEW (list per dataset)
    all_scores_heatmap       = p_all,
    summary_table_scores     = summary_table_scores,
    summary_table_loadings   = summary_table_loadings,
    auc_results              = auc_results,
    pairs_scores_by_factor               = pairs_scores_by_factor,           # NEW
    pairs_loadings_by_dataset_and_factor = pairs_loadings_by_dataset_and_factor  # NEW
    
  )
}

benchx <- benchmark_factor_methods(benchmark_res)


#============================================
# Create scores dataframes
#============================================

