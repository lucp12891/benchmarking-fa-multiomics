# ------------- COSINE SIMILARITY SAMPLE SCORES ------------------------
#--- your objects (assuming they exist exactly as you showed) ---
fabia_sample_scores <- data.frame(data_scores$score_FABIA1)
mofa_sample_scores  <- data.frame(data_scores$score_MOFA1)
mfa_sample_scores   <- data.frame(data_scores$score_MFA1)
gfa_sample_scores   <- data.frame(data_scores$score_GFA2)

colnames(fabia_sample_scores) <- c("F1 (FABIA)")
colnames(mofa_sample_scores)  <- c("F1 (MOFA)")
colnames(mfa_sample_scores)   <- c("F1 (MFA)")
colnames(gfa_sample_scores)   <- c("F1 (GFA)")

score_all_rad = data.frame(
  sample = data_scores$samples,
  F1_FABIA = fabia_sample_scores$`F1 (FABIA)`,
  F1_MOFA = mofa_sample_scores$`F1 (MOFA)`,
  F1_MFA = mfa_sample_scores$`F1 (MFA)`,
  F1_GFA = gfa_sample_scores$`F1 (GFA)`
)


align_to_reference <- function(df, ref_col, target_cols = NULL,
         cor_method = "pearson") {
  stopifnot(ref_col %in% names(df))
  num_ok <- vapply(df, is.numeric, logical(1))
  if (is.null(target_cols)) {
    target_cols <- setdiff(names(df)[num_ok], ref_col)
  } else {
    target_cols <- intersect(target_cols, names(df)[num_ok])
  }
  r_ref <- df[[ref_col]]
  flips <- setNames(logical(length(target_cols)), target_cols)
  cors  <- setNames(rep(NA_real_, length(target_cols)), target_cols)
  
  for (nm in target_cols) {
    ok   <- is.finite(r_ref) & is.finite(df[[nm]])
    if (sum(ok) >= 3) {
      r <- suppressWarnings(stats::cor(r_ref[ok], df[[nm]][ok], method = cor_method))
      cors[nm] <- r
      if (is.finite(r) && r < 0) {
        df[[nm]] <- -df[[nm]]
        flips[nm] <- TRUE
      }
    }
  }
  list(data = df, flipped = flips, corr = cors)
}

# flip all to match FABIA direction first:
score_all_rad_aligned <- align_to_reference(
  score_all_rad, ref_col = "F1_FABIA",
  target_cols = c("F1_MOFA","F1_MFA","F1_GFA")
)$data

colnames(score_all_rad_aligned) <- c("sample", "FABIA", "MOFA", "MFA", "GFA")

# reorder: feature, FABIA, MOFA, GFA, MFA
score_all_rad_aligned <- score_all_rad_aligned[
  , c("sample", "FABIA", "MOFA", "GFA", "MFA")
]

# --- cosine similarity ---
cosine_sim <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y)) stop("length mismatch")
  nx <- sqrt(sum(x^2)); ny <- sqrt(sum(y^2))
  if (nx == 0 || ny == 0) stop("cosine undefined for zero vector")
  sum(x * y) / (nx * ny)
}

# --- extract the score matrix (drop sample column) ---
M <- as.matrix(score_all_rad_aligned[, c("FABIA","MOFA","GFA","MFA")])

# --- pairwise cosine similarity matrix ---
methods <- colnames(M)
cos_mat <- outer(seq_along(methods), seq_along(methods),
                 Vectorize(function(i, j) cosine_sim(M[, i], M[, j])))
dimnames(cos_mat) <- list(methods, methods)

cos_mat

# keep lower triangle (including diagonal)
cos_lower <- cos_mat
cos_lower[upper.tri(cos_lower)] <- NA

# round to two decimals
cos_lower <- round(cos_lower, 2)

library(ggplot2)

df_lower <- as.data.frame(as.table(cos_lower))
colnames(df_lower) <- c("Method1", "Method2", "Cosine")

ggplot(df_lower, aes(Method1, Method2, fill = Cosine)) +
  geom_tile(color = "white", linewidth = 0.3, na.rm = TRUE) +
  geom_text(
    aes(label = ifelse(is.na(Cosine), "", sprintf("%.2f", Cosine))),
    color = "white",
    size = 4
  ) +
  coord_fixed() +
  scale_fill_gradient(
    low = "red",
    high = "#08306B",
    limits = c(0, 1),
    na.value = "white"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  )

# ------------- COSINE SIMILARITY MRNA DATA ------------------------

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
  mRNA = "mRNA",
  proteomics = "proteomics"
)

# Apply to each metho
loading_split_rad <- split_by_pattern(
  df = data_weights,
  column = "view",
  patterns = patterns
)


proteomics_rad_tbl = data.frame(
  feature = loading_split_rad$proteomics$feature,
  F1_FABIA = loading_split_rad$proteomics$loading_FABIA1,
  F1_MOFA = loading_split_rad$proteomics$loading_MOFA1,
  F1_MFA = loading_split_rad$proteomics$loading_MFA1,
  F1_GFA = loading_split_rad$proteomics$loading_GFA1
)

# flip all to match FABIA direction first:
proteomics_rad_aligned <- align_to_reference(
  proteomics_rad_tbl, ref_col = "F1_FABIA",
  target_cols = c("F1_MOFA","F1_MFA","F1_GFA")
)$data

colnames(proteomics_rad_aligned) <- c("feature", "FABIA", "MOFA", "MFA", "GFA")
# reorder: feature, FABIA, MOFA, GFA, MFA
proteomics_rad_aligned <- proteomics_rad_aligned[
  , c("feature", "FABIA", "MOFA", "GFA", "MFA")
]

# --- cosine similarity ---
cosine_sim <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y)) stop("length mismatch")
  nx <- sqrt(sum(x^2)); ny <- sqrt(sum(y^2))
  if (nx == 0 || ny == 0) stop("cosine undefined for zero vector")
  sum(x * y) / (nx * ny)
}

# --- extract the score matrix (drop sample column) ---
proteomicsM <- as.matrix(proteomics_rad_aligned[, c("FABIA","MOFA","GFA","MFA")])

# --- pairwise cosine similarity matrix ---
methods <- colnames(proteomicsM)
proteomics_cos_mat <- outer(seq_along(methods), seq_along(methods),
                 Vectorize(function(i, j) cosine_sim(proteomicsM[, i], proteomicsM[, j])))
dimnames(proteomics_cos_mat) <- list(methods, methods)

proteomics_cos_mat

# keep lower triangle (including diagonal)
proteomics_cos_lower <- proteomics_cos_mat
proteomics_cos_lower[upper.tri(proteomics_cos_lower)] <- NA

# round to two decimals
proteomics_cos_lower <- round(proteomics_cos_lower, 2)

library(ggplot2)

df_lower <- as.data.frame(as.table(proteomics_cos_lower))
colnames(df_lower) <- c("Method1", "Method2", "Cosine")

ggplot(df_lower, aes(Method1, Method2, fill = Cosine)) +
  geom_tile(color = "white", linewidth = 0.3, na.rm = TRUE) +
  geom_text(
    aes(label = ifelse(is.na(Cosine), "", sprintf("%.2f", Cosine))),
    color = "white",
    size = 4
  ) +
  coord_fixed() +
  scale_fill_gradient(
    low = "red",
    high = "#08306B",
    limits = c(0, 1),
    na.value = "white"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  )

# ------------- COSINE SIMILARITY mRNA DATA ------------------------

mRNA_rad_tbl = data.frame(
    feature = loading_split_rad$mRNA$feature,
    F1_FABIA = loading_split_rad$mRNA$loading_FABIA1,
    F1_MOFA = loading_split_rad$mRNA$loading_MOFA1,
    F1_MFA = loading_split_rad$mRNA$loading_MFA1,
    F1_GFA = loading_split_rad$mRNA$loading_GFA1
  )

# flip all to match FABIA direction first:
mRNA_rad_aligned <- align_to_reference(
  mRNA_rad_tbl, ref_col = "F1_FABIA",
  target_cols = c("F1_MOFA","F1_MFA","F1_GFA")
)$data

colnames(mRNA_rad_aligned) <- c("feature", "FABIA", "MOFA", "MFA", "GFA")

# reorder: feature, FABIA, MOFA, GFA, MFA
mRNA_rad_aligned <- mRNA_rad_aligned[
  , c("feature", "FABIA", "MOFA", "GFA", "MFA")
]

# --- cosine similarity ---
cosine_sim <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y)) stop("length mismatch")
  nx <- sqrt(sum(x^2)); ny <- sqrt(sum(y^2))
  if (nx == 0 || ny == 0) stop("cosine undefined for zero vector")
  sum(x * y) / (nx * ny)
}

# --- extract the score matrix (drop sample column) ---
mRNAM <- as.matrix(mRNA_rad_aligned[, c("FABIA","MOFA","GFA","MFA")])

# --- pairwise cosine similarity matrix ---
methods <- colnames(mRNAM)
mRNA_cos_mat <- outer(seq_along(methods), seq_along(methods),
                       Vectorize(function(i, j) cosine_sim(mRNAM[, i], mRNAM[, j])))
dimnames(mRNA_cos_mat) <- list(methods, methods)

mRNA_cos_mat

# keep lower triangle (including diagonal)
mRNA_cos_lower <- mRNA_cos_mat
mRNA_cos_lower[upper.tri(mRNA_cos_lower)] <- NA

# round to two decimals
mRNA_cos_lower <- round(mRNA_cos_lower, 2)

library(ggplot2)

df_lower <- as.data.frame(as.table(mRNA_cos_lower))
colnames(df_lower) <- c("Method1", "Method2", "Cosine")

ggplot(df_lower, aes(Method1, Method2, fill = Cosine)) +
  geom_tile(color = "white", linewidth = 0.3, na.rm = TRUE) +
  geom_text(
    aes(label = ifelse(is.na(Cosine), "", sprintf("%.2f", Cosine))),
    color = "white",
    size = 4
  ) +
  coord_fixed() +
  scale_fill_gradient(
    low = "red",
    high = "#08306B",
    limits = c(0, 1),
    na.value = "white"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  )

