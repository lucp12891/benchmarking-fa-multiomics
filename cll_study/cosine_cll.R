# ------------- COSINE SIMILARITY SAMPLE SCORES ------------------------
#--- your objects (assuming they exist exactly as you showed) ---
fabia_sample_scores <- data.frame(sample_scores$FABIA)
mofa_sample_scores  <- data.frame(sample_scores$MOFA)
mfa_sample_scores   <- data.frame(sample_scores$MFA)
gfa_sample_scores   <- data.frame(sample_scores$GFA)

colnames(fabia_sample_scores) <- c("F1 (FABIA)", "F2 (FABIA)")
colnames(mofa_sample_scores)  <- c("F1 (MOFA)",  "F2 (MOFA)")
colnames(mfa_sample_scores)   <- c("F1 (MFA)",   "F2 (MFA)")
colnames(gfa_sample_scores)   <- c("F1 (GFA)",   "F2 (GFA)")

score_all = data.frame(
  sample = rownames(fabia_sample_scores),
  F1_FABIA = fabia_sample_scores$`F1 (FABIA)`,
  F1_MOFA = mofa_sample_scores$`F1 (MOFA)`,
  F1_MFA = mfa_sample_scores$`F1 (MFA)`,
  F1_GFA = gfa_sample_scores$`F1 (GFA)`
)

# flip all to match FABIA direction first:
score_all_aligned <- align_to_reference(
  score_all, ref_col = "F1_FABIA",
  target_cols = c("F1_MOFA","F1_MFA","F1_GFA")
)$data

colnames(score_all_aligned) <- c("sample", "FABIA", "MOFA", "MFA", "GFA")

# reorder: feature, FABIA, MOFA, GFA, MFA
score_all_aligned <- score_all_aligned[
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
M <- as.matrix(score_all_aligned[, c("FABIA","MOFA","GFA","MFA")])

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

# ------------- COSINE SIMILARITY DRUGS DATA ------------------------

drugs_cll_tbl = data.frame(
  feature = fabia_split$drugs$feature,
  F1_FABIA = fabia_split$drugs$`F1 (FABIA)`,
  F1_MOFA = mofa_split$drugs$`F1 (MOFA)`,
  F1_MFA = mfa_split$drugs$`F1 (MFA)`,
  F1_GFA = gfa_split$drugs$`F1 (GFA)`
)

# flip all to match FABIA direction first:
drugs_cll_aligned <- align_to_reference(
  drugs_cll_tbl, ref_col = "F1_FABIA",
  target_cols = c("F1_MOFA","F1_MFA","F1_GFA")
)$data

colnames(drugs_cll_aligned) <- c("feature", "FABIA", "MOFA", "MFA", "GFA")
# reorder: feature, FABIA, MOFA, GFA, MFA
drugs_cll_aligned <- drugs_cll_aligned[
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
drugsM <- as.matrix(drugs_cll_aligned[, c("FABIA","MOFA","GFA","MFA")])

# --- pairwise cosine similarity matrix ---
methods <- colnames(drugsM)
drugs_cos_mat <- outer(seq_along(methods), seq_along(methods),
                 Vectorize(function(i, j) cosine_sim(drugsM[, i], drugsM[, j])))
dimnames(drugs_cos_mat) <- list(methods, methods)

drugs_cos_mat

# keep lower triangle (including diagonal)
drugs_cos_lower <- drugs_cos_mat
drugs_cos_lower[upper.tri(drugs_cos_lower)] <- NA

# round to two decimals
drugs_cos_lower <- round(drugs_cos_lower, 2)

library(ggplot2)

df_lower <- as.data.frame(as.table(drugs_cos_lower))
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

# ------------- COSINE SIMILARITY METHYLATION DATA ------------------------

methyl_cll_tbl = data.frame(
  feature = fabia_split$methyl$feature,
  F1_FABIA = fabia_split$methyl$`F1 (FABIA)`,
  F1_MOFA = mofa_split$methyl$`F1 (MOFA)`,
  F1_MFA = mfa_split$methyl$`F1 (MFA)`,
  F1_GFA = gfa_split$methyl$`F1 (GFA)`
)

# flip all to match FABIA direction first:
methyl_cll_aligned <- align_to_reference(
  methyl_cll_tbl, ref_col = "F1_FABIA",
  target_cols = c("F1_MOFA","F1_MFA","F1_GFA")
)$data

colnames(methyl_cll_aligned) <- c("feature", "FABIA", "MOFA", "MFA", "GFA")

# reorder: feature, FABIA, MOFA, GFA, MFA
methyl_cll_aligned <- methyl_cll_aligned[
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
methylM <- as.matrix(methyl_cll_aligned[, c("FABIA","MOFA","GFA","MFA")])

# --- pairwise cosine similarity matrix ---
methods <- colnames(methylM)
methyl_cos_mat <- outer(seq_along(methods), seq_along(methods),
                       Vectorize(function(i, j) cosine_sim(methylM[, i], methylM[, j])))
dimnames(methyl_cos_mat) <- list(methods, methods)

methyl_cos_mat

# keep lower triangle (including diagonal)
methyl_cos_lower <- methyl_cos_mat
methyl_cos_lower[upper.tri(methyl_cos_lower)] <- NA

# round to two decimals
methyl_cos_lower <- round(methyl_cos_lower, 2)

library(ggplot2)

df_lower <- as.data.frame(as.table(methyl_cos_lower))
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

