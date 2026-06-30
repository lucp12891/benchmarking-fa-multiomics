#===========================================
# Libraries
#===========================================
library(fabia); library(MOFA2); library(MOFAdata);library(MOFAdata);library(data.table); 
library(gridExtra); library(ggplot2); library(tidyverse); library(dplyr); library(GFA);
library(dplyr); library(reshape2); library(patchwork); library(grid)

# H I P P O C A M P U S     A N A L Y S I S (D A Y 1 7)
#===========================================
# Read datasets: 
#===========================================
setwd("C:/Users/bosangir/OneDrive - Studiecentrum voor Kernenergie/Desktop/Ausan - DELL/Output - Paper II/R project/paper_II/Sep 24, 2025/2. Data")
rad_metadata <- read.csv("rad_metadata.csv", header = TRUE, check.names = FALSE, sep = ";")
mRNA <- read.csv("D17_mRNA_Hippocampus_normalized_tmm.csv", header = TRUE, check.names = FALSE, row.names = 1)
proteomics <- read.csv("D17_Protein_Hippocampus_normalized_quant.csv", header = TRUE, check.names = FALSE, row.names = 1)
rad_metadata <- rad_metadata[-17,]
# drugs
mRNAx <- data.frame(mRNA)
image(c(1:dim(t(mRNA))[1]),c(1:dim(t(mRNA))[2]),t(mRNA), ylab="mRNA",xlab="samples")

# methylation
proteomicsx <- data.frame(proteomics)
image(c(1:dim(t(proteomics))[1]),c(1:dim(t(proteomics))[2]),t(proteomics), ylab="proteomics",xlab="samples")

# Create the list
rad_datax <- list(mRNA = as.matrix(mRNA), proteomics = as.matrix(proteomics))
lapply(rad_datax,dim)
rad_metadata <- rad_metadata

# Create the MOFA object and train the model
MOFAobject <- create_mofa(rad_datax)
MOFAobject

# Plot data overview
plot_data_overview(MOFAobject)

# Define MOFA options
## Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts

## Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 2
model_opts

## Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts

## Train the MOFA model
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)
outfile_object_rad = paste0(getwd(),"model_object_rad.hdf5")
MOFAobject.rad.trained <- run_mofa(MOFAobject, outfile_object_rad)

# Overview of the trained MOFA model
## Slots - The MOFA object consists of multiple slots where relevant data and information is stored. For descriptions, you can read the documentation using ?MOFA. The most important slots are:
# data: input data used to train the model (features are centered at zero mean)
# samples_metadata: sample metadata information
# expectations: expectations of the posterior distributions for the Weights and the Factors
slotNames(MOFAobject)

## Add sample metadata to the model
samples_metadata(MOFAobject) <- rad_metadata
samples_metadata(MOFAobject.rad.trained) <- rad_metadata

## Correlation between factors # A good sanity check is to verify that the Factors are largely uncorrelated. 
plot_factor_cor(MOFAobject.rad.trained)

## Plot factor values
plot_factor(MOFAobject.rad.trained, 
            factors = 2, 
            color_by = "Factor2")

## Plot feature weights
plot_weights(MOFAobject.rad.trained,
             view = "mRNA",
             factor = 1,
             nfeatures = 15,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_weights(MOFAobject.rad.trained,
             view = "proteomics",
             factor = 1,
             nfeatures = 15,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject.rad.trained,
                 view = "proteomics",
                 factor = 1,
                 nfeatures = 15,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

plot_factor(MOFAobject.rad.trained, 
            factors = 1, 
            color_by = "Condition",
            add_violin = TRUE,
            dodge = TRUE
)

plot_factor(MOFAobject.rad.trained, 
            factors = 1, 
            color_by = "Treatment",
            add_violin = TRUE,
            dodge = TRUE
)

## Plot molecular signatures in the input data
plot_data_scatter(MOFAobject.rad.trained, 
                  view = "mRNA",
                  factor = 1,  
                  features = 4,
                  sign = "positive",
                  color_by = "Treatment"
) + labs(y="mRNA")

plot_data_heatmap(MOFAobject.rad.trained, 
                  view = "mRNA",
                  factor = 1,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)
## Plot variance decomposition
### Variance decomposition by Factor
plot_variance_explained(MOFAobject.rad.trained, max_r2=15)

### Total variance explained per view
plot_variance_explained(MOFAobject.rad.trained, plot_total = T)[[2]]
MOFAobject.rad.trained@cache[["variance_explained"]][["r2_per_factor"]][["group1"]]
MOFAobject.rad.trained@cache[["variance_explained"]][["r2_total"]][["group1"]]

# Factor scores
factor_mofa <- get_factors(MOFAobject.rad.trained, factors = 'all', as.data.frame = T)
factor_mofa <- factor_mofa %>%
  spread(factor, value)
factor_mofa2 = factor_mofa[, !(colnames(factor_mofa) %in% c('group'))]
colnames(factor_mofa2) =  c('samples', 'score_MOFA1', 'score_MOFA2')
#factor_mofa2 = factor_mofa2[, -3]

# Weights/Loading 
loading_mofa<- get_weights(MOFAobject.rad.trained, factors = 'all', as.data.frame = T)
loading_mofa2 <- loading_mofa %>%
  spread(factor, value)
colnames(loading_mofa2) =  c('feature', 'view', 'loading_MOFA1', 'loading_MOFA2')
#loading_mofa2 = loading_mofa[, -4]

# FABIA

# merge the two datasets
rad_fabia <- rbind(mRNA, proteomics)

# FABIA MODEL
set.seed(123)
FABIAobject.rad.trained <- fabia(rad_fabia,
                         p = 2,           # number of hidden factors = number of biclusters, default = 5
                         alpha = 0.01,     # sparseness loadings (0 - 1.0); default = 0.1
                         cyc = 1000,      # number of iterations; default = 500
                         spl = 0.5,       # sparseness prior loadings (0 - 2.0); default = 0 (Laplace)
                         spz = 0.5,       # sparseness factors (0.5 - 2.0); default = 0.5 (Laplace)
                         random = 1.0,    # random initialization of loadings in [-random,random]; default = 1.0.
                         center = 2,      # data centering: 1 (mean), 2 (median), > 2 (mode), 0 (no); default = 2
                         norm = 2,        # data normalization: 1 (0.75-0.25 quantile), >1 (var=1), 0 (no); default = 1
                         lap = 1.0,       # minimal value of the variational parameter, default = 1
                         nL = 1           # maximal number of biclusters at which a row element can participate; default = 0 (no limit)
)

fabia.rad.scaledData <- FABIAobject.rad.trained@X; 
fabia.rad.scaledData <- as.matrix(fabia.rad.scaledData)

# samples
samples_info.rad = data.frame(rownames(t(fabia.rad.scaledData))); colnames(samples_info.rad) <- c("samples")

# features
features_info.rad = data.frame(rownames(fabia.rad.scaledData)) 
colnames(features_info.rad) <- c("features") 
features_info.rad$ID <- as.integer(row.names(features_info.rad))
features_info.rad$dataview_info <- NA  # Create a new column with NAs
features_info.rad$dataview_info[1:18329] <- 'mRNA'  # Assign 'x' to rows 1 to 50
features_info.rad$dataview_info[18330:21340] <- 'proteomics'  # Assign 'y' to rows 51 to 100

# fabia scores
data_factors_fabia.rad = data.frame(t(FABIAobject.rad.trained@Z)) 
data_factors_fabia.rad$samples = rownames(data_factors_fabia.rad); data_factors_fabia.rad$samples <-as.factor(data_factors_fabia.rad$samples)
colnames(data_factors_fabia.rad) <- c("score_FABIA1", "score_FABIA2", "samples")
colnames(rad_metadata)[[1]] <- c("sample")
data_factors_fabia.rad <- merge(x=data_factors_fabia.rad, y=rad_metadata, by="samples", all = TRUE)
factor_fabia2 <- data_factors_fabia.rad

# Weights/Loading 
loading_fabia = data.frame(FABIAobject.rad.trained@L)
loading_fabia$feature = rownames(loading_fabia)
colnames(loading_fabia) <- c("loading_FABIA1", "loading_FABIA2", "feature")
loading_fabia2 = loading_fabia

# MFA
set.seed(123)
rad_mfa <- rbind(mRNA, proteomics) 
rad_mfa <- as.data.frame(rad_mfa)
#rad_mfa$feature <- rownames(as.data.frame(rad_mfa))
rad_mfa_tr <- data.frame(t(rad_mfa))
MFAobject.rad.trained <- FactoMineR::MFA(rad_mfa_tr, group = c(18329, 3011), type = c("s", "s"),
                                         name.group = c("mRNA", "proteomics"),
                                         graph = FALSE)

# Factor scores
factor_mfa = data.frame(MFAobject.rad.trained$ind$coord)
colnames(factor_mfa) <- c("score_MFA1", "score_MFA2", "score_MFA3", "score_MFA4", "score_MFA5")
factor_mfa$samples = rownames(factor_mfa)
factor_mfa2 = factor_mfa[,-c(3:5)]

# Weights/Loading 
loading_mfa = data.frame(MFAobject.rad.trained$quanti.var$coord)
colnames(loading_mfa) <- c("loading_MFA1", "loading_MFA2", "loading_MFA3", "loading_MFA4", "loading_MFA5")
loading_mfa$feature = rownames(loading_mfa)
loading_mfa2 = loading_mfa[, -c(3:5)]

# GFA
set.seed(123)
GFA_rad.data = t(rad_datax) ; str(GFA_rad.data)
norm <- normalizeData(GFA_rad.data, type="center")
train_t <- lapply(norm$train, t)   # now 16 x 18329 and 16 x 3011
lapply(train_t, dim)
model_option <- getDefaultOpts()
model_option$iter.max <- 1000
model_option$iter.burnin <- 500
GFAobject.rad.trained <- gfa(cbind(train_t), K = 2, opts = model_option)

# Factor scores
factor_gfa = as.data.frame(GFAobject.rad.trained$X)
colnames(factor_gfa) <- c("score_GFA1", "score_GFA2")
factor_gfa$samples = rownames(factor_gfa)

# Weights/Loading 
loading_gfa = as.data.frame(GFAobject.rad.trained$W)
colnames(loading_gfa) <- c("loading_GFA1","loading_GFA2")
loading_gfa$feature = rownames(loading_gfa)

# Merging Factor scores data for FABIA, MOFA, MFA, GFA
merged_fabia_mofa <- merge(factor_fabia2, factor_mofa2, by = 'samples', all = TRUE) # Merge factor_fabia2 and factor_mofa2
merged_fabia_mofa_mfa <- merge(merged_fabia_mofa, factor_mfa2, by = 'samples', all = TRUE) # Merge merged_fabia_mofa with factor_mfa2
data_scores <- merge(merged_fabia_mofa_mfa, factor_gfa, by = 'samples', all = TRUE) # Merge merged_fabia_mofa_mfa with factor_gfa

data_scores <- merge(merge(merge(factor_fabia2, factor_mofa2, by = 'samples', all = TRUE), factor_mfa2, by = 'samples', all = TRUE), factor_gfa, by = 'samples', all = TRUE) # Merge merged_fabia_mofa_mfa with factor_gfa

score_cols <- grep("^score", colnames(data_scores), value = TRUE)
score_cols

rescale_m1_1 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))  # constant vector case
  -1 + (x - rng[1]) / diff(rng) * 2
}

data_scores[, score_cols] <- lapply(
  data_scores[, score_cols, drop = FALSE],
  rescale_m1_1
)

# Merging Weights/loading data for FABIA, MOFA, MFA, GFA
merged_fabia_mofa_loading <- merge(loading_fabia2, loading_mofa2, by = 'feature', all = TRUE) # Merge loading_fabia2 and loading_mofa2
merged_fabia_mofa_mfa_loading <- merge(merged_fabia_mofa_loading, loading_mfa2, by = 'feature', all = TRUE) # Merge merged_fabia_mofa_loading with loading_mfa2
data_weights <- merge(merged_fabia_mofa_mfa_loading, loading_gfa, by = 'feature', all = TRUE) # Merge merged_fabia_mofa_mfa with loading_gfa

data_weights <- merge(merge(merge(loading_fabia2, loading_mofa2, by = 'feature', all = TRUE), loading_mfa2, by = 'feature', all = TRUE), loading_gfa, by = 'feature', all = TRUE) # Merge merged_fabia_mofa_mfa with factor_gfa


# Visuaization of scores using correlation matrix plot
library(GGally); library(psych)

#data_scores_cor <- data_scores
#ggpairs(data_scores_cor[c(2,11,12,13)])

#============================================
# Plot based on metadata
#============================================
#============================================
# Updated helpers for wide data_scores_cor
#============================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)   # unit()

# ---- keep your build_palette() as you had it ----
build_palette <- function(pal = NULL) {
  default <- c("FALSE"="#E37449", "TRUE"="#00366C")
  if (is.null(pal)) return(default)
  if (is.list(pal)) pal <- unlist(pal, use.names = TRUE)
  pal <- as.character(pal)
  if (is.null(names(pal))) stop("Custom palette must be a *named* vector with names 'TRUE','FALSE',...")
  default[names(pal)] <- pal
  default
}

# ---- keep your plot_scores_panel() as you had it ----
plot_scores_panel <- function(
    df, y_lab, panel_tag = NULL,
    tag_position = c("tl_in","none"),
    palette = NULL
) {
  tag_position <- match.arg(tag_position)
  pal <- build_palette(palette)
  df$.group <- factor(df$.group, levels = c("TRUE","FALSE"))
  
  p <- ggplot(df, aes(x = .sample_idx, y = .score, fill = .group, color = .group)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(shape = 21, size = 2.7, stroke = 0.8, alpha = 0.95) +
    scale_fill_manual(values = pal, drop = FALSE, name = "Treatment") +
    scale_color_manual(values = pal, drop = FALSE, guide = "none") +
    labs(x = "Samples", y = y_lab) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.x = element_text(size = 14, face = "plain", colour = "black"),
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

# ---------------- utilities for wide df ----------------
resolve_score_col_df <- function(df, method, factor_idx = 1) {
  # Expect names like score_FABIA1, score_MOFA2, ...
  meth <- toupper(method)
  cand <- paste0("score_", meth, factor_idx)
  if (cand %in% names(df)) return(cand)
  # be forgiving about case
  hit <- grep(paste0("^score_?", meth, factor_idx, "$"), names(df), ignore.case = TRUE, value = TRUE)
  if (length(hit)) return(hit[1])
  stop("Could not find column for ", meth, " factor ", factor_idx,
       ". Looked for ", cand, ".")
}

get_sample_col <- function(df) {
  for (s in c("samples","sample","sample.1")) if (s %in% names(df)) return(s)
  NULL
}

scores_df_for_method_df <- function(df, method, factor_idx = 1, group_col = "Treatment") {
  # drop unnamed (NA) columns defensively
  df <- df[, !is.na(names(df)), drop = FALSE]
  
  score_col <- resolve_score_col_df(df, method, factor_idx)
  samp_col  <- get_sample_col(df)
  
  grp_raw <- if (group_col %in% names(df)) df[[group_col]] else NA
  grp_chr <- dplyr::case_when(
    is.na(grp_raw) ~ "NA",
    is.logical(grp_raw) ~ ifelse(grp_raw, "TRUE", "FALSE"),
    as.character(grp_raw) %in% c("TRUE","FALSE","True","False") ~ toupper(as.character(grp_raw)),
    TRUE ~ "NA"
  )
  
  tibble(
    sample      = if (!is.null(samp_col)) df[[samp_col]] else seq_len(nrow(df)),
    .sample_idx = seq_len(nrow(df)),
    .score      = as.numeric(df[[score_col]]),
    .group      = grp_chr,
    method      = toupper(method)
  )
}

# ---------------- 2x2 grid for MOFA/FABIA/MFA/GFA ----------------
plot_scores_scatter_grid_df <- function(
    df,
    factor_idx      = 1,
    group_col       = "Treatment",
    methods         = c("MOFA","FABIA","MFA","GFA"),
    panel_tags      = NULL,                 # e.g., c("A","B","C","D") or NULL
    tag_position    = c("tl_in","none"),
    title_text      = "Factor scores (Treatment TRUE/FALSE)",
    palette         = NULL,                 # passed to plot_scores_panel()
    legend_position = "right",              # "right" or "bottom"
    title_height    = 0.12
) {
  stopifnot(length(methods) == 4)
  tag_position <- match.arg(tag_position)
  
  # build panel data frames
  dfs <- lapply(methods, function(m) scores_df_for_method_df(df, m, factor_idx, group_col))
  names(dfs) <- methods
  
  # normalize panel tags
  if (is.null(panel_tags)) tags <- rep("", 4) else tags <- rep_len(panel_tags, 4)
  
  # build panels
  p_list <- lapply(seq_along(methods), function(i) {
    plot_scores_panel(
      dfs[[methods[i]]],
      paste0("Score (", toupper(methods[i]), ")"),
      panel_tag    = if (tag_position == "tl_in") tags[i] else NULL,
      tag_position = if (tag_position == "tl_in") "tl_in" else "none",
      palette      = palette
    )
  })
  # legend settings on each panel
  leg_theme <- theme(legend.position = legend_position,
                     legend.title = element_text(face = "bold"))
  p_list <- lapply(p_list, `+`, leg_theme)
  
  # title
  title_grob <- ggplot() +
    annotate("text", x = 0, y = 0, label = title_text,
             fontface = 2, size = 5.5, colour = "#6a1b9a") +
    theme_void() + theme(plot.margin = margin(0,0,6,0))
  
  # 2x2 grid
  panel_grid <- (p_list[[1]] | p_list[[2]]) / (p_list[[3]] | p_list[[4]])
  
  title_grob / panel_grid +
    plot_layout(heights = c(title_height, 1), guides = "collect")
}

#============================================
# Example use (Treatment has levels TRUE/FALSE)
#============================================
p_trt <- plot_scores_scatter_grid_df(
  df          = data_scores,
  factor_idx  = 1,                        # or 2, as available
  group_col   = "Treatment",
  methods     = c("MOFA","FABIA","MFA","GFA"),
  panel_tags  = c("A","B","C","D"),
  tag_position= "none"
  )
print(p_trt)

## Define new variable names
#new_var_names <- c("Score (FABIA)", "Score (MOFA)", "Score (MFA)", "Score (GFA)")
#colnames(data_scores_cor)[c(2,11,12,13)] <- new_var_names
#pairs.panels(data_scores_cor[c(2,11,12,13)], hist.col="lightgoldenrod", breaks = 16, pch = 8, cex = 0.5, ellipses=TRUE, lm=TRUE)

# ------------------------------------- VISUALIZATION OF SCORES -------------------------------------------

# scaling
scale_matrix <- function(x) {
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)
  return(2 * (x - min_x) / (max_x - min_x) - 1)
}

scale_data_scores_cor <- data_scores_cor#[c(1,2,9,11,12,13)]
scale_data_scores_cor$`Score (FABIA)` = scale_matrix(scale_data_scores_cor$score_FABIA1)
scale_data_scores_cor$`Score (MOFA)` = -1*(scale_matrix(scale_data_scores_cor$score_MOFA1))
scale_data_scores_cor$`Score (MFA)` = scale_matrix(scale_data_scores_cor$score_MFA1)
scale_data_scores_cor$`Score (GFA)` = scale_matrix(scale_data_scores_cor$score_GFA1)


#============================================
# Example use (Treatment has levels TRUE/FALSE)
#============================================
# Ensure the scaled 'Score (METHOD)' columns exist (skip if already present)
if (!exists("scale_matrix")) scale_matrix <- function(x) as.numeric(scale(x))

need <- c("Score (MOFA)", "Score (FABIA)", "Score (MFA)", "Score (GFA)")
if (!all(need %in% names(scale_data_scores_cor))) {
  if ("score_MOFA1"  %in% names(scale_data_scores_cor)) scale_data_scores_cor$`Score (MOFA)`  <- scale_matrix(scale_data_scores_cor$score_MOFA1)
  if ("score_FABIA1" %in% names(scale_data_scores_cor)) scale_data_scores_cor$`Score (FABIA)` <- scale_matrix(scale_data_scores_cor$score_FABIA1)
  if ("score_MFA1"   %in% names(scale_data_scores_cor)) scale_data_scores_cor$`Score (MFA)`   <- scale_matrix(scale_data_scores_cor$score_MFA1)
  if ("score_GFA1"   %in% names(scale_data_scores_cor)) scale_data_scores_cor$`Score (GFA)`   <- scale_matrix(scale_data_scores_cor$score_GFA1)
}

# Plot (uses your plot_scores_scatter_grid_df + resolve_score_col_df that targets "Score (METHOD)")
p_trt <- plot_scores_scatter_grid_df(
  df           = scale_data_scores_cor,
  factor_idx   = 1,                  # kept for signature; resolver uses "Score (METHOD)"
  group_col    = "Treatment",        # TRUE/FALSE
  methods      = c("MOFA","FABIA","MFA","GFA"),
  panel_tags   = c("A","B","C","D"),
  tag_position = "none",
  title_text   = "" #Scaled factor-1 scores by Treatment
)
print(p_trt)

library(ggrepel)
library(ggplot2)

score_rad_tbl = data.frame(
  sample = scale_data_scores_cor$samples,
  F1_FABIA = scale_data_scores_cor$`Score (FABIA)`,
  F1_MOFA = scale_data_scores_cor$`Score (MOFA)`,
  F1_MFA = scale_data_scores_cor$`Score (MFA)`,
  F1_GFA = scale_data_scores_cor$`Score (GFA)`
)

# ===== pairs: diag = method + hist+density, lower = scatter+lm, upper = |r| =====
plot_scores_pairs_hist_scatter <- function(
    df,
    method_order = c("FABIA","MOFA","GFA","MFA"),  # order on diagonal
    factor_prefix = "F1_",                        # columns like F1_FABIA, ...
    abs_corr = TRUE,
    diag_col = "#00366C",                         # method label color
    hist_col = "grey85",
    dens_col = "black",
    line_col = "red",
    pch_pt = 5, cex_pt = 0.6
) {
  # columns to use (in that order)
  cols <- paste0(factor_prefix, method_order)
  cols <- cols[cols %in% names(df)]
  if (length(cols) < 2L) stop("Need at least two score columns. Got: ", paste(cols, collapse=", "))
  X <- as.data.frame(df[cols], check.names = FALSE)
  
  # method labels (exactly as requested, no F1_ prefix anywhere)
  method_labels <- method_order[method_order %in% sub(paste0("^", factor_prefix), "", cols)]
  
  # --- upper triangle: correlation number
  panel_cor_num <- function(x, y, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs"))
    if (isTRUE(abs_corr)) r <- abs(r)
    text(0.5, 0.5, if (is.finite(r)) sprintf("%.2f", r) else "NA", cex = 2.2, font = 2)
  }
  
  # --- lower triangle: scatter + straight (lm) line
  panel_scatter_lm <- function(x, y, ...) {
    points(x, y, pch = pch_pt, cex = cex_pt)
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= 3) abline(lm(y ~ x), col = line_col, lwd = 1)
  }
  
  # --- diagonal: histogram + density + (only) method name
  diag_hist_label <- local({
    i <- 0L
    function(x, ...) {
      i <<- i + 1L
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(usr[1:2], 0, 1.5))
      h <- hist(x, plot = FALSE)
      y <- if (max(h$counts, na.rm = TRUE) == 0) h$counts else h$counts / max(h$counts, na.rm = TRUE)
      rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], y, col = hist_col, border = "white")
      dx <- try(density(x, na.rm = TRUE), silent = TRUE)
      if (!inherits(dx, "try-error") && is.finite(max(dx$y))) {
        lines(dx$x, dx$y / max(dx$y, na.rm = TRUE), lwd = 1, col = dens_col)
      }
      rx <- range(h$breaks, na.rm = TRUE)
      text(rx[1] + 0.05 * diff(rx), 1.32, method_labels[i], cex = 1.7, font = 2, col = diag_col, adj = c(-0.15, 0.65))
    }
  })
  
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(oma = c(0, 0, 0, 0), mar = c(2.5, 2.5, 1.2, 1), xaxs = "i", yaxs = "i")
  pairs(
    X,
    labels      = rep("", ncol(X)),  # <- suppress default "F1_METHOD" text
    diag.panel  = diag_hist_label,
    upper.panel = panel_cor_num,
    lower.panel = panel_scatter_lm
  )
}

# --- run it on your table -----------------------------------------------------
# score_rad_tbl columns: sample, F1_FABIA, F1_MOFA, F1_MFA, F1_GFA
plot_scores_pairs_hist_scatter(
  score_rad_tbl,
  method_order = c("FABIA","MOFA","GFA","MFA"),
  factor_prefix = "F1_",
  abs_corr = TRUE
)
### PLOT LOADING
data_weights_nonull <- data_weights
data_weights_nonull[is.na(data_weights_nonull)] <- 0

# Split by view: "mRNA" vs "Proteomics"
split_by_view <- function(
    df,
    view_col     = "view",
    mrna_pat     = "^mRNA",
    prot_pat     = "^(Prot|Proteomics)",
    ignore_case  = TRUE
) {
  v <- as.character(df[[view_col]])
  is_mrna <- grepl(mrna_pat, v, ignore.case = ignore_case)
  is_prot <- grepl(prot_pat, v, ignore.case = ignore_case)
  
  list(
    mRNA       = df[is_mrna, , drop = FALSE],
    Proteomics = df[is_prot, , drop = FALSE],
    Other      = df[!(is_mrna | is_prot) & !is.na(v), , drop = FALSE]  # for any leftovers
  )
}

# Use it on your NA->0 data
weights_rad_split <- split_by_view(data_weights_nonull)

# Get each piece
df_mRNA_weights  <- weights_rad_split$mRNA
df_Proteomics    <- weights_rad_split$Proteomics

library(dplyr)

# helper: keep id cols + all columns like "loading_*1"
keep_first_loading <- function(df, id_cols = c("feature", "view")) {
  keep <- c(intersect(id_cols, names(df)),
            grep("^loading_.*1$", names(df), value = TRUE))
  df %>% select(all_of(keep)) %>% as_tibble()
}

# apply to each split
df_mRNA_weights_tbl     <- keep_first_loading(df_mRNA_weights)
df_Proteomics_tbl       <- keep_first_loading(df_Proteomics)

# quick checks
names(df_mRNA_weights_tbl)
names(df_Proteomics_tbl)

# ===== pairs: diag = method + hist+density, lower = scatter+lm, upper = |r| =====
plot_weights_pairs_hist_scatter <- function(
    df,
    method_order = c("FABIA","MOFA","GFA","MFA"),
    prefix       = "loading_",
    suffix_pat   = "1$",                 # looks for loading_METHOD1 by default
    abs_corr     = TRUE,
    diag_col     = "#00366C",
    hist_col     = "grey85",
    dens_col     = "black",
    line_col     = "red",
    pch_pt       = 5, cex_pt = 0.6,
    diag_cex     = 1.7,                  # method label size (diagonal)
    corr_cex     = 2.0                   # correlation font size (upper)
) {
  # find columns per method, prefer prefix+method+suffix_pat, else exact prefix+method
  find_cols <- function(nms, methods, prefix, suffix_pat) {
    out <- character()
    for (m in methods) {
      rx <- paste0("^", prefix, m, suffix_pat)
      hit <- grep(rx, nms, value = TRUE)
      if (length(hit) == 0) {
        alt <- paste0(prefix, m)
        if (alt %in% nms) hit <- alt
      }
      if (length(hit)) out <- c(out, hit)
    }
    out
  }
  
  cols <- find_cols(names(df), method_order, prefix, suffix_pat)
  if (length(cols) < 2L) stop("Need ≥2 loading columns; found: ", paste(cols, collapse = ", "))
  X <- as.data.frame(df[cols], check.names = FALSE)
  # ensure numeric
  for (j in seq_along(X)) if (!is.numeric(X[[j]])) X[[j]] <- suppressWarnings(as.numeric(X[[j]]))
  
  # labels for diagonal (just method names, no prefixes)
  method_labels <- sub(paste0("^", prefix), "", cols)
  method_labels <- sub(suffix_pat, "", method_labels)
  
  # upper: correlation number
  panel_cor_num <- function(x, y, ...) {
    usr <- par("usr"); on.exit(par(usr)); par(usr = c(0,1,0,1))
    r <- suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs"))
    if (isTRUE(abs_corr)) r <- abs(r)
    text(0.5, 0.5, if (is.finite(r)) sprintf("%.2f", r) else "NA", cex = corr_cex, font = 2)
  }
  
  # lower: scatter + straight lm line
  panel_scatter_lm <- function(x, y, ...) {
    points(x, y, pch = pch_pt, cex = cex_pt)
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= 3) abline(lm(y ~ x), col = line_col, lwd = 1)
  }
  
  # diagonal: histogram + density + method label (no F1_, no suffix)
  diag_hist_label <- local({
    i <- 0L
    function(x, ...) {
      i <<- i + 1L
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(usr[1:2], 0, 1.5))
      h <- hist(x, plot = FALSE)
      y <- if (max(h$counts, na.rm = TRUE) == 0) h$counts else h$counts / max(h$counts, na.rm = TRUE)
      rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], y, col = hist_col, border = "white")
      dx <- try(density(x, na.rm = TRUE), silent = TRUE)
      if (!inherits(dx, "try-error") && is.finite(max(dx$y))) {
        lines(dx$x, dx$y / max(dx$y, na.rm = TRUE), lwd = 1, col = dens_col)
      }
      rx <- range(h$breaks, na.rm = TRUE)
      text(rx[1] + 0.05 * diff(rx), 1.32, method_labels[i],
           cex = diag_cex, font = 2, col = diag_col, adj = c(-0.15, 0.65))
    }
  })
  
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(oma = c(0,0,0,0), mar = c(2.5,2.5,1.2,1), xaxs = "i", yaxs = "i")
  pairs(
    X,
    labels      = rep("", ncol(X)),  # suppress column names
    diag.panel  = diag_hist_label,
    upper.panel = panel_cor_num,
    lower.panel = panel_scatter_lm
  )
}

# --- Run it on your two tables (they have loading_*1 columns) ---
plot_weights_pairs_hist_scatter(
  df            = df_mRNA_weights_tbl,
  method_order  = c("FABIA","MOFA","GFA","MFA"),
  prefix        = "loading_",
  suffix_pat    = "1$"   # uses loading_METHOD1
)
scale_df_mRNA_weights_tbl = df_mRNA_weights_tbl
scale_df_mRNA_weights_tbl$loading_FABIA1 = scale_matrix(scale_df_mRNA_weights_tbl$loading_FABIA1)
scale_df_mRNA_weights_tbl$loading_MOFA1 = (-1)*scale_matrix(scale_df_mRNA_weights_tbl$loading_MOFA1)
scale_df_mRNA_weights_tbl$loading_MFA1 = scale_matrix(scale_df_mRNA_weights_tbl$loading_MFA1)
scale_df_mRNA_weights_tbl$loading_GFA1 = scale_matrix(scale_df_mRNA_weights_tbl$loading_GFA1)
plot_weights_pairs_hist_scatter(
  df            = scale_df_mRNA_weights_tbl,
  method_order  = c("FABIA","MOFA","GFA","MFA"),
  prefix        = "loading_",
  suffix_pat    = "1$"   # uses loading_METHOD1
)

scale_df_Proteomics_tbl = df_Proteomics_tbl
scale_df_Proteomics_tbl$loading_FABIA1 = scale_matrix(scale_df_Proteomics_tbl$loading_FABIA1)
scale_df_Proteomics_tbl$loading_MOFA1 = (-1)*scale_matrix(scale_df_Proteomics_tbl$loading_MOFA1)
scale_df_Proteomics_tbl$loading_MFA1 = scale_matrix(scale_df_Proteomics_tbl$loading_MFA1)
scale_df_Proteomics_tbl$loading_GFA1 = scale_matrix(scale_df_Proteomics_tbl$loading_GFA1)
plot_weights_pairs_hist_scatter(
  df            = scale_df_Proteomics_tbl,
  method_order  = c("FABIA","MOFA","GFA","MFA"),
  prefix        = "loading_",
  suffix_pat    = "1$"   # uses loading_METHOD1
)

# Non-zero FABIA
library(dplyr)

df_mRNA_weights_tbl_nz  <- scale_df_mRNA_weights_tbl  %>% filter(loading_FABIA1 != 0)
df_Proteomics_tbl_nz  <- scale_df_Proteomics_tbl  %>% filter(loading_FABIA1 != 0)

plot_weights_pairs_hist_scatter(
  df            = df_Proteomics_tbl,
  method_order  = c("FABIA","MOFA","GFA","MFA"),
  prefix        = "loading_",
  suffix_pat    = "1$"
)

varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.60)[[1]]
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}

# ------------------------------ JACCARD INDEX ------------------------------------
# Jaccard index
jaccard.similarity.index <- function(data){
  # Extract index columns
  index_columns <- grep("_index$", names(data), value = TRUE)
  # Initialize an empty dataframe with a single row of NA values
  jaccard_df <- data.frame(
    Column1 = NA_character_,
    Column2 = NA_character_,
    similar_count = NA_real_,
    dissimilar_count = NA_real_,
    sum_column1 = NA_real_,
    sum_column2 = NA_real_,
    jaccard_index = NA_real_,
    stringsAsFactors = FALSE
  )
  jaccard_df <- jaccard_df[0, ]  # Remove the initial row
  # Calculate Jaccard similarity for each pair of index columns
  for (i in 1:(length(index_columns) - 1)) {
    for (j in (i + 1):length(index_columns)) {
      column1 <- index_columns[i]
      column2 <- index_columns[j]
      if (length(data[[column1]]) != length(data[[column2]])) {
        stop("Columns must have the same length.")
      }
      # Calculate Jaccard similarity index
      similar_count <- sum(data[[column1]] == 1 & data[[column2]] == 1)
      dissimilar_count <- sum(data[[column1]] != data[[column2]])
      # Modify this part to correctly calculate 
      sum_column1 <- sum(data[[column1]] == 1)
      sum_column2 <- sum(data[[column2]] == 1)
      #total_sum
      jaccard_index <- similar_count / (sum_column1 + sum_column2 - similar_count)
      ji = ifelse(jaccard_index == -1, 1, jaccard_index)
      # Append the results to the dataframe
      ji <- round(ji, 2)
      jaccard_df <- rbind(jaccard_df, c(column1, column2, ji, similar_count, dissimilar_count, sum_column1, sum_column2))
    }
  }
  # Rename the columns in the final result
  colnames(jaccard_df) <- c("Column1", "Column2", "ji", "similar_count", "dissimilar_count", "sum_column1", "sum_column2")
  return(jaccard_df)
}

# Install necessary package
if (!require(ggVennDiagram)) install.packages("ggVennDiagram")

# Load required libraries
library(ggVennDiagram)

# ---------------------------------- Overall ------------------------------------

# Define your list of binary vectors
x = list(
  A = data_scores_cor_index$FABIA_score_index,  # FABIA set
  B = data_scores_cor_index$MOFA_score_index, # MOFA set
  C = data_scores_cor_index$MFA_score_index,  # MFA set,
  D = data_scores_cor_index$GFA_score_index  # GFA set
)

# Convert the binary vectors into positions for the Venn diagram
x <- list(
  A = which(x$A == 1), 
  B = which(x$B == 1),
  C = which(x$C == 1),
  D = which(x$D == 1)
)

# Generate the 4D Venn diagram
ggVennDiagram(x, relative_width = 0.1) + 
  scale_fill_gradient(low = "#BBBFBF", high = "#41709C")

# --------------------------------------- FEATURES ------------------------------------------

# ------------------ JACCARD INDEX ---------------------------
# Jaccard index for weights
ji_data_weights_cor_mRNA <- df_mRNA_weights_tbl
data_weights_mRNA_index <- calculate_threshold_and_index(ji_data_weights_cor_mRNA)
weights_index_columns_methyl.rad <- grep("_index$", names(data_weights_cor_methyl_index), value = TRUE)

ji_data_weights_cor_drugs <- data_weights.drugs
data_weights_cor_drugs_index <- calculate_threshold_and_index(ji_data_weights_cor_drugs)
weights_index_columns_drugs.rad <- grep("_index$", names(data_weights_cor_drugs_index), value = TRUE)

jaccard_weights_data_methyl.rad = jaccard.similarity.index(data_weights_cor_methyl_index)
jaccard_weights_data_drugs.rad = jaccard.similarity.index(data_weights_cor_drugs_index)

#
# Install necessary package
if (!require(ggVennDiagram)) install.packages("ggVennDiagram")

# Load required libraries
library(ggVennDiagram)


# Define your list of binary vectors
x = list(
  A = data_weights_cor_drugs_index$FABIA_index,  # FABIA set
  B = data_weights_cor_drugs_index$MOFA_index # MOFA set
)

# Convert the binary vectors into positions for the Venn diagram
x <- list(
  A = which(x$A == 1), 
  B = which(x$B == 1)
)

# Generate the 4D Venn diagram
ggVennDiagram(x, relative_width = 0.1) + 
  scale_fill_gradient(low = "#BBBFBF", high = "#41709C")

# ---------------------------------- Overall ------------------------------------

# Define your list of binary vectors
x = list(
  A = data_weights_cor_drugs_index$FABIA_index,  # FABIA set
  B = data_weights_cor_drugs_index$MOFA_index, # MOFA set
  C = data_weights_cor_drugs_index$MFA_index,  # MFA set,
  D = data_weights_cor_drugs_index$GFA_index  # GFA set
)

# Convert the binary vectors into positions for the Venn diagram
x <- list(
  A = which(x$A == 1), 
  B = which(x$B == 1),
  C = which(x$C == 1),
  D = which(x$D == 1)
)

# Generate the 4D Venn diagram
ggVennDiagram(x, relative_width = 0.1) + 
  scale_fill_gradient(low = "#BBBFBF", high = "#41709C")




















# Cosine similarity
#install.packages("lsa")
library(lsa)
cosine(mRNA_data_17$loading_FABIA, mRNA_data_17$loading_GFA)
cosine(mRNA_data_17$loading_FABIA, mRNA_data_17$loading_MFA)
cosine(mRNA_data_17$loading_FABIA, mRNA_data_17$loading_MOFA)
cosine(mRNA_data_17$loading_MFA, mRNA_data_17$loading_MOFA)
cosine(mRNA_data_17$loading_GFA, mRNA_data_17$loading_MOFA)

cosine(Protein_data_17$loading_FABIA, Protein_data_17$loading_GFA)
cosine(Protein_data_17$loading_FABIA, Protein_data_17$loading_MFA)
cosine(Protein_data_17$loading_FABIA, Protein_data_17$loading_MOFA)
cosine(Protein_data_17$loading_MFA, Protein_data_17$loading_MOFA)
cosine(Protein_data_17$loading_GFA, Protein_data_17$loading_MOFA)

#######################################################################################
# Visualization
mRNA_data_17$ID = row_number(mRNA_data_17)
ggplot(mRNA_data_17, aes(x=ID)) + 
  geom_point(aes(y = loading_FABIA), color = "darkred") + 
  geom_point(aes(y = loading_MFA), color="steelblue")  +
  geom_point(aes(y = (-1)*loading_MOFA), color="blue")+
  geom_point(aes(y = loading_GFA), color="yellow3")+
  labs(x = "Features", y = "Loading") +
  theme_bw()

# ------------------------------------------ END VISUALIZATION --------------------------------------------
# Visuaization of weights using correlation matrix plot
library(GGally); library(psych)

data_weights_cor <- data_weights
ggpairs(data_weights_cor[c(2,4,6,7)])


# Visualize scores using scatter plot
# Plot
data_scores_cor$IGHV <- as.factor(data_scores_cor$IGHV)
ggplot(data_scores_cor, aes(x = sample, y = `Score (GFA)`)) +
  geom_point(aes(fill = IGHV), size = 4, shape = 21) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("GFA: F 1") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c("0" = "yellow","1" = "steelblue1", "NA" = "grey")) +  # Manually assign colors to Condition
  labs(x = "Samples", y = "Score")

# Top samples
# Sort the values and get the top 50 values
sorted_values <- sort(abs(data_scores_cor$`Score (FABIA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = 15)

# Identify the corresponding feature names for the top 50 values
top_50_samples <- data_scores_cor$sample[match(top_50_values, abs(data_scores_cor$`Score (FABIA)`))]

print(top_50_samples)

# Visualize weights/loading using scatter plot

# Filter the data for only View=drugs
data_weights_cor_drugs <- data_weights_cor %>%
  filter(view == "Drugs")

# Top drugs
# Sort the values and get the top 50 values
sorted_values <- sort(abs(data_weights_cor_drugs$`Loading (FABIA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = 15)

N =100
# Sort the values and get the top 50 values
sorted_values <- sort(abs(data_weights_cor_drugs$`Loading (FABIA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_FABIA <- data_weights_cor_drugs$feature[match(top_50_values, abs(data_weights_cor_drugs$`Loading (FABIA)`))]

sorted_values <- sort(abs(data_weights_cor_drugs$`Loading (MOFA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_MOFA <- data_weights_cor_drugs$feature[match(top_50_values, abs(data_weights_cor_drugs$`Loading (MOFA)`))]

sorted_values <- sort(abs(data_weights_cor_drugs$`Loading (MFA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_MFA <- data_weights_cor_drugs$feature[match(top_50_values, abs(data_weights_cor_drugs$`Loading (MFA)`))]

sorted_values <- sort(abs(data_weights_cor_drugs$`Loading (GFA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_GFA <- data_weights_cor_drugs$feature[match(top_50_values, abs(data_weights_cor_drugs$`Loading (GFA)`))]

# Intersection
find_intersection <- function(vec1, vec2, vec3, vec4) {
  # Find intersection of vec1 and vec2
  intersection_vec <- intersect(vec1, vec2)
  
  # Find intersection of intersection_vec and vec3
  intersection_vec <- intersect(intersection_vec, vec3)
  
  # Find intersection of intersection_vec and vec4
  intersection_vec <- intersect(intersection_vec, vec4)
  
  return(intersection_vec)
}

result_intersection <- find_intersection(top_50_features_FABIA,top_50_features_MOFA, top_50_features_GFA, top_50_features_MFA)
print(result_intersection)

# Identify the corresponding feature names for the top 50 values
top_50_features <- data_weights_cor_drugs$feature[match(top_50_values, abs(data_weights_cor_drugs$`Loading (FABIA)`))]

print(top_50_features)

# Filter the data for only View=methylation
data_weights_cor_meth <- data_weights_cor %>%
  filter(view == "Methylation")

# Top methylation
N =200
# Sort the values and get the top 50 values
sorted_values <- sort(abs(data_weights_cor_meth$loading_FABIA1), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_FABIA <- data_weights_cor_meth$feature[match(top_50_values, abs(data_weights_cor_meth$loading_FABIA1))]

sorted_values <- sort(abs(data_weights_cor_meth$loading_MOFA1), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_MOFA <- data_weights_cor_meth$feature[match(top_50_values, abs(data_weights_cor_meth$loading_MOFA1))]

sorted_values <- sort(abs(data_weights_cor_meth$loading_MFA1), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_MFA <- data_weights_cor_meth$feature[match(top_50_values, abs(data_weights_cor_meth$loading_MFA1))]

sorted_values <- sort(abs(data_weights_cor_meth$loading_GFA1), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_GFA <- data_weights_cor_meth$feature[match(top_50_values, abs(data_weights_cor_meth$loading_GFA1))]

# Intersection
find_intersection <- function(vec1, vec2, vec3){#, vec4) {
  # Find intersection of vec1 and vec2
  intersection_vec <- intersect(vec1, vec2)
  
  # Find intersection of intersection_vec and vec3
  intersection_vec <- intersect(intersection_vec, vec3)
  
  # Find intersection of intersection_vec and vec4
  #intersection_vec <- intersect(intersection_vec, vec4)
  
  return(intersection_vec)
}

result_intersection <- find_intersection(top_50_features_MOFA, top_50_features_MFA, top_50_features_FABIA)
print(result_intersection)

# Intersection
find_intersection <- function(vec1, vec2, vec3, vec4) {
  # Find intersection of vec1 and vec2
  intersection_vec <- intersect(vec1, vec2)
  
  # Find intersection of intersection_vec and vec3
  intersection_vec <- intersect(intersection_vec, vec3)
  
  # Find intersection of intersection_vec and vec4
  intersection_vec <- intersect(intersection_vec, vec4)
  
  return(intersection_vec)
}

result_intersection <- find_intersection(top_50_features_FABIA, top_50_features_MOFA, top_50_features_MFA, top_50_features_GFA)
print(result_intersection)




# Correlation plots
# Define new variable names
new_var_names <- c("Loading (FABIA)", "Loading (MOFA)", "Loading (MFA)", "Loading (GFA)")
colnames(data_weights_cor_meth)[c(2,4,6,7)] <- new_var_names
pairs.panels(data_weights_cor_meth[c(2,4,6,7)], hist.col="deepskyblue4", breaks = 16, pch = 8, cex = 0.5, ellipses=TRUE, lm=TRUE)

# Plot using a single palette
plot1 = ggplot(data_weights_cor_drugs, aes(x = `Loading (MFA)`, y = `Loading (GFA)`)) + 
  geom_point(color = "dodgerblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Loading (MFA)", y = "Loading (GFA)") +
  ggtitle("Drugs") +
  theme(
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 1))+
  theme_minimal()

# Second plot
plot2 = ggplot(data_weights_cor_meth, aes(x = `Loading (MFA)`, y = `Loading (GFA)`)) + 
  geom_point(color = "deepskyblue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Loading (MFA)", y = "Loading (MFA)") +
  ggtitle("Methylation") +
  theme_minimal()

# Arrange plots side by side
grid.arrange(plot1, plot2, ncol = 2)

# Functions 
# Threshold calculator
calculate_threshold_and_index <- function(data) {
  for (i in seq_along(data)) {
    if (is.numeric(data[[i]])) {
      threshold <- varphi.data(data[[i]])
      
      # Create corresponding index column name
      index_column <- paste0(names(data)[i], "_index")
      
      # Create index column based on the threshold
      data[[index_column]] <- ifelse(data[[i]] > threshold, 1, 0)
    }
  }
  return(data)
}

# Jaccard Similarity Index
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.90)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}

# Function to re-scale a vector to (-1,1) scale
rescale_minusone_to_one <- function(vector){
  scaled_vector <- vector / max(abs(vector), 1e-6) # Using 1e-6 to avoid division by zero
  return(scaled_vector)
}

# subset data to only scores
data_scores_cor_subset <- data_scores_cor[,c(2,11,13,14)]

# Obtain thresholds
data_scores_cor_subset.i <- data_scores_cor_subset %>%
  mutate_if(is.numeric, ~ rescale_minusone_to_one(.))

data_scores_cor_subset.a <- data_scores_cor_subset.i %>%
  mutate_if(is.numeric, ~ abs(.))

varphi.data(abs(data_scores_cor_subset.a$`Score (FABIA)`))
varphi.data(abs(data_scores_cor_subset.a$`Score (MOFA)`))
varphi.data(abs(data_scores_cor_subset.a$`Score (MFA)`))
varphi.data(abs(data_scores_cor_subset.a$`Score (GFA)`))

# calculate_threshold_and_index
data_scores_cor_subset.b <- calculate_threshold_and_index(data_scores_cor_subset.a)
data_scores_cor_index.b <- grep("_index$", names(data_scores_cor_subset.b), value = TRUE)

# Subset the dataframe based on selected columns
data_scores_cor_subset.c <- data_scores_cor_subset.b[, c(data_scores_cor_index.b)]
jaccard_scores = jaccard.index(data_scores_cor_subset.c)

# 2. Loading/Weight
# Jaccard Similarity Index
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.60)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}

# subset data to only loading
data_weights_cor_drugs_subset <- data_weights_cor_drugs[,c(2,4,6,7)]

# Obtain thresholds
data_weights_cor_drugs_subset.i <- data_weights_cor_drugs_subset %>%
  mutate_if(is.numeric, ~ rescale_minusone_to_one(.))

data_weights_cor_drugs_subset.a <- data_weights_cor_drugs_subset.i %>%
  mutate_if(is.numeric, ~ abs(.))

varphi.data(abs(data_weights_cor_drugs_subset.a$`Loading (FABIA)`))
varphi.data(abs(data_weights_cor_drugs_subset.a$`Loading (MOFA)`))
varphi.data(abs(data_weights_cor_drugs_subset.a$`Loading (MFA)`))
varphi.data(abs(data_weights_cor_drugs_subset.a$`Loading (GFA)`))

# calculate_threshold_and_index
data_weights_cor_drugs_subset.b <- calculate_threshold_and_index(data_weights_cor_drugs_subset.a)
data_weights_cor_drugs_index.b <- grep("_index$", names(data_weights_cor_drugs_subset.b), value = TRUE)

# Subset the dataframe based on selected columns
data_weights_cor_drugs_subset.c <- data_weights_cor_drugs_subset.b[, c(data_weights_cor_drugs_index.b)]
jaccard_scores_drugs = jaccard.index(data_weights_cor_drugs_subset.c)

# (ii) Methylation
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.95)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}
# subset data to only loading
data_weights_cor_meth_subset <- data_weights_cor_meth[,c(2,4,6,7)]

# Obtain thresholds
data_weights_cor_meth_subset.i <- data_weights_cor_meth_subset %>%
  mutate_if(is.numeric, ~ rescale_minusone_to_one(.))

data_weights_cor_meth_subset.a <- data_weights_cor_meth_subset.i %>%
  mutate_if(is.numeric, ~ abs(.))

varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (FABIA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (MOFA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (MFA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (GFA)`))

# calculate_threshold_and_index
data_weights_cor_meth_subset.b <- calculate_threshold_and_index(data_weights_cor_meth_subset.a)
data_weights_cor_meth_index.b <- grep("_index$", names(data_weights_cor_meth_subset.b), value = TRUE)

# Subset the dataframe based on selected columns
data_weights_cor_meth_subset.c <- data_weights_cor_meth_subset.b[, c(data_weights_cor_meth_index.b)]
jaccard_scores_meth = jaccard.index(data_weights_cor_meth_subset.c)


# Functions 
# Threshold calculator
calculate_threshold_and_index <- function(data) {
  for (i in seq_along(data)) {
    if (is.numeric(data[[i]])) {
      threshold <- varphi.data(data[[i]])
      
      # Create corresponding index column name
      index_column <- paste0(names(data)[i], "_index")
      
      # Create index column based on the threshold
      data[[index_column]] <- ifelse(data[[i]] > threshold, 1, 0)
    }
  }
  return(data)
}

# Jaccard Similarity Index
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.60)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}

# Function to re-scale a vector to (-1,1) scale
rescale_minusone_to_one <- function(vector){
  scaled_vector <- vector / max(abs(vector), 1e-6) # Using 1e-6 to avoid division by zero
  return(scaled_vector)
}

data_scores_cor_f2_subset.i <- data_scores_cor_F1[,c("score_FABIA1", "score_MOFA2", "score_MFA2", "score_GFA2")]
# subset data to only scores
data_scores_cor_f2_subset.ii <- df_normalized[,c("score_FABIA1", "score_MOFA2", "score_MFA2", "score_GFA2")]

# Obtain thresholds
#data_scores_cor_f1_subset.i <- data_scores_cor_f1_subset %>%
#  mutate_if(is.numeric, ~ rescale_minusone_to_one(.))

data_scores_cor_f2_subset.a <- data_scores_cor_f2_subset.ii %>%
  mutate_if(is.numeric, ~ abs(.))

varphi.data(abs(data_scores_cor_f2_subset.a$score_FABIA1))
varphi.data(abs(data_scores_cor_f2_subset.a$score_MOFA2))
varphi.data(abs(data_scores_cor_f2_subset.a$score_MFA2))
varphi.data(abs(data_scores_cor_f2_subset.a$score_GFA2))


# Jaccard Similarity Index
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.60)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}
# (ii) Methylation
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.95)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}
# subset data to only loading
data_weights_cor_meth_subset <- data_weights_cor_meth[,c(2,4,6,7)]

# Obtain thresholds
data_weights_cor_meth_subset.i <- data_weights_cor_meth_subset %>%
  mutate_if(is.numeric, ~ rescale_minusone_to_one(.))

data_weights_cor_meth_subset.a <- data_weights_cor_meth_subset.i %>%
  mutate_if(is.numeric, ~ abs(.))

varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (FABIA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (MOFA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (MFA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (GFA)`))

# calculate_threshold_and_index
data_weights_cor_meth_subset.b <- calculate_threshold_and_index(data_weights_cor_meth_subset.a)
data_weights_cor_meth_index.b <- grep("_index$", names(data_weights_cor_meth_subset.b), value = TRUE)

# Subset the dataframe based on selected columns
data_weights_cor_meth_subset.c <- data_weights_cor_meth_subset.b[, c(data_weights_cor_meth_index.b)]
jaccard_scores_meth = jaccard.index(data_weights_cor_meth_subset.c)

# ------------------------------------------------ FUNCTIONS -----------------------------------------------------------------

# Jaccard index
jaccard.similarity.index <- function(data){
  # Extract index columns
  index_columns <- grep("_index$", names(data), value = TRUE)
  # Initialize an empty dataframe with a single row of NA values
  jaccard_df <- data.frame(
    Column1 = NA_character_,
    Column2 = NA_character_,
    similar_count = NA_real_,
    dissimilar_count = NA_real_,
    sum_column1 = NA_real_,
    sum_column2 = NA_real_,
    jaccard_index = NA_real_,
    stringsAsFactors = FALSE
  )
  jaccard_df <- jaccard_df[0, ]  # Remove the initial row
  # Calculate Jaccard similarity for each pair of index columns
  for (i in 1:(length(index_columns) - 1)) {
    for (j in (i + 1):length(index_columns)) {
      column1 <- index_columns[i]
      column2 <- index_columns[j]
      if (length(data[[column1]]) != length(data[[column2]])) {
        stop("Columns must have the same length.")
      }
      # Calculate Jaccard similarity index
      similar_count <- sum(data[[column1]] == 1 & data[[column2]] == 1)
      dissimilar_count <- sum(data[[column1]] != data[[column2]])
      sum_column1 <- sum(data[[column1]] == 1)
      sum_column2 <- sum(data[[column2]] == 1)
      jaccard_index <- similar_count / (sum_column1 + sum_column2 - similar_count)
      ji = ifelse(jaccard_index == -1, 1, jaccard_index)
      # Append the results to the dataframe
      ji <- round(ji, 2)
      jaccard_df <- rbind(jaccard_df, c(column1, column2, ji, similar_count, dissimilar_count, sum_column1, sum_column2))
    }
  }
  # Rename the columns in the final result
  colnames(jaccard_df) <- c("Column1", "Column2", "ji", "similar_count", "dissimilar_count", "sum_column1", "sum_column2")
  return(jaccard_df)
}
# ---------------------------------------------- CREATING SCORES DATAFRAME ---------------------------------------------------
# merge with metadata
df_sampleInfo_Hip_D17b$samples = rownames(df_sampleInfo_Hip_D17b)#$Internal.ID
data_score_methods.17 = merge(data_score_df.17, df_sampleInfo_Hip_D17b, by = 'samples', all = TRUE)

# Jaccard index for scores
ji_data_scores_cor <- data_scores_cor[c(1,2,11,12,13)]
new_var_names <- c('sample', 'FABIA_score', 'MOFA_score', 'MFA_score', 'GFA_score')
colnames(ji_data_scores_cor) <- new_var_names
data_scores_cor_index <- calculate_threshold_and_index(ji_data_scores_cor)
score_index_columns.rad <- grep("_index$", names(data_scores_cor_index), value = TRUE)

# Jaccard Index
jaccard_score_data.rad = jaccard.similarity.index(data_scores_cor_index)

# Venn
# install.packages("ggVennDiagram")
library(ggVennDiagram)
x = data_scores_cor_index
# List of items
x <- list(x$FABIA_score_index, x$MOFA_score_index)

# 2D Venn diagram
ggVennDiagram(x) 

################################################################
score_tbl_index = calculate_threshold_and_index(score_tbl)
jac_index = jaccard.similarity.index(score_tbl_index)

