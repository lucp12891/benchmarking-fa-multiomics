  # clean the datasets

  fabload_split_fab<- split_by_pattern(
    df = fabload,
    column = "feature",
    patterns = patterns
  )
  
  fabload_split_fab$drugs$feature  <- gsub("^drugs::",  "", fabload_split_fab$drugs$feature)
  fabload_split_fab$methyl$feature <- gsub("^methyl::", "", fabload_split_fab$methyl$feature)
  

  plot(proteomics_rad_tbl$F1_FABIA)
  plot(mRNA_rad_tbl$F1_FABIA)


  # ---- Packages ----
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(ggrepel)
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  plot_factor_loadings_gg <- function(
    tbl,
    feature_col = "feature",
    loading_col = "F1",
    top_n = 10,
    show_sd = FALSE,        # NEW: optional
    sd_k = 2,
    show_quantile = FALSE, # NEW: optional
    quantile_p = 0.95,
    flip = TRUE,
    title = NULL
  ) {
    stopifnot(is.data.frame(tbl), loading_col %in% names(tbl))
    
    df <- tbl %>%
      { if (!(feature_col %in% names(.)))
        rownames_to_column(., var = feature_col) else . } %>%
      mutate(
        loading = .data[[loading_col]],
        abs_loading = abs(loading)
      )
    
    # Top-N by |loading|
    top_feats <- df %>%
      arrange(desc(abs_loading)) %>%
      slice_head(n = top_n) %>%
      pull(.data[[feature_col]])
    
    df <- df %>%
      mutate(
        is_topN = .data[[feature_col]] %in% top_feats,
        feature_label = ifelse(is_topN, as.character(.data[[feature_col]]), NA_character_),
        !!feature_col := factor(
          .data[[feature_col]],
          levels = .data[[feature_col]][order(abs_loading)]
        )
      )
    
    p <- ggplot(df, aes(x = .data[[feature_col]], y = loading)) +
      geom_point(aes(color = abs_loading), size = 1.8) +
      scale_color_gradient("Loading", low = "steelblue", high = "darkorange") +
      geom_text_repel(
        data = filter(df, is_topN),
        aes(label = feature_label),
        size = 5,
        box.padding = 1.2,
        point.padding = 0.4,
        max.overlaps = Inf
      ) +
      labs(
        title = title %||% sprintf("Factor loadings (%s)", loading_col),
        x = NULL,
        y = sprintf("Loading (Factor-1)", loading_col) #%s
      ) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    
    # ---- OPTIONAL cutoffs ----
    if (show_sd) {
      sd_thr <- sd(df$loading, na.rm = TRUE) * sd_k
      p <- p + geom_hline(
        yintercept = c(-sd_thr, sd_thr),
        linetype = "dashed",
        linewidth = 0.5,
        color = "red"
      )
    }
    
    if (show_quantile) {
      q_thr <- quantile(df$abs_loading, quantile_p, na.rm = TRUE)
      p <- p + geom_hline(
        yintercept = c(-q_thr, q_thr),
        linetype = "dashed",
        linewidth = 0.5,
        color = "blue"
      )
    }
    
    if (flip) p <- p + coord_flip()
    p
  }
  
  # Proteomics – Factor 1
  p_prot_F1 <- plot_factor_loadings_gg(
    proteomics_rad_tbl,
    feature_col = "feature",
    loading_col = "F1_FABIA",
    top_n = 10,
    title = "" #FABIA loadings (Drugs) – Factor 1
  )
  p_prot_F1
  
  # Proteomics – Factor 1
  p_mRNA_F1 <- plot_factor_loadings_gg(
    mRNA_rad_tbl,
    feature_col = "feature",
    loading_col = "F1_FABIA",
    top_n = 10,
    title = "" #FABIA loadings (Drugs) – Factor 1
  )
  p_mRNA_F1
  