
#====================================================
# Intersection of objects
#====================================================

# Core pkgs
library(ComplexHeatmap)
library(grid)          # for grid.text
library(circlize)      # color utilities (optional)

# lt: named list of character vectors (elements in each set), or logical/binary matrix
# mode: "distinct" (default), "intersect", or "union"
# universal_set: optional vector for adding the complement (e.g., all possible elements)
build_comb_mat <- function(lt, mode = "distinct", universal_set = NULL) {
  make_comb_mat(lt, mode = mode, universal_set = universal_set)
}

plot_upset_pub <- function(m,
                           bar_col = "black",
                           set_bar_col = "black",
                           add_numbers = TRUE,
                           base_size = 14,
                           title = NULL,
                           comb_order = order(comb_degree(m), -comb_size(m)),
                           set_order  = order(set_size(m))) {
  
  cs <- comb_size(m)
  ss <- set_size(m)
  
  UpSet(
    m,
    comb_order = comb_order,
    set_order  = set_order,
    top_annotation = upset_top_annotation(
      m,
      gp = gpar(fill = bar_col, col = NA),
      add_numbers = add_numbers,
      height = unit(28, "mm")
    ),
    right_annotation = upset_right_annotation(
      m,
      gp = gpar(fill = set_bar_col, col = NA),
      width = unit(28, "mm"),
      add_numbers = add_numbers
    ),
    show_row_names = TRUE,
    row_names_side = "left",
    column_title = title
  ) |>
    draw(merge_legend = TRUE)
}

plot_upset_by_degree <- function(m,
                                 pal = c("#2C6BA2", "#2A9D8F", "#D85C4B", "#8E6BBF", "#C49A00"),
                                 base_size = 14,
                                 title = NULL) {
  deg <- comb_degree(m)
  col_by_deg <- setNames(pal[seq_len(max(deg))], seq_len(max(deg)))
  
  UpSet(
    m,
    comb_col = col_by_deg[as.character(deg)],
    top_annotation = upset_top_annotation(
      m,
      gp = gpar(fill = col_by_deg[as.character(deg)], col = NA),
      height = unit(28, "mm"),
      add_numbers = TRUE,
      annotation_name_side = "left",
      annotation_name_rot = 0
    ),
    right_annotation = upset_right_annotation(
      m,
      gp = gpar(fill = "grey20", col = NA),
      width = unit(24, "mm"),
      add_numbers = TRUE
    ),
    column_title = title
  ) |>
    draw(merge_legend = TRUE)
}

plot_upset_fraction <- function(m,
                                ylim = c(0, 1),
                                base_size = 12,
                                title = NULL) {
  frac <- comb_size(m) / sum(comb_size(m))
  
  UpSet(
    m,
    top_annotation = HeatmapAnnotation(
      "Relative\nfraction" = anno_barplot(
        frac,
        ylim = ylim,
        gp = gpar(fill = "black"),
        border = FALSE,
        height = unit(28, "mm")
      ),
      annotation_name_side = "left",
      annotation_name_rot = 0
    ),
    right_annotation = upset_right_annotation(
      m,
      gp = gpar(fill = "grey20", col = NA),
      width = unit(24, "mm"),
      add_numbers = TRUE
    ),
    column_title = title
  ) |>
    draw()
}

# Keep only intersections with size >= min_size
filter_by_size <- function(m, min_size = 1) m[comb_size(m) >= min_size]

# Keep only intersections of a given degree (e.g., 2-way overlaps only)
filter_by_degree <- function(m, degree = 2) m[comb_degree(m) == degree]

# Add complement set by specifying the universe
add_complement <- function(lt, universe) make_comb_mat(lt, universal_set = universe)

# nm must be one of comb_name(m), e.g., "101" or "1101"
get_elements_in_combination <- function(m, nm) {
  extract_comb(m, nm)
}

plot_upset_modes <- function(lt) {
  m1 <- make_comb_mat(lt, mode = "distinct")
  m2 <- make_comb_mat(lt, mode = "intersect")
  m3 <- make_comb_mat(lt, mode = "union")
  
  UpSet(m1, row_title = "distinct") %v%
    UpSet(m2, row_title = "intersect") %v%
    UpSet(m3, row_title = "union") |>
    draw(merge_legend = TRUE)
}

#====================================================
# Intersection of objects: Scores
#====================================================

# x4 should be a named list of index vectors, e.g.:
xscore_rad <- list(
  A = which(score_rad_tbl_index$F1_FABIA_index == 1), 
  B = which(score_rad_tbl_index$F1_GFA_index == 1),
  C = which(score_rad_tbl_index$F1_MOFA_index == 1),
  D = which(score_rad_tbl_index$F1_MFA_index == 1)
)

m <- build_comb_mat(xscore_rad)                    # default mode = "distinct"
m_sz <- filter_by_size(m, min_size = 2)    # optional: keep intersections size >= 2
m_d2 <- filter_by_degree(m, degree = 2)    # optional: only 2-way overlaps

# 1) Publication style (counts + numbers)
plot_upset_pub(m, title = "UpSet — distinct mode")

# 2) Color by degree
plot_upset_by_degree(m, title = "")


library(ggVennDiagram)

pal4 <- c("#2C6BA2", "#2A9D8F", "#D85C4B", "#8E6BBF")  # blue, teal, brick, violet

ggVennDiagram(
  xscore_rad,
  category.names = names(x4),
  show_intersect = FALSE,
  set_color = pal4,            # different outline color per set
  set_size  = 1.2,
  label     = "both",
  label_alpha = 0,
  label_geom  = "label",
  label_color = "black",
  label_size  = 5,
  edge_lty = "solid",
  edge_size = 1,
  force_upset = FALSE,
  order.intersect.by = "size",
  order.set.by = "name",
  relative_height = 3,
  relative_width  = 0.3
) +
  scale_fill_gradient(low = "#EAF2EE", high = "#2A9D8F") +   # soft → deep teal for overlaps
  theme_void() +
  labs(title = "") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    venn.set.label    = element_text(size = 16, face = "bold"),
    venn.region.label = element_text(size = 14)
  )


#====================================================
# Intersection of objects: DNA methylation
#====================================================

# x4 should be a named list of index vectors, e.g.:
# Upset plot - DNA methylation
xmrna <- list(
  A = which(mrna_tbl_index$F1_FABIA_index == 1), 
  B = which(mrna_tbl_index$F1_GFA_index == 1),
  C = which(mrna_tbl_index$F1_MOFA_index == 1),
  D = which(mrna_tbl_index$F1_MFA_index == 1)
)

m <- build_comb_mat(xmrna)                    # default mode = "distinct"
m_sz <- filter_by_size(m, min_size = 2)    # optional: keep intersections size >= 2
m_d2 <- filter_by_degree(m, degree = 2)    # optional: only 2-way overlaps

# 1) Publication style (counts + numbers)
plot_upset_pub(m, title = "UpSet — distinct mode")

# 2) Color by degree
plot_upset_by_degree(m, title = "")

library(ggVennDiagram)

pal4 <- c("#2C6BA2", "#2A9D8F", "#D85C4B", "#8E6BBF")  # blue, teal, brick, violet

ggVennDiagram(
  xmrna,
  category.names = names(x4),
  show_intersect = FALSE,
  set_color = pal4,            # different outline color per set
  set_size  = 1.2,
  label     = "both",
  label_alpha = 0,
  label_geom  = "label",
  label_color = "black",
  label_size  = 5,
  edge_lty = "solid",
  edge_size = 1,
  force_upset = FALSE,
  order.intersect.by = "size",
  order.set.by = "name",
  relative_height = 3,
  relative_width  = 0.3
) +
  scale_fill_gradient(low = "#EAF2EE", high = "#2A9D8F") +   # soft → deep teal for overlaps
  theme_void() +
  labs(title = "") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    venn.set.label    = element_text(size = 16, face = "bold"),
    venn.region.label = element_text(size = 14)
  )

#====================================================
# Intersection of objects: Drugs profiles
#====================================================

# x4 should be a named list of index vectors, e.g.:
xproteins <- list(
  A = which(proteins_tbl_index$F1_FABIA_index == 1), 
  B = which(proteins_tbl_index$F1_GFA_index == 1),
  C = which(proteins_tbl_index$F1_MOFA_index == 1),
  D = which(proteins_tbl_index$F1_MFA_index == 1)
)

m <- build_comb_mat(xdrugs)                    # default mode = "distinct"
m_sz <- filter_by_size(m, min_size = 2)    # optional: keep intersections size >= 2
m_d2 <- filter_by_degree(m, degree = 2)    # optional: only 2-way overlaps

# 1) Publication style (counts + numbers)
plot_upset_pub(m, title = "UpSet — distinct mode")

# 2) Color by degree
plot_upset_by_degree(m, title = "")

library(ggVennDiagram)

pal4 <- c("#2C6BA2", "#2A9D8F", "#D85C4B", "#8E6BBF")  # blue, teal, brick, violet

ggVennDiagram(
  xproteins,
  category.names = names(x4),
  show_intersect = FALSE,
  set_color = pal4,            # different outline color per set
  set_size  = 1.2,
  label     = "both",
  label_alpha = 0,
  label_geom  = "label",
  label_color = "black",
  label_size  = 5,
  edge_lty = "solid",
  edge_size = 1,
  force_upset = FALSE,
  order.intersect.by = "size",
  order.set.by = "name",
  relative_height = 3,
  relative_width  = 0.3
) +
  scale_fill_gradient(low = "#EAF2EE", high = "#2A9D8F") +   # soft → deep teal for overlaps
  theme_void() +
  labs(title = "") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    venn.set.label    = element_text(size = 16, face = "bold"),
    venn.region.label = element_text(size = 14)
  )
