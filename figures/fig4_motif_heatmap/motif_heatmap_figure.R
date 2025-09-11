# Library path and quiet loading
.libPaths(c("~/Library/R/4.2/library", .libPaths()))

suppressPackageStartupMessages({
  library(ggraph)
  library(tidygraph)
  library(graphlayouts)
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
  library(cowplot)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(pheatmap)
  library(tidyverse) 
  library(clusterProfiler)
  library(org.Hs.eg.db)
})
# Define color scheme for motif hierarchy 
parent_color <- '#93C5FD'  # Light blue for parent motifs 
child_color <- '#BBD6B8'   # Light green for child motifs 
orphan_color <- '#D1D5DB'  # Light gray for orphan motifs 

# Create the motif hierarchy diagram
create_motif_hierarchy <- function() {
  graph_definition <- "
  digraph motif_hierarchy {
    # Graph settings for publication quality
    graph [rankdir = TB,
           bgcolor = white,
           nodesep = 1.0,
           ranksep = 1.5,
           compound = true,
           splines = ortho,
           pad = 0.2]
    
    # Default node styling 
    node [shape = ellipse,
          style = 'filled,bold',
          color = '#1F2937',
          penwidth = 1.0,
          fontsize = 120,
          fontname = 'Helvetica-Bold',
          fontcolor = black,
          fixedsize = false,
          height = 0.50
          ]
    
    # Default edge styling
    edge [penwidth = 4.0,
          color = '#374151',
          arrowsize = 3.5,
          arrowhead = normal,
          style = solid,
          dir = forward]
    
    # Invisible top node for layout
    invisible_top [label = '', style = invis, height = 0.01, width = 0.01]
    
    # Parent nodes (Core motifs)
    subgraph cluster_0 {
      style = invis
      rank = same
      
      IK [label = 'IK (665)', fillcolor = '#93C5FD']
      KxD [label = 'KxD (419)', fillcolor = '#93C5FD']
      KxE [label = 'KxE (1512)', fillcolor = '#93C5FD']
      VK [label = 'VK (733)', fillcolor = '#93C5FD']
      ExK [label = 'ExK (800)', fillcolor = '#93C5FD']
      KP [label = 'KP (910)', fillcolor = '#93C5FD']
      DxK [label = 'DxK (6730)', fillcolor = '#93C5FD']
    }
    
    # Orphan nodes
    subgraph cluster_1 {
      style = invis
      rank = same
      
      KL [label = 'KL (610)', fillcolor = '#D1D5DB']
      KxS [label = 'KxS (425)', fillcolor = '#D1D5DB']
      KxP [label = 'KxP (301)', fillcolor = '#D1D5DB']
      KxN [label = 'KxN (237)', fillcolor = '#D1D5DB']
      KxQ [label = 'KxQ (291)', fillcolor = '#D1D5DB']
      KxG [label = 'KxG (360)', fillcolor = '#D1D5DB']
    }
    
    # Child nodes
    # IK children
    IKQE [label = 'IKQE (73)', fillcolor = '#BBD6B8']
    IKxE [label = 'IKxE (360)', fillcolor = '#BBD6B8']
    IKxD [label = 'IKxD (61)', fillcolor = '#BBD6B8']
    
    # KxD children
    LKxD [label = 'LKxD (85)', fillcolor = '#BBD6B8']
    
    # KxE children
    {rank = same; 
     TKxE [label = 'TKxE (111)', fillcolor = '#BBD6B8']
     FKxE [label = 'FKxE (77)', fillcolor = '#BBD6B8']
     MKxE [label = 'MKxE (78)', fillcolor = '#BBD6B8']
     PKxE [label = 'PKxE (146)', fillcolor = '#BBD6B8']
    }
    
    {rank = same;
     VKxE [label = 'VKxE (409)', fillcolor = '#BBD6B8']
     LKxE [label = 'LKxE (312)', fillcolor = '#BBD6B8']
     KxEP [label = 'KxEP (273)', fillcolor = '#BBD6B8']
    }
    
    # Other children
    ExKP [label = 'ExKP (181)', fillcolor = '#BBD6B8']
    HTGEKPYK [label = 'HTGEKPYK (224)', fillcolor = '#BBD6B8']
    DxKP [label = 'DxKP (161)', fillcolor = '#BBD6B8']
    
    # Second generation
    IKxEP [label = 'IKxEP (84)', fillcolor = '#BBD6B8']
    LKxEP [label = 'LKxEP (35)', fillcolor = '#BBD6B8']
    PKxExxE [label = 'PKxExxE (36)', fillcolor = '#BBD6B8']
    VKxExxE [label = 'VKxExxE (94)', fillcolor = '#BBD6B8']
    
    # Edges
    invisible_top -> KxE [style = invis]
    
    # IK relationships
    IK -> IKQE
    IK -> IKxE
    IK -> IKxD
    
    # KxD relationships
    KxD -> IKxD
    KxD -> LKxD
    
    # KxE relationships
    KxE -> TKxE
    KxE -> FKxE
    KxE -> MKxE
    KxE -> PKxE
    KxE -> VKxE
    KxE -> LKxE
    KxE -> KxEP
    KxE -> IKxE
    
    # VK relationships
    VK -> VKxE
    
    # ExK relationships
    ExK -> ExKP
    
    # KP relationships
    KP -> ExKP
    KP -> DxKP
    KP -> HTGEKPYK
    
    # DxK relationships
    DxK -> DxKP
    
    # Second level relationships
    IKxE -> IKxEP
    LKxE -> LKxEP
    KxEP -> LKxEP
    KxEP -> IKxEP
    VKxE -> VKxExxE
    PKxE -> PKxExxE
    
    # Invisible edge for layout
    KxE -> KxP [style = invis, weight = 10]
  }
  "
  
  # Create the graph
  motif_graph <- grViz(graph_definition)
  
  # Export as SVG
  motif_svg <- DiagrammeRsvg::export_svg(motif_graph)
  
  return(motif_svg)
}

# Generate the heatmap
generate_heatmap <- function() {
  data_file <- "example_input/data_filtered.csv"
  
  if (!file.exists(data_file)) {
    
    # motifs list
    motifs <- c('DxK', 'DxKP', 'ExK', 'ExKP', 'IK', 'IKxD', 'IKxE', 'IKxEP',
                'KL', 'KP', 'KxD', 'KxE', 'KxEP', 'KxG', 'KxN', 'KxP', 'KxQ', 'KxS', 
                'LKxD', 'LKxE', 'LKxEP', 'MKxE', 'PKxE', 'TKxE', 'VK', 'VKxE')
    
    # Perform GO enrichment analysis
    dfs <- list()
    for (motif in motifs) {
      fileName <- paste0('example_input/', motif, ".txt")
      if (file.exists(fileName)) {
        gene_ids <- readLines(fileName)
        gene_ids <- gsub("^\\w+\\|(.+)\\|\\w+$", "\\1", gene_ids)
        
        universe_ids <- readLines("example_input/Gold_Silver_Bronze.txt")
        universe_ids <- gsub("^\\w+\\|(.+)\\|\\w+$", "\\1", universe_ids)
        
        # GO over-representation test
        enrich_result <- enrichGO(
          gene = gene_ids,
          universe = universe_ids,
          OrgDb = org.Hs.eg.db,
          keyType = "UNIPROT",
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.5,
          qvalueCutoff = 0.5,
          readable = TRUE
        )
        
        # Simplify results
        result_simplified <- clusterProfiler::simplify(
          enrich_result,
          by = "p.adjust",
          cutoff = 0.3,
          select_fun = min
        )
        
        GO_term_dataframe <- as.data.frame(result_simplified) 
        if (nrow(GO_term_dataframe) > 0) {
          GO_term_dataframe$motif_group <- motif
          dfs[[motif]] <- GO_term_dataframe
        }
      }
    }
    
    # Combine data
    df_all <- bind_rows(dfs)
    
    # Filter for significant terms
    data_filtered <- df_all %>% 
      group_by(ID) %>% 
      filter(any(p.adjust <= 0.05)) %>% 
      ungroup() %>%
      mutate(neg_log_p.adjust = -log10(p.adjust)) %>%
      select(Description, motif_group, neg_log_p.adjust) %>%
      spread(motif_group, neg_log_p.adjust) %>%
      column_to_rownames('Description')
    
    write.csv(data_filtered, file = data_file, row.names = TRUE)
  } else {
    data_filtered <- read.csv(data_file, row.names = 1)
  }
  
  # Replace specific GO term description
  rownames(data_filtered) <- gsub("regulation of intracellular steroid hormone receptor signaling pathway", 
                                   "regulation of steroid hormone receptor signalling pathway", 
                                   rownames(data_filtered))
  
  # Read gene counts file
  counts_file <- "example_input/gene_counts.csv"
  display_numbers <- NULL
  
  if (file.exists(counts_file)) {
    Counts <- read.csv(counts_file, row.names = 1)
    
    # Replace specific GO term description in counts as well
    rownames(Counts) <- gsub("regulation of intracellular steroid hormone receptor signaling pathway", 
                             "regulation of steroid hormone receptor signalling pathway", 
                             rownames(Counts))
    
    # Match row and column names
    if (all(rownames(data_filtered) %in% rownames(Counts)) && 
        all(colnames(data_filtered) %in% colnames(Counts))) {
      Counts <- Counts[rownames(data_filtered), colnames(data_filtered)]
      display_numbers <- as.matrix(Counts)
    }
  }
  
  # Replace NA values for custom handling in heatmap
  data_filtered_with_NA <- data_filtered
  data_filtered_with_NA[is.na(data_filtered_with_NA)] <- -1  # Assign -1 for NA values
  
  # Create display numbers matrix with NA handling
  if (!is.null(display_numbers)) {
    display_numbers <- ifelse(
      data_filtered_with_NA == -1, 
      "NA",                      
      ifelse(
        data_filtered_with_NA >= -log10(0.05), 
        display_numbers,           # Keep gene count if significant enrichment
        ""                         # Mask non-significant values
      )
    )
  }
  
  # Heatmap color palette 
  my_palette <- c(
    "grey", # Grey for NA values
    "#ffffcc", # yellow for insignificant values
    colorRampPalette(c("#FFD580", "#FFCC80", "#FFB366", "#FF9955", "#FF8844", 
                       "#FF7733", "#FF6622", "#FF5511", "#FF4400", "#CC2200", 
                       "#990000", "#800026"), space = "Lab")(68)
  )
  
  # Adjust breaks for the color palette
  breaks <- c(
    -1.1, # Value for NA
    -1,   # Boundary for NA
    -log10(0.05),  # Boundary for insignificant enrichment
    seq(-log10(0.05), max(data_filtered_with_NA[data_filtered_with_NA > 0], na.rm = TRUE), length.out = 68)
  )
  
  # Make breaks unique
  breaks <- unique(breaks)
  
  number_color <- matrix("black", nrow = nrow(data_filtered_with_NA), ncol = ncol(data_filtered_with_NA))
  
  # Set specific cells to white for better visibility on dark backgrounds
  if(nrow(number_color) >= 56 && ncol(number_color) >= 1) {
    number_color[56, 1] <- "white"
  }
  
  for(i in 1:nrow(data_filtered_with_NA)) {
    for(j in 1:ncol(data_filtered_with_NA)) {
      if(!is.na(data_filtered_with_NA[i, j]) && data_filtered_with_NA[i, j] > 15) {
        number_color[i, j] <- "white"
      }
    }
  }
  
  # Generate the heatmap as a grob
  heatmap_grob <- pheatmap(
    data_filtered_with_NA, 
    scale = 'none', 
    clustering_distance_rows = 'euclidean', 
    clustering_method = 'complete', 
    show_rownames = TRUE,
    display_numbers = display_numbers, 
    fontsize_row = 8,  
    fontsize_col = 7,    
    fontsize_number = 6, 
    fontsize = 7,       
    legend_breaks = c(0, 5, 10, 15, 20, 25), 
    legend_labels = c("0", "5", "10", "15", "20", "25"),  # Legend labels
    angle_col = 45, 
    color = my_palette,
    breaks = breaks,
    border_color = "grey60",
    number_color = number_color, 
    treeheight_row = 15,  # Dendrogram raw height
    treeheight_col = 15,  # Column dendrogram height
    cellwidth = NA,       # Expand cells to fill available width
    cellheight = NA,      # Expand cells to fill available height
    fontfamily = "sans", 
    silent = TRUE  # Don't display the heatmap immediately
  )
  
  return(heatmap_grob)
}

# Create the combined figure 
create_publication_figure <- function() {
  
  # Generate motif hierarchy
  motif_svg <- create_motif_hierarchy()

  motif_png_data <- rsvg::rsvg_png(charToRaw(motif_svg), width = 12000, height = 3000)
  motif_img <- png::readPNG(motif_png_data)
  motif_grob <- rasterGrob(motif_img, interpolate = TRUE, width = unit(0.95, "npc"), height = unit(0.85, "npc"))
  
  # Generate heatmap
  heatmap_obj <- generate_heatmap()
  
  # Convert heatmap to grob
  heatmap_grob <- heatmap_obj$gtable
  
  # Add panel labels 
  panel_a_labeled <- gTree(children = gList(
    motif_grob,
    textGrob("A", x = 0.02, y = 0.85, 
             just = c("left", "top"), 
             gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "sans"))
  ))
  
  panel_b_labeled <- gTree(children = gList(
    heatmap_grob,
    textGrob("B", x = 0.02, y = 0.99, 
             just = c("left", "top"), 
             gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "sans"))
  ))
  
  # Combine panels
  final_plot <- arrangeGrob(
    panel_a_labeled,
    panel_b_labeled,
    ncol = 1,
    heights = c(0.20, 0.80), 
    padding = unit(0, "line") 
  )
  
  # Get the grob dimensions
  plot_width_inches <- 180 / 25.4  # Convert mm to inches

  
  n_rows <- length(heatmap_obj$tree_row$labels)
  
  n_cols <- length(heatmap_obj$tree_col$labels)
  
  # Set fixed height 
  fixed_height_mm <- 200
  
  # Create output directory if it doesn't exist
  if (!dir.exists("output")) {
    dir.create("output")
  }
  # Save figure
  ggsave("output/integrated_motif_heatmap_figure.png",
         plot = final_plot,
         width = 180,
         height = fixed_height_mm,
         units = "mm",
         dpi = 600,
         bg = "white")
  
  ggsave("output/integrated_motif_heatmap_figure.pdf",
         plot = final_plot,
         width = 180,
         height = fixed_height_mm,
         units = "mm",
         dpi = 600,
         bg = "white")
  

  
  return(final_plot)
}

# Run main function
final_figure <- create_publication_figure()

