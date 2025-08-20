#===============================================================================
# LIGHTNING-FAST VISUALIZATION FUNCTIONS
#===============================================================================

#' Create transcript plot data using pre-computed GTF cache (INSTANT)
#' 
create_fast_transcript_plot_data <- function(gene_id, gene_symbol, as_events = NULL) {
  # Load pre-computed GTF data (0.002 seconds vs 15+ seconds)
  gtf_data <- load_gtf_visualization_data(gene_id)
  
  if (!gtf_data$success) {
    return(list(
      success = FALSE,
      message = gtf_data$message,
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = gtf_data$message, 
                size = 5, hjust = 0.5) +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
    ))
  }
  
  # Extract pre-computed data
  exons_by_transcript <- gtf_data$exons_by_transcript
  cds_by_transcript <- gtf_data$cds_by_transcript
  transcript_ids <- gtf_data$transcript_ids
  
  # Check if we have any transcripts
  if (length(transcript_ids) == 0) {
    return(list(
      success = FALSE,
      message = "No transcripts found for this gene",
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = "No transcripts found for this gene", 
                size = 5, hjust = 0.5) +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
    ))
  }
  
  # Use pre-computed gene boundaries (no calculation needed)
  padding <- 5000
  gene_start <- gtf_data$gene_start - padding
  gene_end <- gtf_data$gene_end + padding
  chromosome <- gtf_data$chromosome
  strand_display <- "5'>3'"
  
  # Create plot data frames efficiently
  transcript_df <- data.frame(
    transcript = transcript_ids,
    y_position = seq_along(transcript_ids),
    stringsAsFactors = FALSE
  )
  
  # Pre-compute exon data
  exon_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    exon_number = integer(),
    stringsAsFactors = FALSE
  )
  
  # Pre-compute CDS data
  cds_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    cds_number = integer(),
    stringsAsFactors = FALSE
  )
  
  # Process transcripts efficiently
  for (i in seq_along(transcript_ids)) {
    tx <- transcript_ids[i]
    
    # Add exons (already filtered by transcript)
    if (tx %in% names(exons_by_transcript) && length(exons_by_transcript[[tx]]) > 0) {
      tx_exons <- exons_by_transcript[[tx]]
      tx_exons <- tx_exons[order(start(tx_exons))]
      
      for (j in seq_along(tx_exons)) {
        exon_df <- rbind(exon_df, data.frame(
          transcript = tx,
          y_position = i,
          start = start(tx_exons[j]),
          end = end(tx_exons[j]),
          exon_number = j,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add CDS (already filtered by transcript)
    if (tx %in% names(cds_by_transcript) && length(cds_by_transcript[[tx]]) > 0) {
      tx_cds <- cds_by_transcript[[tx]]
      tx_cds <- tx_cds[order(start(tx_cds))]
      
      for (j in seq_along(tx_cds)) {
        cds_df <- rbind(cds_df, data.frame(
          transcript = tx,
          y_position = i,
          start = start(tx_cds[j]),
          end = end(tx_cds[j]),
          cds_number = j,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Create the plot (same as before but with pre-computed data)
  max_y_pos <- max(transcript_df$y_position) + 1
  
  p <- ggplot() +
    # Add transcript lines
    geom_segment(data = transcript_df, 
                aes(x = gene_start, xend = gene_end, 
                    y = y_position, yend = y_position),
                linewidth = 0.5, color = "grey70") +
    
    # Add exon blocks
    geom_rect(data = exon_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.3, ymax = y_position + 0.3,
                 fill = "Transcript"),
             color = "black", alpha = 0.8) +
    
    # Add CDS overlay
    geom_rect(data = cds_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.25, ymax = y_position + 0.25,
                 fill = "CDS"),
             color = "black", alpha = 0.9) +
    
    # Add transcript labels
    geom_text(data = transcript_df,
             aes(x = gene_start - 100, y = y_position, label = transcript),
             hjust = 1, size = 3.5)
  
  # Add direction indicator
  direction_df <- data.frame(
    x = c(gene_start + padding, gene_end - padding),
    y = c(max_y_pos, max_y_pos),
    label = c("5'", "3'"),
    stringsAsFactors = FALSE
  )
  
  p <- p + 
    geom_text(data = direction_df[1,], aes(x = x, y = y, label = label), size = 4, fontface = "bold") +
    geom_text(data = direction_df[2,], aes(x = x, y = y, label = label), size = 4, fontface = "bold") +
    geom_segment(data = data.frame(x = gene_start + padding + 50, xend = gene_end - padding - 50, 
                                  y = max_y_pos, yend = max_y_pos),
                aes(x = x, xend = xend, y = y, yend = yend),
                arrow = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"),
                color = "black", linewidth = 0.7)
  
  # Add fill scale for transcript structure - only show relevant elements
  has_cds <- nrow(cds_df) > 0
  fill_values <- c(
    "Transcript" = "rgba(77, 175, 74, 0.8)",
    "CDS" = "rgba(255, 221, 0, 0.8)"
  )
  
  p <- p + 
    scale_fill_manual(
      values = fill_values,
      breaks = c("Transcript", "CDS"),
      labels = c("Transcript", "CDS"),
      name = ""
    ) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(6, 6, 6, 6)
    ) +
    labs(
      x = paste0("Genomic Position (chromosome ", chromosome, ") - ", strand_display), 
      title = paste0("Transcript Structure - ", gene_symbol, " (", gene_id, ")"),
      subtitle = paste0("Chromosome ", chromosome, ": ", gene_start, " - ", gene_end)
    ) +
    coord_cartesian(xlim = c(gene_start, gene_end), 
                   ylim = c(0.5, max_y_pos + 0.5))
  
  return(list(
    success = TRUE,
    plot = p,
    chromosome = chromosome,
    has_cds = has_cds
  ))
}

#' Create peptide plot data using pre-computed GTF cache (INSTANT)
#' 
create_fast_peptide_plot_data <- function(gene_id, gene_symbol, processed_data, protease) {
  # Load pre-computed GTF data (0.002 seconds vs 15+ seconds)
  gtf_data <- load_gtf_visualization_data(gene_id)
  
  if (!gtf_data$success) {
    return(list(
      success = FALSE,
      message = gtf_data$message,
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = gtf_data$message, 
                size = 5, hjust = 0.5) +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
    ))
  }
  
  # Extract pre-computed data
  exons_by_transcript <- gtf_data$exons_by_transcript
  cds_by_transcript <- gtf_data$cds_by_transcript
  transcript_ids <- gtf_data$transcript_ids
  
  # Use pre-computed gene boundaries
  padding <- 5000
  gene_start <- gtf_data$gene_start - padding
  gene_end <- gtf_data$gene_end + padding
  chromosome <- gtf_data$chromosome
  strand_display <- "5'>3'"
  
  # Create plot data frames (same structure as before)
  transcript_df <- data.frame(
    transcript = transcript_ids,
    y_position = seq_along(transcript_ids),
    stringsAsFactors = FALSE
  )
  
  # Pre-compute all data frames
  exon_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    exon_number = integer(),
    stringsAsFactors = FALSE
  )
  
  cds_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    cds_number = integer(),
    stringsAsFactors = FALSE
  )
  
  peptide_df <- data.frame(
    transcript = character(),
    y_position = numeric(),
    start = numeric(),
    end = numeric(),
    peptide = character(),
    is_junction_spanning = logical(),
    hover_text = character(),
    stringsAsFactors = FALSE
  )
  
  # Process efficiently (same logic as before but with pre-computed data)
  for (i in seq_along(transcript_ids)) {
    tx <- transcript_ids[i]
    
    # Add exons
    if (tx %in% names(exons_by_transcript) && length(exons_by_transcript[[tx]]) > 0) {
      tx_exons <- exons_by_transcript[[tx]]
      tx_exons <- tx_exons[order(start(tx_exons))]
      
      for (j in seq_along(tx_exons)) {
        exon_df <- rbind(exon_df, data.frame(
          transcript = tx,
          y_position = i,
          start = start(tx_exons[j]),
          end = end(tx_exons[j]),
          exon_number = j,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add CDS
    if (tx %in% names(cds_by_transcript) && length(cds_by_transcript[[tx]]) > 0) {
      tx_cds <- cds_by_transcript[[tx]]
      tx_cds <- tx_cds[order(start(tx_cds))]
      
      for (j in seq_along(tx_cds)) {
        cds_df <- rbind(cds_df, data.frame(
          transcript = tx,
          y_position = i,
          start = start(tx_cds[j]),
          end = end(tx_cds[j]),
          cds_number = j,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add peptides
    peptide_data <- get_transcript_peptides_genomic(tx, processed_data, protease)
    if (!is.null(peptide_data) && length(peptide_data$genomic_ranges) > 0) {
      tx_peptides <- peptide_data$genomic_ranges
      
      # Count peptide occurrences
      peptide_counts <- table(tx_peptides$peptide)
      junction_spanning <- names(peptide_counts[peptide_counts > 1])
      
      # Add all peptide ranges
      peptide_df <- rbind(peptide_df, data.frame(
        transcript = tx,
        y_position = i,
        start = start(tx_peptides),
        end = end(tx_peptides),
        peptide = tx_peptides$peptide,
        is_junction_spanning = tx_peptides$peptide %in% junction_spanning,
        hover_text = clean_hover_text(paste0(
          "Peptide: ", tx_peptides$peptide,
          "<br>Position: ", start(tx_peptides), "-", end(tx_peptides),
          "<br>Transcript: ", tx,
          "<br>Type: ", ifelse(tx_peptides$peptide %in% junction_spanning, "Junction-spanning", "Regular")
        )),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Create the plot (same as before)
  p <- ggplot() +
    geom_segment(data = transcript_df, 
                aes(x = gene_start, xend = gene_end, 
                    y = y_position, yend = y_position),
                linewidth = 0.5, color = "grey70") +
    
    geom_rect(data = exon_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.3, ymax = y_position + 0.3,
                 fill = "Transcript"),
             color = "black", alpha = 0.8) +
    
    geom_rect(data = cds_df,
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.25, ymax = y_position + 0.25,
                 fill = "CDS"),
             color = "black", alpha = 0.9) +
    
    # Add peptides with fill aesthetic for legend consistency
    geom_rect(data = peptide_df[!peptide_df$is_junction_spanning, ],
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.15, ymax = y_position + 0.15,
                 text = hover_text, fill = "peptide"),
             color = "black", alpha = 0.9) +
    geom_rect(data = peptide_df[peptide_df$is_junction_spanning, ],
             aes(xmin = start, xmax = end, 
                 ymin = y_position - 0.15, ymax = y_position + 0.15,
                 text = hover_text, fill = "junction spanning peptide"),
             color = "black", alpha = 0.9) +
    
    geom_text(data = transcript_df,
             aes(x = gene_start - 100, y = y_position, label = transcript),
             hjust = 1, size = 3.5)
  
  # Add direction indicator and styling (same as before)
  direction_df <- data.frame(
    x = c(gene_start + padding, gene_end - padding),
    y = c(max(transcript_df$y_position) + 1, max(transcript_df$y_position) + 1),
    label = c("5'", "3'"),
    stringsAsFactors = FALSE
  )
  
  p <- p + 
    geom_text(data = direction_df[1,], aes(x = x, y = y, label = label), size = 4, fontface = "bold") +
    geom_text(data = direction_df[2,], aes(x = x, y = y, label = label), size = 4, fontface = "bold") +
    geom_segment(data = data.frame(x = gene_start + padding + 50, xend = gene_end - padding - 50, 
                                  y = max(transcript_df$y_position) + 0.5, yend = max(transcript_df$y_position) + 0.5),
                aes(x = x, xend = xend, y = y, yend = yend),
                arrow = arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"),
                color = "black", linewidth = 0.7) +
    
    # Use the exact same unified fill scale as transcript structure plot
    {
      fill_values <- c(
        "Transcript" = "rgba(77, 175, 74, 0.8)",
        "CDS" = "rgba(255, 221, 0, 0.8)",
        "peptide" = "rgba(52, 152, 219, 0.9)",
        "junction spanning peptide" = "rgba(231, 76, 60, 0.9)"
      )
      
      scale_fill_manual(
        values = fill_values,
        breaks = c("Transcript", "CDS", "peptide", "junction spanning peptide"),
        labels = c("Transcript", "CDS", "peptide", "junction spanning peptide"),
        name = ""
      )
    } +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(6, 6, 6, 6)
    ) +
    labs(
      x = paste0("Genomic Position (chromosome ", chromosome, ") - ", strand_display), 
      title = paste0("Peptide Visualization - ", gene_symbol, " (", gene_id, ")"),
      subtitle = paste0("Chromosome ", chromosome, ": ", gene_start, " - ", gene_end, 
                       " | Enzyme: ", protease)
    ) +
    coord_cartesian(xlim = c(gene_start, gene_end), 
                   ylim = c(0.5, max(transcript_df$y_position) + 0.5))
  
  has_peptides <- nrow(peptide_df) > 0
  if (!has_peptides) {
    p <- p + annotate("text", x = (gene_start + gene_end)/2, y = max(transcript_df$y_position)/2, 
                     label = paste0("No mapped peptides found for enzyme: ", protease),
                     size = 5, color = "red")
  }
  
  return(list(
    success = TRUE,
    plot = p,
    chromosome = chromosome,
    has_peptides = has_peptides
  ))
}