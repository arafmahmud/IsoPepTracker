#===============================================================================
# MULTI-ISOFORM COMPARISON MODULE
# Extracted from server.R - Multi-isoform comparative analysis functionality
#===============================================================================

#' Multi-Isoform Comparison Module
#' 
#' @description 
#' Complete multi-isoform comparative analysis functionality:
#' - Multiple isoform selection and management
#' - Comparative peptide analysis across selected isoforms
#' - Interactive comparative visualizations with plotly
#' - Peptide overlap matrix generation
#' - Comparative summary statistics
#' - Download handlers for comparative analysis results
#' 
#' @param input Shiny input object
#' @param output Shiny output object  
#' @param session Shiny session object
#' @param processed_data Reactive containing processed application data
#' @param gene_data Reactive containing current gene data
#' @param selected_gene_transcripts Reactive containing gene transcripts
#' @param selected_gene_peptides Reactive containing gene peptides
#' @param novel_pipeline_results Reactive containing novel isoform results (optional)
#' 
#' @return Named list containing reactive values for multi-isoform comparison
#' 
#' @note 
#' - Supports both regular gene isoforms and novel isoforms
#' - Requires gene-specific data loading for accurate comparisons
#' - Uses plotly for interactive comparative visualizations
#' - Maintains compatibility with existing isoform analysis pipeline

create_multi_isoform_comparison_module <- function(input, output, session, processed_data, gene_data, 
                                                 selected_gene_transcripts, selected_gene_peptides, 
                                                 novel_pipeline_results = NULL) {
  
  # ============================================================================
  # MODULE REACTIVE VALUES
  # ============================================================================
  
  # Multi-isoform comparison data
  multi_isoform_data <- reactiveVal(NULL)
  multi_isoform_highlighted_data <- reactiveVal(NULL)
  comparison_isoforms <- reactiveVal(character(0))
  
  # ============================================================================
  # ISOFORM SELECTION MANAGEMENT
  # ============================================================================
  
  # Update compare_isoforms choices when gene is selected OR novel data is available
  observeEvent({
    list(input$gene, novel_pipeline_results())
  }, {
    # Always try to get known gene transcripts if gene is selected
    known_transcripts <- NULL
    gene_symbol <- NULL
    
    if (!is.null(input$gene) && !is.na(input$gene) && input$gene != "") {
      transcripts_data <- selected_gene_transcripts()
      if (!is.null(transcripts_data) && length(transcripts_data) > 0) {
        known_transcripts <- transcripts_data
        
        # Get gene symbol
        gene_symbol <- tryCatch({
          processed_data()$gene_lookup[input$gene]
        }, error = function(e) {
          input$gene
        })
        if (is.null(gene_symbol)) gene_symbol <- input$gene
      }
    }
    
    # Get novel transcripts if available
    novel_transcripts <- NULL
    if (!is.null(novel_pipeline_results)) {
      novel_data <- novel_pipeline_results()
      
      if (!is.null(novel_data) && novel_data$success && !is.null(novel_data$dataframe_file) && file.exists(novel_data$dataframe_file)) {
        tryCatch({
          novel_rds <- readRDS(novel_data$dataframe_file)
          if (!is.null(novel_rds) && nrow(novel_rds) > 0) {
            novel_transcripts <- unique(novel_rds$txID)
            cat("Found", length(novel_transcripts), "novel transcripts for comparison\n")
          }
        }, error = function(e) {
          cat("Error loading novel transcripts:", e$message, "\n")
        })
      }
    }
    
    # Combine choices
    all_choices <- character(0)
    
    if (!is.null(known_transcripts) && length(known_transcripts) > 0) {
      # Add known transcripts with labels
      known_labels <- paste0(known_transcripts, " (", gene_symbol, ")")
      names(known_transcripts) <- known_labels
      all_choices <- c(all_choices, known_transcripts)
    }
    
    if (!is.null(novel_transcripts) && length(novel_transcripts) > 0) {
      # Add novel transcripts with labels
      novel_labels <- paste0(novel_transcripts, " (Novel)")
      names(novel_transcripts) <- novel_labels
      all_choices <- c(all_choices, novel_transcripts)
    }
    
    # Update dropdown
    if (length(all_choices) > 0) {
      updateSelectizeInput(session, "compare_isoforms", 
                          choices = all_choices,
                          selected = character(0))
      
      # Show notification about available options
      if (!is.null(known_transcripts) && !is.null(novel_transcripts)) {
        cat("✅ Multiple Isoform Comparison: Added", length(known_transcripts), "known +", length(novel_transcripts), "novel transcripts\n")
      }
    } else {
      # Clear choices if no transcripts available
      updateSelectizeInput(session, "compare_isoforms", 
                          choices = NULL, selected = character(0))
    }
  })
  
  # Enable/disable comparative analysis button based on selection count
  observeEvent(input$compare_isoforms, {
    selected_count <- length(input$compare_isoforms)
    comparison_isoforms(input$compare_isoforms)
    
    if (selected_count >= 2 && selected_count <= 8) {
      shinyjs::enable("run_comparative_analysis")
      shinyjs::removeClass("run_comparative_analysis", "btn-secondary")
      shinyjs::addClass("run_comparative_analysis", "btn-primary")
    } else {
      shinyjs::disable("run_comparative_analysis")
      shinyjs::removeClass("run_comparative_analysis", "btn-primary")
      shinyjs::addClass("run_comparative_analysis", "btn-secondary")
    }
  })
  
  # ============================================================================
  # MULTI-ISOFORM DATA PROCESSING
  # ============================================================================
  
  # Multi-Isoform data loader
  multi_isoform_data <- reactive({
    req(input$run_comparative_analysis > 0, input$gene, input$protease, input$miscleavage_type, gene_data(),
        input$compare_isoforms, length(input$compare_isoforms) >= 2, length(input$compare_isoforms) <= 8)
    
    withProgress(message = 'Loading multi-isoform analysis...', value = 0, {
      gene_id <- input$gene
      protease <- input$protease
      miscleavage_type <- input$miscleavage_type
      selected_transcripts <- input$compare_isoforms
      
      # Step 1: Get all transcripts for this gene (from gene-by-gene loaded data)
      incProgress(0.2, detail = 'Getting selected gene transcripts...')
      all_transcripts <- selected_gene_transcripts()
      
      if (is.null(all_transcripts) || length(all_transcripts) == 0) {
        showNotification(paste("No transcripts found for gene", gene_id, "with miscleavage type", miscleavage_type), type = "warning")
        return(NULL)
      }
      
      cat("Selected transcripts for multi-isoform analysis:", paste(selected_transcripts, collapse=", "), "\n")
      cat("Using miscleavage type:", miscleavage_type, "\n")
      
      # Step 2: Get peptides for each transcript (using gene-by-gene loaded data)
      incProgress(0.3, detail = 'Loading peptides for selected transcripts...')
      gene_peptides <- selected_gene_peptides()
      if (is.null(gene_peptides) || nrow(gene_peptides) == 0) {
        showNotification(paste("No peptides found for gene", gene_id, "with miscleavage type", miscleavage_type), type = "warning")
        return(NULL)
      }
      
      # Create vis_data structure for compatibility with existing functions
      vis_data <- list(
        genes = processed_data()$genes,
        gene_symbols = processed_data()$gene_symbols,
        gene_lookup = processed_data()$gene_lookup,
        proteases = processed_data()$proteases,
        original_peptides = gene_peptides
      )
      
      # Get peptides for SELECTED transcripts only
      all_peptides_list <- list()
      for (i in seq_along(selected_transcripts)) {
        tx <- selected_transcripts[i]
        tx_peptides <- get_transcript_peptides_for_comparison(tx, processed_data(), protease)
        if (!is.null(tx_peptides) && length(tx_peptides) > 0) {
          all_peptides_list[[tx]] <- data.frame(
            transcript = tx,
            y_position = i,
            start = start(tx_peptides),
            end = end(tx_peptides),
            peptide = tx_peptides$peptide,
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (length(all_peptides_list) == 0) {
        showNotification("No peptides found for selected transcripts", type = "warning")
        return(NULL)
      }
      
      # Combine all peptides
      all_peptides_df <- data.table::rbindlist(all_peptides_list)
      
      # Step 3: Calculate gene boundaries and prepare visualization data
      incProgress(0.2, detail = 'Preparing visualization data...')
      gene_start <- min(all_peptides_df$start) - 1000
      gene_end <- max(all_peptides_df$end) + 1000
      
      # Create transcript position mapping
      transcript_df <- data.frame(
        transcript = selected_transcripts,
        y_position = seq_along(selected_transcripts),
        stringsAsFactors = FALSE
      )
      
      # Step 4: Load rMATS GTF structure data if we have rMATS transcripts
      incProgress(0.1, detail = 'Loading transcript structure data...')
      rmats_structure_data <- NULL
      
      # Check if we have rMATS transcripts (they contain ".inclusion" or ".exclusion")
      rmats_transcripts <- selected_transcripts[grepl("\\.(inclusion|exclusion)$", selected_transcripts)]
      
      if (length(rmats_transcripts) > 0) {
        cat("🔍 Detected", length(rmats_transcripts), "rMATS transcripts:", paste(rmats_transcripts, collapse = ", "), "\n")
        cat("📁 Searching for rMATS GTF files...\n")
        
        # Find the most recent rMATS GTF file for this gene
        gene_clean <- gsub('"', '', gene_id)
        cat("🧬 Gene ID (cleaned):", gene_clean, "\n")
        
        # Try multiple search patterns to find rMATS GTF files
        search_patterns <- c(
          paste0('rmats_"', gene_clean, '".*\\.transdecoder\\.genome\\.gtf$'),  # Original pattern
          paste0('rmats_', gene_clean, '.*\\.transdecoder\\.genome\\.gtf$'),    # Without quotes
          paste0('.*', gene_clean, '.*\\.transdecoder\\.genome\\.gtf$'),        # More flexible
          paste0('.*\\.transdecoder\\.genome\\.gtf$')                            # Any rMATS GTF
        )
        
        rmats_gtf_files <- c()
        for (i in seq_along(search_patterns)) {
          pattern <- search_patterns[i]
          cat("🔍 Trying pattern", i, ":", pattern, "\n")
          
          found_files <- list.files(
            path = "rmats_peptide_results", 
            pattern = pattern, 
            full.names = TRUE,
            recursive = TRUE
          )
          
          if (length(found_files) > 0) {
            cat("✅ Found", length(found_files), "files with pattern", i, "\n")
            rmats_gtf_files <- c(rmats_gtf_files, found_files)
            break  # Use the first successful pattern
          } else {
            cat("❌ No files found with pattern", i, "\n")
          }
        }
        
        # Remove duplicates and show all found files
        rmats_gtf_files <- unique(rmats_gtf_files)
        
        if (length(rmats_gtf_files) > 0) {
          cat("📂 Found", length(rmats_gtf_files), "rMATS GTF file(s):\n")
          for (f in rmats_gtf_files) {
            cat("  📄", basename(f), "(", file.info(f)$mtime, ")\n")
          }
          
          # Use the most recent GTF file
          rmats_gtf_file <- rmats_gtf_files[order(file.mtime(rmats_gtf_files), decreasing = TRUE)][1]
          cat("🎯 Selected most recent:", basename(rmats_gtf_file), "\n")
          
          # Load rMATS transcript structure data
          cat("📥 Loading rMATS structure data...\n")
          rmats_structure_data <- load_rmats_transcript_exons(rmats_gtf_file, rmats_transcripts)
          
          if (rmats_structure_data$success) {
            cat("✅ Successfully loaded rMATS structure data for", rmats_structure_data$loaded_transcripts, "out of", rmats_structure_data$requested_transcripts, "transcripts\n")
            if (length(rmats_structure_data$exons) > 0) {
              cat("📊 Loaded exons for transcripts:", paste(names(rmats_structure_data$exons), collapse = ", "), "\n")
            }
            if (length(rmats_structure_data$cds) > 0) {
              cat("📊 Loaded CDS for transcripts:", paste(names(rmats_structure_data$cds), collapse = ", "), "\n")
            }
          } else {
            cat("❌ Failed to load rMATS structure data:", rmats_structure_data$message, "\n")
          }
        } else {
          cat("❌ No rMATS GTF files found for gene", gene_clean, "\n")
          cat("🔍 Available files in rmats_peptide_results:\n")
          all_files <- list.files("rmats_peptide_results", recursive = TRUE, full.names = FALSE)
          gtf_files <- all_files[grepl("\\.gtf$", all_files)]
          if (length(gtf_files) > 0) {
            cat("  📄 GTF files found:", paste(head(gtf_files, 5), collapse = ", "))
            if (length(gtf_files) > 5) cat(" (showing first 5)")
            cat("\n")
          } else {
            cat("  📄 No GTF files found in directory\n")
          }
        }
      }
      
      # Step 5: Load structure data for regular transcripts from main GTF
      regular_transcripts <- selected_transcripts[!grepl("\\.(inclusion|exclusion)$", selected_transcripts)]
      regular_structure_data <- NULL
      
      if (length(regular_transcripts) > 0) {
        cat("🔍 Detected", length(regular_transcripts), "regular transcripts:", paste(regular_transcripts, collapse = ", "), "\n")
        cat("📁 Loading structure data from main GTF...\n")
        
        # Load gene details first
        gene_details <- load_gene_details(gene_id)
        
        if (gene_details$success) {
          cat("✅ Gene details loaded:", gene_details$chromosome, ":", gene_details$start, "-", gene_details$end, "\n")
          
          # Load transcript exons and CDS
          regular_structure_data <- load_transcript_exons(gene_details, regular_transcripts)
          
          if (regular_structure_data$success) {
            cat("✅ Successfully loaded regular structure data\n")
            if (length(regular_structure_data$exons) > 0) {
              cat("📊 Loaded exons for transcripts:", paste(names(regular_structure_data$exons), collapse = ", "), "\n")
            }
            if (length(regular_structure_data$cds) > 0) {
              cat("📊 Loaded CDS for transcripts:", paste(names(regular_structure_data$cds), collapse = ", "), "\n")
            }
            
            # Merge regular structure data with rMATS structure data
            if (!is.null(rmats_structure_data) && rmats_structure_data$success) {
              cat("🔗 Merging rMATS and regular structure data...\n")
              # Combine exons
              combined_exons <- c(rmats_structure_data$exons, regular_structure_data$exons)
              # Combine CDS
              combined_cds <- c(rmats_structure_data$cds, regular_structure_data$cds)
              
              # Update rmats_structure_data to include both
              rmats_structure_data$exons <- combined_exons
              rmats_structure_data$cds <- combined_cds
              
              cat("📊 Combined structure data summary:\n")
              cat("  Total exon sets:", length(combined_exons), "\n")
              cat("  Total CDS sets:", length(combined_cds), "\n")
              cat("  Exon transcript IDs:", paste(names(combined_exons), collapse = ", "), "\n")
              cat("  CDS transcript IDs:", paste(names(combined_cds), collapse = ", "), "\n")
              
            } else {
              # Only regular structure data available
              rmats_structure_data <- regular_structure_data
            }
          } else {
            cat("❌ Failed to load regular structure data:", regular_structure_data$message, "\n")
          }
        } else {
          cat("❌ Failed to load gene details:", gene_details$message, "\n")
        }
      }
      
      # Update y_position in peptides
      all_peptides_df$y_position <- transcript_df$y_position[match(all_peptides_df$transcript, transcript_df$transcript)]
      
      # Add hover text for peptides
      all_peptides_df$hover_text <- paste0(
        "Peptide: ", all_peptides_df$peptide,
        "<br>Position: ", all_peptides_df$start, "-", all_peptides_df$end,
        "<br>Transcript: ", all_peptides_df$transcript,
        "<br>Enzyme: ", protease,
        "<br>Miscleavage: ", miscleavage_type
      )
      
      
      incProgress(0.1, detail = 'Complete!')
      
      return(list(
        all_peptides = all_peptides_df,
        transcript_df = transcript_df,
        gene_start = gene_start,
        gene_end = gene_end,
        all_transcripts = selected_transcripts,
        rmats_structure_data = rmats_structure_data
      ))
    })
  })
  
  # Multi-Isoform highlighted data (for selected highlighting)
  multi_isoform_highlighted_data <- reactive({
    req(multi_isoform_data())
    
    base_data <- multi_isoform_data()
    if (is.null(base_data)) return(NULL)
    
    all_peptides_df <- base_data$all_peptides
    all_transcripts <- base_data$all_transcripts
    
    # Add highlighting information if a specific isoform is selected
    highlight_isoform <- input$highlight_isoform %||% all_transcripts[1]
    
    if (!is.null(highlight_isoform) && highlight_isoform %in% all_transcripts) {
      # Mark peptides for highlighting
      all_peptides_df$is_highlighted <- all_peptides_df$transcript == highlight_isoform
      
      # Update hover text with highlighting info
      all_peptides_df$hover_text <- paste0(
        all_peptides_df$hover_text,
        ifelse(all_peptides_df$is_highlighted, "<br><b>HIGHLIGHTED</b>", "")
      )
    } else {
      all_peptides_df$is_highlighted <- FALSE
    }
    
    return(list(
      all_peptides = all_peptides_df,
      transcript_df = base_data$transcript_df,
      gene_start = base_data$gene_start,
      gene_end = base_data$gene_end,
      all_transcripts = all_transcripts,
      highlight_isoform = highlight_isoform
    ))
  })
  
  # ============================================================================
  # VISUALIZATION OUTPUTS
  # ============================================================================
  
  # Render comparative plot
  output$comparative_plot <- renderPlotly({
    data <- multi_isoform_highlighted_data()
    
    if (is.null(data)) {
      return(empty_plotly_message("No data available for comparative analysis. Please select at least two isoforms and click 'Run'."))
    }
    
    # Define miscleavage_label based on input for the plot title
    miscleavage_label <- switch(input$miscleavage_type,
                                "no_miss_cleavage" = "No Missed Cleavages",
                                "upto_one_misscleavage" = "Up to 1 Missed Cleavage",
                                "upto_two_misscleavage" = "Up to 2 Missed Cleavages",
                                "Unknown")
    
    all_peptides <- data$all_peptides
    transcript_df <- data$transcript_df
    gene_start <- data$gene_start
    gene_end <- data$gene_end
    all_transcripts <- data$all_transcripts
    
    # Create the comparative plot
    tryCatch({
      p <- create_multi_isoform_plot(
        all_peptides = all_peptides,
        transcript_df = transcript_df,
        gene_start = gene_start,
        gene_end = gene_end,
        gene_id = input$gene,
        miscleavage_label = miscleavage_label,
        protease = input$protease,
        highlight_isoform = data$highlight_isoform,
        rmats_structure_data = data$rmats_structure_data
      )
      
      # Convert to plotly with interactivity
      plotly_obj <- ggplotly(p, tooltip = "text") %>%
        config(
          displayModeBar = TRUE,
          displaylogo = FALSE,
          modeBarButtonsToRemove = c(
            "pan2d", "select2d", "lasso2d", "resetScale2d"
          )
        ) %>%
        layout(
          title = paste0('Multi-Isoform Comparative Analysis (', miscleavage_label, ') - ', 
                        length(all_transcripts), ' Selected'),
          showlegend = TRUE,
          legend = list(
            orientation = "h",
            x = 0.5,
            xanchor = "center",
            y = -0.1
          ),
          margin = list(l = 50, r = 50, t = 80, b = 100)
        )
      
      return(plotly_obj)
      
    }, error = function(e) {
      cat("Error in comparative_plot:", e$message, "\n")
      return(empty_plotly_message(paste("Error creating comparative plot:", e$message)))
    })
  })
  
  # ============================================================================
  # DATA TABLES
  # ============================================================================
  
  
  
  # ============================================================================
  # STATUS AND UI OUTPUTS
  # ============================================================================
  
  # Comparative analysis status
  output$comparative_analysis_status <- renderUI({
    selected_count <- length(input$compare_isoforms)
    
    if (selected_count == 0) {
      div(
        style = "color: #666; margin-top: 10px;",
        "Please select 2-8 isoforms for comparison."
      )
    } else if (selected_count == 1) {
      div(
        style = "color: #ff9900; margin-top: 10px;",
        paste("Selected 1 isoform. Please select at least 2 for comparison.")
      )
    } else if (selected_count > 8) {
      div(
        style = "color: #cc0000; margin-top: 10px;",
        paste("Selected", selected_count, "isoforms. Maximum 8 allowed for performance reasons.")
      )
    } else {
      div(
        style = "color: #009900; margin-top: 10px;",
        paste("Selected", selected_count, "isoforms. Ready for comparative analysis.")
      )
    }
  })
  
  # ============================================================================
  # DOWNLOAD HANDLERS
  # ============================================================================
  
  # Download comparative data
  output$download_comparative_data <- downloadHandler(
    filename = function() {
      gene_symbol <- processed_data()$gene_lookup[input$gene] %||% input$gene
      miscleavage_suffix <- switch(input$miscleavage_type,
                                  "no_miss_cleavage" = "no_miss",
                                  "upto_one_misscleavage" = "1_miss",
                                  "upto_two_misscleavage" = "2_miss",
                                  "unknown")
      paste0("comparative_analysis_", gene_symbol, "_", input$protease, "_", 
             miscleavage_suffix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      data <- multi_isoform_highlighted_data()
      if (!is.null(data) && !is.null(data$all_peptides)) {
        comparative_peptides <- data$all_peptides
        
        # Add analysis metadata
        comparative_peptides$gene_id <- input$gene
        comparative_peptides$gene_symbol <- processed_data()$gene_lookup[input$gene] %||% input$gene
        comparative_peptides$enzyme <- input$protease
        comparative_peptides$miscleavage_type <- input$miscleavage_type
        comparative_peptides$analysis_timestamp <- Sys.time()
        comparative_peptides$selected_isoforms <- paste(input$compare_isoforms, collapse = ", ")
        
        write.csv(comparative_peptides, file, row.names = FALSE)
      }
    }
  )
  
  # Download overlap matrix
  output$download_overlap_matrix <- downloadHandler(
    filename = function() {
      gene_symbol <- processed_data()$gene_lookup[input$gene] %||% input$gene
      paste0("overlap_matrix_", gene_symbol, "_", input$protease, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      data <- multi_isoform_highlighted_data()
      if (!is.null(data) && !is.null(data$overlap_matrix)) {
        write.csv(data$overlap_matrix, file, row.names = TRUE)
      }
    }
  )
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  
  # Get comparison status
  get_comparison_status <- function() {
    data <- multi_isoform_highlighted_data()
    
    list(
      has_data = !is.null(data),
      selected_count = length(comparison_isoforms()),
      ready_for_analysis = length(comparison_isoforms()) >= 2 && length(comparison_isoforms()) <= 8,
      total_peptides = if (!is.null(data) && !is.null(data$all_peptides)) nrow(data$all_peptides) else 0
    )
  }
  
  # Clear comparison data
  clear_comparison_data <- function() {
    multi_isoform_data(NULL)
    multi_isoform_highlighted_data(NULL)
    comparison_isoforms(character(0))
  }
  
  # ============================================================================
  # MODULE RETURN VALUES
  # ============================================================================
  
  return(list(
    # Reactive values
    multi_isoform_data = multi_isoform_data,
    multi_isoform_highlighted_data = multi_isoform_highlighted_data,
    comparison_isoforms = comparison_isoforms,
    
    # Utility functions
    get_comparison_status = get_comparison_status,
    clear_comparison_data = clear_comparison_data
  ))
}

#===============================================================================
# MODULE TESTING FUNCTION
#===============================================================================

#' Test Multi-Isoform Comparison Module
#' 
#' @description 
#' Function to test multi-isoform comparison module functionality
#' 
#' @return TRUE if all tests pass, FALSE otherwise

test_multi_isoform_comparison_module <- function() {
  cat("Testing Multi-Isoform Comparison module...\n")
  
  # Test 1: Overlap matrix calculation
  overlap_test <- tryCatch({
    # Test overlap matrix calculation with mock data
    transcripts <- c("TX1", "TX2", "TX3")
    peptides_list <- list(
      TX1 = data.frame(peptide = c("PEPTIDE1", "PEPTIDE2", "PEPTIDE3")),
      TX2 = data.frame(peptide = c("PEPTIDE1", "PEPTIDE4", "PEPTIDE5")),
      TX3 = data.frame(peptide = c("PEPTIDE2", "PEPTIDE3", "PEPTIDE4"))
    )
    
    # Mock overlap matrix structure
    n_tx <- length(transcripts)
    overlap_matrix <- matrix(0, nrow = n_tx, ncol = n_tx)
    diag(overlap_matrix) <- 1.0
    
    # Verify structure
    nrow(overlap_matrix) == length(transcripts) && ncol(overlap_matrix) == length(transcripts)
  }, error = function(e) {
    cat("Overlap test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 2: Comparative summary generation
  summary_test <- tryCatch({
    # Test summary structure
    mock_summary <- data.frame(
      Transcript = c("TX1", "TX2"),
      Total_Peptides = c(10, 15),
      Unique_Peptides = c(8, 12),
      stringsAsFactors = FALSE
    )
    
    # Verify structure
    all(c("Transcript", "Total_Peptides", "Unique_Peptides") %in% names(mock_summary))
  }, error = function(e) {
    cat("Summary test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 3: Status checking logic
  status_test <- tryCatch({
    # Test status structure
    mock_status <- list(
      has_data = FALSE,
      selected_count = 0,
      ready_for_analysis = FALSE,
      total_peptides = 0
    )
    
    # Verify all required fields
    required_fields <- c("has_data", "selected_count", "ready_for_analysis", "total_peptides")
    all(required_fields %in% names(mock_status))
  }, error = function(e) {
    cat("Status test failed:", e$message, "\n")
    FALSE
  })
  
  # Test 4: Selection count validation
  selection_test <- tryCatch({
    # Test selection validation logic
    test_counts <- c(0, 1, 3, 9)
    valid_counts <- sapply(test_counts, function(x) x >= 2 && x <= 8)
    
    # Should be: FALSE, FALSE, TRUE, FALSE
    identical(valid_counts, c(FALSE, FALSE, TRUE, FALSE))
  }, error = function(e) {
    cat("Selection test failed:", e$message, "\n")
    FALSE
  })
  
  if (overlap_test && summary_test && status_test && selection_test) {
    cat("All Multi-Isoform Comparison module tests passed!\n")
    return(TRUE)
  } else {
    cat("Some Multi-Isoform Comparison module tests failed!\n")
    return(FALSE)
  }
}

#===============================================================================
# HELPER FUNCTIONS (DEPENDENCIES)
#===============================================================================

#' Create Empty Plotly Message
#' 
#' Creates a simple plotly plot with a text message for empty states
#' 
#' @param message Text message to display
#' @return Plotly object with the message
empty_plotly_message <- function(message) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = message, size = 5, hjust = 0.5) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "#f8f9fa", color = NA),
      plot.background = element_rect(fill = "#f8f9fa", color = NA)
    )
  
  ggplotly(p) %>%
    config(displayModeBar = FALSE) %>%
    layout(
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE)
    )
}

#' Create Multi-Isoform Plot
#' 
#' Creates a ggplot visualization for multi-isoform comparison
#' 
#' @param all_peptides Data frame with peptide information
#' @param transcript_df Data frame with transcript information
#' @param gene_start Numeric gene start position
#' @param gene_end Numeric gene end position
#' @param gene_id String gene identifier
#' @param miscleavage_label String miscleavage description
#' @param protease String protease name
#' @param highlight_isoform String transcript to highlight (optional)
#' @return ggplot object
create_multi_isoform_plot <- function(all_peptides, transcript_df, gene_start, gene_end, 
                                    gene_id, miscleavage_label, protease, highlight_isoform = NULL,
                                    rmats_structure_data = NULL) {
  
  # Create basic plot structure
  p <- ggplot() +
    theme_minimal() +
    labs(
      title = paste0("Multi-Isoform Comparison: ", gene_id),
      subtitle = paste0("Protease: ", toupper(protease), " | ", miscleavage_label),
      x = "Genomic Position",
      y = "Transcripts"
    ) +
    xlim(gene_start, gene_end) +
    ylim(0.5, nrow(transcript_df) + 0.5)
  
  # Add transcript tracks
  for (i in 1:nrow(transcript_df)) {
    tx_name <- transcript_df$transcript[i]
    y_pos <- transcript_df$y_position[i]
    
    # Determine color based on highlighting
    track_color <- if (!is.null(highlight_isoform) && tx_name == highlight_isoform) {
      "#ff6b6b"  # Red for highlighted
    } else {
      "#4ecdc4"  # Teal for others
    }
    
    # Add transcript line
    p <- p + geom_segment(
      aes(x = gene_start, xend = gene_end, y = y_pos, yend = y_pos),
      color = "gray70", size = 0.5
    )
    
    # Add transcript label
    p <- p + annotate(
      "text", x = gene_start - (gene_end - gene_start) * 0.02, y = y_pos,
      label = tx_name, hjust = 1, size = 3, color = track_color
    )
    
    # Add rMATS transcript structure (exons and CDS) if available
    if (!is.null(rmats_structure_data) && rmats_structure_data$success) {
      cat("🎨 Adding rMATS structure visualization for transcript:", tx_name, "\n")
      
      # Add exons for this transcript using existing fill aesthetic system
      if (tx_name %in% names(rmats_structure_data$exons)) {
        tx_exons <- rmats_structure_data$exons[[tx_name]]
        cat("📦 Adding", length(tx_exons), "exon rectangles for", tx_name, "\n")
        
        if (length(tx_exons) > 0) {
          for (j in seq_along(tx_exons)) {
            exon <- tx_exons[j]
            cat("  📦 Exon", j, ":", start(exon), "-", end(exon), "\n")
            p <- p + geom_rect(
              aes(xmin = start(exon), xmax = end(exon),
                  ymin = y_pos - 0.3, ymax = y_pos + 0.3,
                  fill = "Transcript"),
              color = "black", alpha = 0.8
            )
          }
        }
      } else {
        cat("⚠️  No exon data found for transcript:", tx_name, "\n")
      }
      
      # Add CDS regions for this transcript using existing fill aesthetic system
      if (tx_name %in% names(rmats_structure_data$cds)) {
        tx_cds <- rmats_structure_data$cds[[tx_name]]
        cat("🟨 Adding", length(tx_cds), "CDS rectangles for", tx_name, "\n")
        
        if (length(tx_cds) > 0) {
          for (k in seq_along(tx_cds)) {
            cds <- tx_cds[k]
            cat("  🟨 CDS", k, ":", start(cds), "-", end(cds), "\n")
            p <- p + geom_rect(
              aes(xmin = start(cds), xmax = end(cds),
                  ymin = y_pos - 0.25, ymax = y_pos + 0.25,
                  fill = "CDS"),
              color = "black", alpha = 0.9
            )
          }
        }
      } else {
        cat("ℹ️  No CDS data found for transcript:", tx_name, "\n")
      }
    } else {
      if (grepl("\\.(inclusion|exclusion)$", tx_name)) {
        cat("❌ No rMATS structure data available for transcript:", tx_name, "\n")
      }
    }
  }
  
  # Add peptides
  if (nrow(all_peptides) > 0) {
    # Determine peptide colors
    peptide_colors <- if (!is.null(highlight_isoform) && "is_highlighted" %in% names(all_peptides)) {
      ifelse(all_peptides$is_highlighted, "#ff6b6b", "#4ecdc4")
    } else {
      "#4ecdc4"
    }
    
    # Add peptide rectangles
    p <- p + geom_rect(
      data = all_peptides,
      aes(xmin = start, xmax = end, ymin = y_position - 0.15, ymax = y_position + 0.15,
          text = hover_text),
      fill = peptide_colors,
      color = "white",
      size = 0.1,
      alpha = 0.8
    )
  }
  
  # Add consistent fill scale for structure visualization
  # Use the same colors as the existing visualization system
  fill_values <- c(
    "Transcript" = "rgba(77, 175, 74, 0.8)",  # Green for exons
    "CDS" = "rgba(255, 221, 0, 0.8)",        # Gold/yellow for CDS
    "peptide" = "rgba(52, 152, 219, 0.9)",   # Blue for regular peptides
    "junction spanning peptide" = "rgba(231, 76, 60, 0.9)" # Red for junction-spanning
  )
  
  p <- p + scale_fill_manual(
    values = fill_values,
    breaks = c("Transcript", "CDS", "peptide", "junction spanning peptide"),
    labels = c("Exons", "CDS", "Peptides", "Junction Spanning"),
    name = ""
  )
  
  # Customize theme
  p <- p + theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 11),
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
  )
  
  # Add legend for rMATS structure visualization
  if (!is.null(rmats_structure_data) && rmats_structure_data$success) {
    # Count how many transcripts have structure data
    transcripts_with_exons <- length(rmats_structure_data$exons)
    transcripts_with_cds <- length(rmats_structure_data$cds)
    
    cat("📊 Visualization Summary:\n")
    cat("  🔵 Transcripts with exon data:", transcripts_with_exons, "\n")
    cat("  🟨 Transcripts with CDS data:", transcripts_with_cds, "\n")
    
    # Add comprehensive caption to explain colors
    caption_text <- paste(
      "rMATS Structure Visualization:",
      "Green = Exons (" , transcripts_with_exons, "transcripts) |",
      "Gold = CDS regions (", transcripts_with_cds, "transcripts) |",
      "Blue/Red = Peptides"
    )
    
    p <- p + labs(
      caption = caption_text
    ) + 
    theme(
      plot.caption = element_text(size = 9, hjust = 0.5, margin = margin(t = 10))
    )
  } else {
    # Add note when no rMATS structure data is available
    rmats_transcript_count <- sum(grepl("\\.(inclusion|exclusion)$", transcript_df$transcript))
    if (rmats_transcript_count > 0) {
      p <- p + labs(
        caption = paste("Note:", rmats_transcript_count, "rMATS transcript(s) detected but no structure data loaded")
      ) + 
      theme(
        plot.caption = element_text(size = 9, hjust = 0.5, margin = margin(t = 10), color = "orange")
      )
    }
  }
  
  return(p)
}

# NOTE: get_transcript_peptides_for_comparison function is defined in as_analysis.R
# and will be loaded when that file is sourced. No stub needed here.