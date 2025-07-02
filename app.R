## ---- bootstrap packages (strict minimum) ---------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  # core app + data I/O
  "shiny", "readxl",

  # UI widgets
  "DT", "colourpicker", "rhandsontable",

  # analysis & plotting
  "pheatmap", "matrixStats", "scales"
)

install_if_missing <- function(pkgs) {
  new <- setdiff(pkgs, rownames(installed.packages()))
  if (length(new)) {
    install.packages(new, dependencies = c("Depends", "Imports"))
  }
}

install_if_missing(cran_pkgs)

## load quietly
invisible(lapply(cran_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))
options(shiny.maxRequestSize = 50 * 1024^2) # 50 MB

## --------------------------------------------------------------------------


make_anno <- function(df, cols, labels, id_col) {
  stopifnot(length(cols) == length(labels))
  stopifnot(id_col %in% names(df))
  out <- df[, cols, drop = FALSE]
  colnames(out) <- labels
  out <- data.frame(lapply(out, factor))
  rownames(out) <- df[[id_col]]
  out
}


ui <- fluidPage(
  titlePanel("GeoMx DSP Heatmap Plotter"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload GeoMx xlsx File", accept = ".xlsx"),
      uiOutput("group_select"),
      hr(),
      uiOutput("annot_col_ui"), 
      DTOutput("label_table"), 
      hr(),
      uiOutput("group_filter_ui"), 
      uiOutput("anno_color_ui"),
      hr(),
      uiOutput("filter_select"),
      uiOutput("filter_level_select"),
      checkboxInput("cluster_rows", "Cluster rows", value = TRUE),
      checkboxInput("cluster_cols", "Cluster columns", value = TRUE),
      actionButton("run", "Generate heatmap"),
      downloadButton("dlHeatmap", "Download PNG")
    ),

    mainPanel(
      plotOutput("heatmapPlot", height = "700px")
    )
  )
)

server <- function(input, output, session) {
  geoData <- reactive({
    req(input$file)
    list(
      props  = readxl::read_excel(input$file$datapath, sheet = "SegmentProperties"),
      counts = readxl::read_excel(input$file$datapath, sheet = "TargetCountMatrix")
    )
  })

  output$group_select <- renderUI({
    req(geoData())
    selectizeInput(
      inputId = "groupvar",
      label = "Select ROI identifier (type to search)",
      choices = names(geoData()$props),
      options = list(
        placeholder = "Start typing…",
        maxOptions  = 10 
      )
    )
  })

  output$annot_col_ui <- renderUI({
    req(geoData())
    selectizeInput(
      "annot_cols",
      "Choose annotation columns",
      choices  = names(geoData()$props),
      multiple = TRUE,
      options  = list(placeholder = "Start typing…", maxOptions = 5)
    )
  })

  ## ---------- table ----------------------------------------------
  labels_rv <- reactiveVal(data.frame()) 

  observeEvent(input$annot_cols, {
    labels_rv(data.frame(
      Column = input$annot_cols,
      Label = input$annot_cols,
      stringsAsFactors = FALSE
    ))
  })

  output$label_table <- renderDT(
    {
      req(labels_rv())
      datatable(
        labels_rv(),
        editable = list(
          target  = "cell",
          disable = list(columns = 0) 
        ),
        rownames = FALSE,
        selection = "none",
        options = list(dom = "t")
      )
    },
    server = FALSE
  )

  observeEvent(input$label_table_cell_edit, {
    info <- input$label_table_cell_edit 
    df <- labels_rv()
    df[info$row, "Label"] <- info$value 
    labels_rv(df)
  })

  anno_labels <- reactive(labels_rv()$Label) # Use edited labels

  anno_plot <- reactive({
    req(input$annot_cols, input$groupvar)
    make_anno(geoData()$props,
      cols   = input$annot_cols,
      labels = anno_labels(),
      id_col = input$groupvar
    )
  })

  output$filter_select <- renderUI({
    req(geoData())
    selectizeInput(
      inputId = "filtervar",
      label = "Select filtering variable (optional)",
      # prepend an empty choice named “None”
      choices = c("None" = "", names(geoData()$props)),
      options = list(
        placeholder      = "— none —",
        allowEmptyOption = TRUE,
        maxOptions       = 500
      )
    )
  })

  output$filter_level_select <- renderUI({
    req(nzchar(input$filtervar)) # skip if “None” / ""
    vals <- unique(geoData()$props[[input$filtervar]])
    #  vals <- vals[!is.na(vals)]                 # eliminate NA of selection.
    if (length(vals) >= 2) {
      selectInput("filterval", "Keep samples with", choices = vals)
    }
  })

  # # helper: turn arbitrary strings into safe IDs
  make_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

  observe({
    req(input$annot_cols)

    ui_list <- lapply(input$annot_cols, function(col) {
      map <- lvl_maps[[col]]
      if (is.null(map)) {
        return(NULL)
      } 

      keep_df <- map[as.logical(map$keep), , drop = FALSE]
      if (nrow(keep_df) == 0) {
        return(NULL)
      } 

      wellPanel(
        h4(col),
        lapply(seq_len(nrow(keep_df)), function(i) {
          new_lvl <- keep_df$new[i]

          ## carry over any colour the user already chose (old or new name)
          prev <- NULL
          for (nm in c(keep_df$orig[i], new_lvl)) {
            id_old <- paste0("col_", make_id(col), "_", make_id(nm))
            if (!is.null(input[[id_old]]) && nzchar(input[[id_old]])) {
              prev <- input[[id_old]]
              break
            }
          }

          colourpicker::colourInput(
            inputId = paste0("col_", make_id(col), "_", make_id(new_lvl)),
            label   = new_lvl,
            value   = prev %||% default_cols(nrow(keep_df))[i]
          )
        })
      )
    })

    output$anno_color_ui <- renderUI(tagList(ui_list))
  })

  output$group_filter_ui <- renderUI({
    req(input$annot_cols, geoData())
    props <- geoData()$props

    lapply(input$annot_cols, function(col) {
      lvls <- levels(factor(props[[col]]))
      rHandsontableOutput(paste0("lvltbl_", make_id(col)), height = 120 + 20 * length(lvls))
    })
  })

  lvl_maps <- reactiveValues() # named by annotation column

  make_lvl_df <- function(lvls) {
    data.frame(
      orig = lvls, # original (internal) level
      keep = rep(TRUE, length(lvls)), # checkbox
      new = lvls, # editable text
      stringsAsFactors = FALSE
    )
  }

  # build the tables once per annot column -----------------------------
  observeEvent(input$annot_cols, {
    req(geoData())
    props <- geoData()$props

    for (col in input$annot_cols) {
      lvls <- levels(factor(props[[col]]))
      lvl_maps[[col]] <- make_lvl_df(lvls)
    }
  })

  # render every table --------------------------------------------------
  observe({
    req(input$annot_cols)
    for (col in input$annot_cols) {
      local({
        cc <- col
        tbl <- paste0("lvltbl_", make_id(cc))

        output[[tbl]] <- renderRHandsontable({
          rhandsontable(lvl_maps[[cc]], rowHeaders = NULL) %>%
            hot_col("keep", type = "checkbox") %>% # show as tick box
            hot_col("orig", readOnly = TRUE) %>% # lock original name
            hot_col("orig", title = "Level") %>%
            hot_col("keep", title = "Keep?") %>%
            hot_col("new", title = "Display name")
        })
      })
    }
  })

  # capture edits -------------------------------------------------------
  observe({
    req(input$annot_cols)
    for (col in input$annot_cols) {
      tbl_id <- paste0("lvltbl_", make_id(col))
      if (!is.null(input[[tbl_id]])) {
        lvl_maps[[col]] <- hot_to_r(input[[tbl_id]])
      }
    }
  })



  results <- eventReactive(input$run, {
    req(input$groupvar)
    validate(need(
      input$groupvar %in% names(geoData()$props),
      "'ROI identifier' you selected is not a column in SegmentProperties."
    ))

    data <- geoData()
    props <- data$props
    counts <- data$counts

    ## ─────────────────────────────
    ## 0. drop unselected groups
    ## ─────────────────────────────
    if (length(input$annot_cols)) {
      # drop groups -----------------
      for (col in input$annot_cols) {
        map <- lvl_maps[[col]]
        keep_lvls <- map$orig[as.logical(map$keep)]
        props <- props[props[[col]] %in% keep_lvls, , drop = FALSE]
      }
      # rename levels ---------------
      for (col in input$annot_cols) {
        map <- lvl_maps[[col]]
        ren <- setNames(map$new, map$orig)
        props[[col]] <- ren[props[[col]]]
      }
    }

    ## ─────────────────────────────
    ## 1.  Identify matched samples
    ## ─────────────────────────────
    id_props <- props[[input$groupvar]] 
    sample_cols <- intersect(names(counts)[-1], id_props)

    validate(need(
      length(sample_cols) > 1,
      "TargetCountMatrix has no columns matching the identifier."
    ))

    ## ─────────────────────────────
    ## 2.  Expression matrix
    ## ─────────────────────────────
    expr_mat <- as.matrix(counts[, sample_cols])
    rownames(expr_mat) <- counts$TargetName
    expr_mat <- expr_mat[!grepl("Negative Probe", rownames(expr_mat)), ]
    storage.mode(expr_mat) <- "numeric"

    ## ─────────────────────────────
    ## 3.  Phenotype data-frame
    ## ─────────────────────────────
    pheno <- props[match(sample_cols, id_props), ]
    pheno <- as.data.frame(pheno) 
    rownames(pheno) <- pheno[[input$groupvar]] 
    stopifnot(identical(colnames(expr_mat), rownames(pheno)))

    ## ─────────────────────────────
    ## 3b.  OPTIONAL: filter samples
    ## ─────────────────────────────
    if (nzchar(input$filtervar) && nzchar(input$filterval)) {
      keep_filter <- pheno[[input$filtervar]] == input$filterval
      validate(need(
        sum(keep_filter) > 1,
        "No samples remain after applying the filter."
      ))
      pheno <- pheno[keep_filter, , drop = FALSE]
      expr_mat <- expr_mat[, keep_filter, drop = FALSE]
      stopifnot(identical(colnames(expr_mat), rownames(pheno)))
    }
    ## ─────────────────────────────
    ## 4.  Drop samples with missing vars
    ## ─────────────────────────────
    # build the full set of columns that must be present
    vars_needed <- unique(c(
      input$groupvar,
      input$blockvar,
      input$annot_cols
    ))

    # drop empty strings (in case blockvar is blank)
    vars_needed <- vars_needed[nzchar(vars_needed)]

    # sanity-check: all columns exist in pheno
    missing_vars <- setdiff(vars_needed, names(pheno))
    validate(need(
      length(missing_vars) == 0,
      paste(
        "Variable(s) not found in SegmentProperties:",
        paste(missing_vars, collapse = ", ")
      )
    ))

    # keep rows where *every* requested column is not NA
    keep <- stats::complete.cases(pheno[vars_needed])

    pheno <- pheno[keep, , drop = FALSE]
    expr_mat <- expr_mat[, keep, drop = FALSE]

    stopifnot(identical(colnames(expr_mat), rownames(pheno)))

    list(
      expr_mat = expr_mat,
      pheno = pheno
    )
  })

  clean_pheno <- reactive({
    req(results(), input$annot_cols)
    ph <- results()$pheno

    ph[input$annot_cols] <- lapply(
      ph[input$annot_cols],
      function(x) droplevels(factor(x))
    )
    ph
  })


  ## ─────────────────────────────
  ## 5.  create the clean data and color objects
  ##
  ## ─────────────────────────────
  anno_plot <- reactive({
    req(input$annot_cols, clean_pheno())
    make_anno(
      df     = clean_pheno(), 
      cols   = input$annot_cols,
      labels = anno_labels(),
      id_col = input$groupvar
    )
  })


  `%||%` <- function(a, b) if (is.null(a) || a == "") b else a
  default_cols <- function(n) scales::hue_pal()(n)


  anno_colors <- reactive({
    req(input$annot_cols)
    out <- lapply(input$annot_cols, function(col) {
      map <- lvl_maps[[col]]
      if (is.null(map)) {
        return(NULL)
      }

      keep_df <- map[as.logical(map$keep), , drop = FALSE]
      lvls <- keep_df$new

      setNames(
        vapply(seq_along(lvls), function(i) {
          id <- paste0("col_", make_id(col), "_", make_id(lvls[i]))
          input[[id]] %||% default_cols(length(lvls))[i]
        }, character(1)),
        lvls
      )
    })
    names(out) <- anno_labels() 
    out
  })

  ## ─────────────────────────────
  ## 6.  Heat-map object
  ##    (fires only when Run is pressed)
  ## ─────────────────────────────
  heat_obj <- eventReactive(input$run, {
    ## 6a.  grab the cleaned matrix & annotation -------------------------------
    expr_mat <- results()$expr_mat 
    anno_df <- anno_plot() 

    ## 6b.  log-transform and pick high-CV genes -------------------------------
    log_mat <- log2(expr_mat + 1)

    cv <- matrixStats::rowSds(log_mat) /
      matrixStats::rowMeans2(log_mat)

    cut_cv <- stats::quantile(cv, 0.80, na.rm = TRUE) # 80-th percentile
    keep_g <- which(cv > cut_cv)

    mat_hcv <- log_mat[keep_g, , drop = FALSE]
    mat_hcv <- mat_hcv[order(cv[keep_g], decreasing = TRUE), ]

    ## 6c.  draw pheatmap and *return* the object ------------------------------
    pheatmap::pheatmap(
      mat_hcv,
      scale = "row",
      show_rownames = FALSE,
      show_colnames = FALSE,
      clustering_method = "average",
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      cluster_rows = input$cluster_rows,
      cluster_cols = input$cluster_cols,
      annotation_col = anno_df,
      annotation_colors = anno_colors(),
      border_color = NA,
      breaks = seq(-3, 3, length.out = 121),
      color = grDevices::colorRampPalette(
        c("purple3", "black", "yellow2")
      )(120)
    )
  })

  ## ─────────────────────────────
  ## 7.  Render in the main panel
  ## ─────────────────────────────
  output$heatmapPlot <- renderPlot({
    heat_obj()
  })

  ## ─────────────────────────────
  ## 8.  Download handler
  ## ─────────────────────────────
  output$dlHeatmap <- downloadHandler(
    filename = function() paste0("GeoMx_heatmap_", Sys.Date(), ".png"),
    content = function(file) {
      png(file, width = 2400, height = 1600, res = 200)
      grid::grid.newpage()
      grid::grid.draw(heat_obj()$gtable) # explicitly draw the plot grob
      dev.off()
    }
  )
}

shinyApp(ui, server, options = list(host = "0.0.0.0", port = 3839))
