# Author: Etienne CAMENEN
# Date: 2019
# Contact: arthur.tenenhaus@l2s.centralesupelec.fr
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.)
# and produces textual and graphical outputs (e.g. variables and individuals
# plots).
    
server <- function(input, output, session) {
    ################################################ Render UI ################################################

    
    hide(selector = "#tabset li a[data-value=RGCCA]")
    hide(id = "b_x_custom")
    hide(id = "b_y_custom")
    
    output$tau_custom <- renderUI({
        refresh <- c(input$superblock)
        isolate (setAnalysis())
        setTauUI()
    })

    output$nb_mark_custom <- renderUI({
        refresh <- c(input$blocks_names_custom_x, input$blocks_names_custom_x)
        sliderInput(
            inputId = "nb_mark",
            label = "Number of top variables",
            min = 10,
            max = getMaxCol(),
            value = getDefaultCol(),
            step = 1
        )
    })

    output$connection_custom <- renderUI({
        setUiConnection()
    })

    output$response_custom <- renderUI({
        setUiResponse()
    })

    output$blocks_names_custom_x <- renderUI({
        setNamesInput("x", bool = input$navbar == "Samples")
    })

    output$blocks_names_custom_y <- renderUI({
        setNamesInput("y")
    })

    output$blocks_names_response <- renderUI({
        setNamesInput("response", "Block used as a response")
    })

    output$nb_compcustom <- renderUI({
        refresh <- c(input$superblock, input$each_ncomp)
        isolate (setAnalysis())
        setCompUI()
    })

    refreshAnalysis <- function()
        c(
            input$nb_comp,
            input$block,
            input$sep,
            input$scheme,
            input$scale,
            input$superblock,
            input$supervised
        )

    output$compx_custom <- renderUI({
        refresh <- refreshAnalysis()
        refresh <- input$names_block_x
        isolate(uiComp("x", 1, id_block, input$navbar != "Fingerprint"))
    })

    output$compy_custom <- renderUI({
        refresh <- refreshAnalysis()
        refresh <- input$names_block_y
        uiComp("y", min(getNcomp()), id_block_y)
    })

    output$analysis_type_custom <- renderUI({
        refresh <- c(input$blocks, input$sep)
        selectInput(
            inputId = "analysis_type",
            "Analysis method",
            selected = analysis_type,
            choices = list(
                `One block` = one_block,
                `Two blocks` = two_blocks,
                `Multiblocks` = multiple_blocks,
                `Multiblocks with a superblock` = multiple_blocks_super
            )
        )
    })
    
    output$b_x_custom <- renderUI(b_index("x", "estimate"))
    output$b_y_custom <- renderUI({
            b_index("y", 
                if (tolower(input$analysis_type) == "sgcca")
                    "occurrences"
                else
                    "sign")
        })

    b_index <- function(x, y) {
        
        l_choices <- list(
            `RGCCA weights` = "estimate",
            `Bootstrap-ratio` = "bootstrap_ratio",
            `Mean bootstrap weights` = "mean"
        )
        
        if (tolower(input$analysis_type) == "sgcca")
            l_choices[["Non-zero occurrences"]] <- "occurrences"
        else
            l_choices[["Significant 95% interval"]] <- "sign"

        selectInput(
            inputId = paste0("b_", x),
            paste("Bootstrap indexes for", x, "axis"),
            choices = l_choices,
            selected = y
        )
            
    }

    ################################################ UI function ################################################

    setTauUI <- function(superblock = NULL) {
        refresh <- c(input$superblock, input$supervised, input$tau_opt)

        if (!is.null(input$analysis_type) &&
            input$analysis_type == "SGCCA") {
            par_name <- "Sparsity"
            cond <- "input.tau_opt == false && input.analysis_type == SGCCA"
        } else{
            par_name <- "Tau"
            cond <- "input.tau_opt == false"
        }

        conditionalPanel(condition = cond,
                        lapply(1:(length(blocks)), function(i) {
                            sliderInput(
                                inputId = paste0("tau", i),
                                label = paste(par_name, "for", names(getNames())[i]),
                                min = ifelse(
                                    par_name == "Tau",
                                    0,
                                    ceiling(1 / sqrt(NCOL(blocks[[i]])) * 100) / 100),
                                max = 1,
                                value = ifelse(
                                    is.null(input[[paste0("tau", i)]]),
                                    1, 
                                    input[[paste0("tau", i)]]),
                                step = .01
                            )
                        }))
    }

    setCompUI <- function(superblock = NULL) {
        if (!input$each_ncomp)
            sliderInput(
                inputId = "ncomp",
                label = "Number of components",
                min = 1,
                max = min(getMaxComp()),
                value = 2, #TODO: test with univariate
                step = 1
            )
        else
            lapply(1:(length(blocks)), function(i) {
                sliderInput(
                    inputId = paste0("ncomp", i),
                    label = paste("Number of components for", names(getNames())[i]),
                    min = 1,
                    max = getMaxComp()[i],
                    value = 2,
                    step = 1
                )
            })
    }

    setNamesInput <- function(x, label = NULL, bool = TRUE) {
        refesh <- c(input$superblock, input$supervised, input$analysis_type)

        if (is.null(label)) {
            label <- "Block"

            if (bool)
                label <- paste0("Block for the ", x , "-axis")

        }

        selectInput(
            inputId = paste0("names_block_", x),
            label = label,
            choices = getNames(),
            selected = setBlockNames()
        )
    }

    # Define the names of the blocks and set by default on the last block
    setBlockNames <- function() {
        if (!is.null(input$blocks)) {
            if (!is.null(id_block))
                return(id_block)
            else
                return(round(length(blocks)))
            # Set selected value on the last block

        } else{
            # If any dataset is selected
            return(1)
        }
    }

    uiComp <- function(x, y, id_block, bool = TRUE) {
        label <- "Component"

        if (bool)
            label <- paste0("Component for the ", x, "-axis")

         comp <- getNcomp()
        if (length(comp) > 1)
            comp <- comp[id_block]

        sliderInput(
            inputId = paste0("comp", x),
            label = label,
            min = 1,
            max = comp,
            value = y,
            step = 1
        )
    }

    output$file_custom <- renderUI({
        ui <- fileInput(inputId = "blocks",
                        label = "Blocks",
                        multiple = TRUE)

        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(title = "One or multiple CSV files
                    containing a matrix with : (i) quantitative values only
                    (decimal should be separated by '.'), (ii) the samples in
                    lines (should be labelled in the 1rst column) and (iii)
                    variables in columns (should have a header)")
            )

        return(ui)
    })

    output$sep_custom <- renderUI({
        ui <- radioButtons(
            inputId = "sep",
            label = "Column separator",
            choices = c(
                Comma = ",",
                Semicolon = ";",
                Tabulation = "\t"
            ),
            selected = "\t"
        )

        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(title = "Character used to separate the
                    column in the dataset")
            )

        return(ui)
    })

    output$scale_custom <- renderUI({
        ui <- checkboxInput(inputId = "scale",
                            label = "Scale the blocks",
                            value = TRUE)

        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(title = "A data centering step is always
                    performed. If ticked, each block is normalised and divided
                    by the square root of its number of variables.")
            )

        return(ui)
    })

    
    output$nperm_custom <- renderUI({
        ui <- sliderInput(
            inputId = "nperm",
            label = "Number of permutations",
            min = 5,
            max = 1000,
            value = 10,
            step = 5
        )
        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(title = "To tune the sparsity coefficient (if the model is sparse) or tau
                                     (otherwise), we observe the deviation between the model and a set of models
                                     where the lines of each block are permuted. The model with the best 
                                     combination of parameters is the one with the highest deviation.")
            )
        
        return(ui)
    })

    output$nboot_custom <- renderUI({
        ui <- sliderInput(
            inputId = "nboot",
            label = "Number of boostraps",
            min = 5,
            max = 1000,
            value = 10,
            step = 5
        )
        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(title = "By taking several random samples from the dataset (bootstrap),
                                     the importance of the variables may vary. The variables that are most 
                                     often selected are those that are retained.")
            )
        
        return(ui)
    })

    output$val_custom <- renderUI({
        ui <- radioButtons(
            "val",
            label = "Type of validation",
            choices = c(#`Train-test` = "test",
                `K-fold` = "kfold",
                `Leave-one-out` = "loo"),
            selected = "loo"
        )
        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(title = "To tune the sparsity coefficient (if the model is sparse) or
                                     tau (otherwise), in supervised mode, we observe the performance (RMSE)
                                     of a model from which individuals were randomly drawn. These individuals
                                     can be divided into k folds where the model will be tested on each fold
                                     and trained on the others. For small datasets (<30 samples), it is 
                                     recommended to use as many folds as there are individuals (leave-one-out; 
                                     loo). The best combination of parameters is the one where, on average, 
                                     the samples perform best.")
            )
        
        return(ui)
    })
    
    output$tau_opt_custom <- renderUI({

        if (!is.null(input$analysis_type) && input$analysis_type == "SGCCA") {
            penalty <- "sparsity"
            text <- "A sparsity coefficient varies from the square root of the variable number (the fewest selected variables) to 1 (all the variables are included)"
        } else if (!is.null(input$analysis_type) && input$analysis_type == "RGCCA") {
            penalty <- "tau"
            text <- "A tau near 0 maximize the the correlation whereas a tau near 1 maximize the covariance"
        } else
            penalty <- text <- ""

        ui <- checkboxInput(inputId = "tau_opt",
                            label = paste("Use an optimal", penalty),
                            value = FALSE)

        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(title = text)
            )

        return(ui)
    })

    output$scheme_custom <- renderUI({
        ui <- radioButtons(
            inputId = "scheme",
            label = "Scheme function",
            choices = c(
                Horst = "horst",
                Centroid = "centroid",
                Factorial = "factorial"
            ),
            selected = "factorial"
        )

        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(
                        title =
                            "Link (i.e. scheme) function for covariance
                            maximization is calculated with: the identity
                            function (horst scheme), the absolute values
                            (centroid scheme), the squared values (factorial
                            scheme). Only, the horst scheme penalizes structural
                            negative correlation. The factorial scheme
                            discriminates more strongly the blocks than the
                            centroid one."
                    )
                    )

        return(ui)
    })


    output$superblock_custom <- renderUI({
        ui <- checkboxInput(inputId = "superblock",
                            label = "Use a superblock",
                            value = TRUE)

        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(
                        title =
                            "If ticked, a superblock is introduced. This
                            superblock is defined as a concatenation of all the
                            other blocks. The space spanned by global components
                            is viewed as a compromise space that integrated
                            all the modalities and facilitates the
                            visualization of the results and their
                            interpretation. If unchecked, a connection file
                            could be used. Otherwise, all blocks are assumed
                            to be connected."
                    )
                    )


        return(ui)
    })

    setUiConnection <- function() {
        refresh <- c(input$connection)

        ui <- fileInput(inputId = "connection",
                        label = "Connection design [OPTIONAL]")

        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(title = "The design matrix is a symmetric
                    matrix of the length of the number of blocks describing
                    the connections between them. Two values are accepted :
                    '1' for a connection between two blocks, or '0' otherwise. 
                    By default, all the blocks are connected together.")
            )

        conditionalPanel(
            condition = "!input.superblock && !input.supervised",
            ui)
    }

    setUiResponse <- function() {
        refresh <- c(input$response)

        ui <- fileInput(inputId = "response",
                        label = "Color with a response [OPTIONAL]")
        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(title = "To color the sample plot.
                    A CSV file containing either : (i) an only column with a
                    qualitative or a quantitative variable; (ii) multiple
                    columns corresponding to a disjunctive table")
            )

        return(ui)
    }

    getMaxComp <- function(){
        comp <- sapply(
            blocks,
            function(x) {
                comp <- NCOL(x)
                if (comp > 5)
                    return(5)
                else
                    return(comp)
            }
        )
    }
   
    getNames <- function() {
        # Get the names of the blocks

        if (!is.null(input$blocks)) {
            # Creates a list of nb_blocks dimension, each one containing a id
            # from 1 to nb_blocks and having the same names as the blocks
            return(as.list(
                sapply(names(blocks), function(i)
                    as.integer(which(names(blocks) == i)), USE.NAMES = TRUE)
            ))
        } else
            return(list(" " = 1))

    }

    getMaxCol <- function() {
        # Get the maximum number of columns among the blocks

        if (!is.null(input$blocks)) {
            return(NCOL(blocks[[id_block]]))
        } else
            return(100)

    }

    getDefaultCol <- function() {
        # Set the maximum of biomarkers to the maximum
        # number of column among the blocks but not lower than 100

        max <- getMaxCol()

        if (max < 50)
            return (max)
        else
            return (50)
    }

    showWarn <- function(f, duration = 10, show = TRUE, warn = TRUE) {

        ids <- character(0)

        try(withCallingHandlers({
            res <- f
        }, message = function(m) {
            if (show)
                duration <<- NULL

            id <- showNotification(
                m$message,
                type = "message", 
                duration = duration)
            ids <<- c(ids, id)

        }, warning = function(w) {
            warning(w$message)
            if (show && warn) {
                id <- showNotification(
                    w$message,
                    type = "warning",
                    duration = duration)
                ids <<- c(ids, id)
            }

        }, error = function(e) {
            message(paste("Error:", e$message))
            id <- showNotification(
                e$message,
                type = "error",
                duration = duration)
            ids <<- c(ids, id)
            res <<- class(e)[1]
        }),
        silent = TRUE)

        if (is.null(duration) & length(ids) != 0) {
            for (id in ids)
                removeNotification(id)
        }

        return(res)
    }

    blocksExists <- function() {
        # Test if the blocks are loaded and contain any errors

        if (!is.null(input$blocks))
            if (!is.null(getInfile()))
                return(TRUE)
        return(FALSE)
    }

    setAnalysisMenu <- function() {
        refresh <- c(input$blocks, input$analysis_type)
        assign("one_block", analyse_methods[[1]], .GlobalEnv)
        assign("two_blocks", analyse_methods[[2]], .GlobalEnv)
        assign("multiple_blocks", analyse_methods[[3]], .GlobalEnv)
        assign("multiple_blocks_super", analyse_methods[[4]], .GlobalEnv)
    }

    setIdBlock <- function() {
        assign("id_block", length(blocks), .GlobalEnv)
        assign("id_block_y", length(blocks), .GlobalEnv)

    }

    getDynamicVariables <- function() {
        # Refresh all the plots when any input is changed

        refresh <- c(
            input$sep,
            input$header,
            input$blocks,
            input$superblock,
            input$connection,
            input$scheme,
            input$nb_mark,
            input$scale,
            input$init,
            input$compx,
            input$compy,
            input$tau,
            input$tau_opt,
            input$analysis_type,
            input$connection,
            input$nb_comp,
            input$response,
            input$names_block_x,
            input$names_block_y,
            input$nboot,
            input$nperm,
            input$run_perm,
            input$run_crossval,
            input$kfold,
            input$show_crossval,
            input$text,
            input$names_block_response,
            input$supervised,
            input$run_analysis,
            input$run_boot,
            input$nb_mark_custom,
            input$blocks_names_custom_x
        )
    }


    ################################################ Plots  ################################################

    getExtension <- function(f) {
        if (!is.null(f)) {
            format <- unlist(strsplit(f, '.', fixed = "TRUE"))
            return(paste(format[-length(format)], collapse = "."))
        } else
            return(f)
    }

    samples <- function() {

        if (!input$show_crossval)
            crossval <- NULL

        if (!is.null(crossval)) {
            response_name <- "Response"
            response <- rep(1, NROW(rgcca_out$Y[[1]]))
        } else if (!is.null(input$response))
            response_name <- getExtension(input$response$name)
        else
            response_name <- ""

        if (!is.null(input$compx))
            plot_ind(
                rgcca = rgcca_out,
                resp = response,
                compx = input$compx,
                compy = input$compy,
                i_block = id_block,
                text = if_text,
                i_block_y = id_block_y,
                response_name = response_name,
                predicted = crossval
            )
    }

    corcircle <- function()
        if (!is.null(input$compx))
            plot_var_2D(
                rgcca = rgcca_out,
                compx = compx,
                compy = compy,
                i_block = id_block,
                text = if_text,
                n_mark = nb_mark
            )

    fingerprint <- function(type)
        plot_var_1D(
            rgcca = rgcca_out,
            comp = compx,
            n_mark = nb_mark,
            i_block = id_block,
            type = type
        )

    ave <- function()
        plot_ave(rgcca = rgcca_out)

    design <- function()
        plot_network2(rgcca_out)
    
    design2 <- function()
        plot_network(rgcca_out)

    plotBoot <- function(){
        refresh <- c(input$names_block_x, id_block, input$blocks_names_custom_x)
        plot_bootstrap_1D(
            df_b = selected.var,
            # x = input$b_x,
            # y = input$b_y,
            n_mark = nb_mark
        )
    }

    viewPerm <- function(){
        head(
            order_df(cbind(perm$penalties, `Z-score` = perm$zstat), ncol(perm$penalties) + 1), 
            n = nb_mark)
    }

    ################################################ Analysis ################################################

    getTau <- function() {
        tau <- integer(0)
        for (i in 1:(length(blocks_without_superb) + ifelse(input$superblock, 1, 0)))
            tau <- c(tau, input[[paste0("tau", i)]])

        return(tau)
    }

    getNcomp <- function() {
        if (input$each_ncomp) {
            ncomp <- integer(0)
            cond <- input$superblock && ( toupper(analysis_type) %in% c("PCA", "RGCCA", "SGCCA") ||
                    analysis_type %in% multiple_blocks_super)
            for (i in 1:(length(blocks_without_superb) + ifelse(cond, 1, 0)))
                ncomp <- c(ncomp, input[[paste0("ncomp", i)]])
        } else
            ncomp <- input$ncomp

        return(ncomp)
    }
    
    setParRGCCA <- function(verbose = TRUE) {
        blocks <- blocks_without_superb

        if (is.null(analysis_type) | is.null(input$analysis_type))
            analysis_type <- "RGCCA"
        else
            analysis_type <- input$analysis_type

        # Tau is set to 1 by default
        if (is.null(input$tau_opt))
             tau <- 1
        else{
            # otherwise the tau value fixed by the user is used
            tau <- getTau()
        }

        setAnalysisMenu()

        if (length(blocks) == 1) {
            # if(verbose)
            # showWarn(warning("Only one block is selected. By default, a 
            # PCA is performed."))
            analysis_type <- "PCA"
            assign("two_blocks", NULL, .GlobalEnv)
            assign("multiple_blocks", NULL, .GlobalEnv)
            assign("multiple_blocks_super", NULL, .GlobalEnv)
        } else if (length(blocks) == 2) {
            assign("one_block", NULL, .GlobalEnv)
            assign("multiple_blocks", NULL, .GlobalEnv)
            assign("multiple_blocks_super", NULL, .GlobalEnv)
            if (!tolower(analysis_type) %in% c("cca", "ra", "ifa", "pls")) {
                # showWarn(warning("Only two blocks are selected. By default, 
                # a PLS is performed."))
                analysis_type <- "PLS"
            }
        } else if (length(blocks) > 2) {
            assign("one_block", NULL, .GlobalEnv)
            assign("two_blocks", NULL, .GlobalEnv)
        }

        if ( (!is.null(input$superblock) && input$superblock) && 
                ( toupper(analysis_type) %in% c("PCA", "RGCCA", "SGCCA")) ||
                analysis_type %in% multiple_blocks_super ) {
            blocks <- c(blocks, superblock = list(Reduce(cbind, blocks)))
        }

        assign("tau", tau, .GlobalEnv)
        assign("analysis_type", analysis_type, .GlobalEnv)
        
        if (analysis_type != "SGCCA")
             assign("perm.par", "tau", .GlobalEnv)
        else
            assign("perm.par", "sparsity", .GlobalEnv)

        return(blocks)
    }

    setRGCCA <- function() {
        # Load the analysis

        isolate({
            if (!is.null(cv))
                tau <- cv$bestpenalties
            else if (!is.null(perm))
                tau <- perm$bestpenalties
            else if (length(grep("[SR]GCCA", analysis_type)) == 1)
                tau <- getTau()
        })

        if (!is.null(input$supervised) && input$supervised)
            response <- input$names_block_response
        else
            response <- NULL

        # scheme_power <- input$power
        # if (input$scheme == "factorial")
        #     scheme <- function (x) x^as.integer(scheme_power)
        # else
            scheme <- input$scheme

        assign("rgcca_out",
               showWarn({
                   func <- quote(
                       rgcca(
                           blocks_without_superb,
                            connection = connection,
                            response = response,
                            superblock = (!is.null(input$supervised) &&
                                !is.null(input$superblock) && input$superblock),
                            ncomp = getNcomp(),
                            scheme = scheme,
                            scale = FALSE,
                            scale_block = FALSE,
                            init = input$init,
                            bias = TRUE,
                            type = analysis_type
                        )
                   )
                   if (tolower(analysis_type) %in% c("sgcca", "spca", "spls"))
                       func[["sparsity"]] <- tau
                   else
                       func[["tau"]] <- tau
                   eval(as.call(func))
                }),
                .GlobalEnv)

    }

    getCrossVal <- function(){

            isolate({
            if (length(grep("[SR]GCCA", analysis_type)) == 1)
                tau <- getTau()
        })

        if (!is.null(input$supervised) && input$supervised)
            response <- input$names_block_response
        else
            response <- NULL

        assign(
            "cv", {
            func <- quote(
                rgcca_cv( 
                    blocks,
                    type = analysis_type,
                    response = response,
                    validation = input$val,
                    k = input$kfold,
                    n_cv = input$ncv,
                    n_cores = parallel::detectCores() - 1,
                    superblock = (!is.null(input$supervised) &&
                                    !is.null(input$superblock) && input$superblock),
                    scale = FALSE,
                    scale_block = FALSE,
                    scheme = input$scheme,
                    parallelization = TRUE,
                    init = input$init,
                    new_scaled = TRUE,
                    ncomp = getNcomp()))
                if (tolower(analysis_type) %in% c("sgcca", "spca", "spls")) {
                    func[["sparsity"]] <- tau
                    func[["par"]] <- "sparsity"
                } else {
                    func[["tau"]] <- tau
                    func[["par"]] <- "tau"
                }
                showWarn(eval(as.call(func)))
            },
            .GlobalEnv
        )

        show(id = "navbar")
        show(id = "run_analysis")
        show(selector = "#navbar li a[data-value=Cross-validation]")
        updateTabsetPanel(session, "navbar", selected = "Cross-validation")
    }

    getCrossVal2 <-  function(){
        assign(
            "crossval",
            rgcca_crossvalidation(rgcca_out, validation = input$val, k = input$kfold, n_cores = 1),
            .GlobalEnv
        )
        showWarn(message(paste("CV score:", round(crossval$score, 4))), show = FALSE)
        updateTabsetPanel(session, "navbar", selected = "Samples")
    }

    getPerm <-  function(){
        isolate({
            if (length(grep("[SR]GCCA", analysis_type)) == 1)
                tau <- getTau()
        })

        if (!is.null(input$supervised) && input$supervised)
            response <- input$supervised
        else
            response <- NULL

        assign("perm", {
            func <- quote(
                rgcca_permutation(
                    blocks_without_superb,
                    perm.par = perm.par,
                    nperm = input$nperm,
                    connection = connection,
                    response = input$names_block_response,
                    superblock = (!is.null(input$supervised) &&
                                      !is.null(input$superblock) && input$superblock),
                    scheme = input$scheme,
                    scale = FALSE,
                    scale_block = FALSE,
                    ncomp = getNcomp(),
                    init = input$init,
                    bias = TRUE,
                    type = analysis_type
                )
            )
            if (tolower(analysis_type) %in% c("sgcca", "spca", "spls"))
                func[["sparsity"]] <- tau
            else
                func[["tau"]] <- tau
            showWarn(eval(as.call(func)))
        },
        .GlobalEnv)
 
        show(id = "navbar")
        show(id= "run_analysis")
        show(selector = "#navbar li a[data-value=Permutation]")
        show(selector = "#navbar li a[data-value='Permutation Summary']")
        updateTabsetPanel(session, "navbar", selected = "Permutation")
    }

    getBoot <-  function(){
        assign(
            "boot",
            showWarn(bootstrap(rgcca_out, n_boot = input$nboot)),
            .GlobalEnv
        )
        assign("selected.var", NULL, .GlobalEnv)
        show(selector = "#navbar li a[data-value=Bootstrap]")
        show(selector = "#navbar li a[data-value='Bootstrap Summary']")
    }

    load_responseShiny = function() {
        response <- showWarn(
            load_response(
                blocks = blocks_without_superb,
                file = response_file,
                sep = input$sep,
                header = input$header
            )
        )

        if (length(response) < 1)
            response <- NULL -> response_file

        return(response)

    }

    set_connectionShiny <- function() {
        supervised <- (!is.null(input$supervised) && input$supervised)

        if (!is.null(connection_file)) {
            connection <- load_connection(file = connection_file, sep = input$sep)

            check <- showWarn(check_connection(connection, blocks))

            # Error due to the superblock disabling and the connection have not the same size than the number of blocks
            if (identical(check, "130")) 
                connection <- NULL

        }

        if (is.matrix(connection)) {
            assign("connection", connection, .GlobalEnv)
            cleanup_analysis_par()
        }

    }

    setAnalysis <- function() {
        blocks <- setParRGCCA()

        if (!is.null(blocks)) {
           cleanup_analysis_par()
            assign("blocks", blocks, .GlobalEnv)
            set_connectionShiny()
            setIdBlock()
        }

    }

    ################################################ Events ################################################


    setToggle <- function(id)
            toggle(
                condition = (
                    input$analysis_type %in% c("RGCCA", "SGCCA") 
                    && length(blocks) > 2
            ),
            id = id)


    setToggle2 <- function(id)
            toggle(
                condition = (input$analysis_type %in% c("RA", "RGCCA", "SGCCA")),
                   id = id)


    setToggleSaveButton <- function(id)
            toggle(condition = !is.null(analysis), id = id)


    observe({
        # Event related to input$analysis_type
        for (i in c("tau_custom", "tau_opt", "scheme", "superblock", "connection", "supervised" ))
            setToggle(i)
        setToggle2("blocks_names_response")
        hide(selector = "#tabset li a[data-value=Graphic]")
        toggle(
            condition = (length(input$blocks$datapath) > 1), 
            id = "blocks_names_custom_x")
        toggle(
            condition = (length(input$blocks$datapath) > 1), 
            id = "blocks_names_custom_y")
    })

    # observeEvent(c(input$compx, input$compy), {
    #     for (i in c("Corcircle", "Samples"))
    #         toggle(condition = !is.null(analysis) && !(input$compx < 2 && input$compy < 2),
    #                selector = paste0("#navbar li a[data-value=", i, "]"))
    # })

    observeEvent(c(input$names_block_x), {
        comp <- getNcomp()
        if (length(comp) > 1)
            comp <- comp[id_block]
        condition <- !is.null(analysis) && comp > 1
        if (!condition)
            updateTabsetPanel(session, "navbar", selected = "Samples")
        for (i in c("Corcircle", "Fingerprint"))
            toggle(condition = condition, 
                   selector = paste0("#navbar li a[data-value=", i, "]"))
    })

    observeEvent(c(input$navbar, input$tabset), {
        toggle(
            condition = (input$navbar %in% c("Corcircle", "Fingerprint", "Bootstrap")),
               id = "nb_mark_custom")
        for (i in c("text", "compy_custom"))
            toggle(
                condition = ( !input$navbar %in% c("Fingerprint", "Bootstrap")),
                   id = i)
        toggle(
            condition = (input$navbar == "Samples" && 
                    length(input$blocks$datapath) > 1),
               id = "blocks_names_custom_y")
        toggle(condition = input$navbar == "Samples", id = "response_custom")
        toggle(condition = input$navbar == "Samples" && !is.null(crossval), id = "show_crossval")
        toggle(
            condition = (input$navbar == "Fingerprint"),
            id = "indexes")
        # for (i in c("b_x_custom", "b_y_custom"))
        #     toggle(condition = (input$navbar == "Bootstrap"), id = i)
        toggle(
            condition = (
                !is.null(analysis) && !input$navbar %in% c("Connection", "AVE", "Cross-validation", "'Bootstrap Summary'", "Permutation", "'Permutation Summary'")
            ),
            selector = "#tabset li a[data-value=Graphic]"
        )
    })


    observeEvent(input$navbar, {
        if (!is.null(analysis) && input$navbar %in% c("Connection", "AVE", "Permutation", "Cross-validation"))
            updateTabsetPanel(session, "tabset", selected = "RGCCA")
        else if (!is.null(analysis))
            updateTabsetPanel(session, "tabset", selected = "Graphic")
    })


    observe({
        # Initial events
        for (i in c("Connection", "AVE", "Samples", "Corcircle", "Fingerprint", "Bootstrap", "'Bootstrap Summary'", "Permutation", "'Permutation Summary'", "Cross-validation"))
            hide(selector = paste0("#navbar li a[data-value=", i, "]"))
        for (i in c("run_boot", "nboot_custom", "header", "init", "navbar", "connection_save", "run_crossval_single", "kfold", "save_all", "format"))
            hide(id = i)
        for (i in c("nperm_custom", "run_perm"))
            toggle(id = i, condition = !input$supervised && !is.null(input$tau_opt) && input$tau_opt) 
        for (i in c("run_crossval", "val_custom"))
            toggle(id = i, condition = input$supervised && !is.null(input$tau_opt) && input$tau_opt) 
        toggle(id = "ncv", condition = input$supervised && input$val == "kfold" && !is.null(input$tau_opt) && input$tau_opt)
        # toggle(id = "kfold", condition = input$supervised && input$val == "kfold")
        })

    observeEvent(c(input$tau_opt, input$supervised), {
        assign("perm", NULL, .GlobalEnv)
        assign("cv", NULL, .GlobalEnv)
        toggle(id = "run_analysis", condition = !is.null(input$tau_opt) && (!input$tau_opt || (input$tau_opt && (!is.null(perm) || !is.null(cv)))))
    })

    onclick("sep", function(e) assign("clickSep", TRUE, .GlobalEnv))


    observeEvent(c(input$blocks, input$sep), {
        # blockExists for having dynamic response to input$blocks

        hide(id = "navbar")
        if (blocksExists()) {

        }

    })


    getInfile <- eventReactive(c(input$blocks, input$sep), {
        # Return the list of blocks

        # Load the blocks
        paths <- paste(input$blocks$datapath, collapse = ",")

        if (length(grep("xlsx?", paths)))
            names <- NULL
        else
            names <- paste(input$blocks$name, collapse = ",")

        cleanup_analysis_par()

        assign("blocks_unscaled",
               showWarn(
                    load_blocks(
                        file = paths,
                        names = names,
                        sep = input$sep,
                        header = TRUE
                    )
                ),
                .GlobalEnv)

        if (!is.list(blocks_unscaled))
            return(NULL)
        else {
            show(selector = "#tabset li a[data-value=RGCCA]")
            setToggle("connection")
        }

        assign(
            "blocks_without_superb",
            scaling(
                blocks_unscaled,
                ifelse(is.null(input$scale),
                        TRUE, input$scale),
                TRUE
            ),
            .GlobalEnv
        )

        # reactualiser l'analyse
        assign("nb_comp", 2, .GlobalEnv)
        assign("analysis_type", NULL, .GlobalEnv)
        assign("analysis", NULL, .GlobalEnv)
        assign("cv", NULL, .GlobalEnv)
        assign("perm", NULL, .GlobalEnv)
        assign("response", NULL, .GlobalEnv)
        assign("connection", NULL, .GlobalEnv)
        assign("response_file", NULL, .GlobalEnv)
        assign("response", load_responseShiny(), .GlobalEnv)

        assign("id_block_resp",
                length(blocks_without_superb),
                .GlobalEnv)
        blocks <- setParRGCCA(FALSE)
        assign("blocks", blocks, .GlobalEnv)
        assign("connection_file", NULL, .GlobalEnv)
        set_connectionShiny()
        setIdBlock()
        updateTabsetPanel(session, "navbar", selected = "Connection")

        return(blocks)
    })


    observeEvent(input$scale, {
        if (blocksExists()) {
            assign(
                "blocks_without_superb",
                scaling(blocks_unscaled, scale=input$scale, scale_block=TRUE),
                .GlobalEnv
            )
            setAnalysis()
            hide(id = "navbar")
        }
    })

    observeEvent(input$connection, {
        hide(id = "navbar")
        if (blocksExists()) {
            assign("connection_file",
                    input$connection$datapath,
                    .GlobalEnv)
            set_connectionShiny()
            setUiConnection()
            showWarn(message("Connection file loaded."), show = FALSE)
            assign("connection_file", NULL, .GlobalEnv)
            cleanup_analysis_par()
        }
    })
    
    cleanup_analysis_par <- function(){
        assign("analysis", NULL, .GlobalEnv)
        assign("boot", NULL, .GlobalEnv)
        assign("selected.var", NULL, .GlobalEnv)
        for (i in c("run_boot", "nboot_custom", "connection_save"))
            hide(id = i)
        for (i in c("Connection", "AVE", "Samples", "Corcircle", "Fingerprint", "Bootstrap", "'Bootstrap Summary'", "Permutation", "'Permutation Summary'", "Cross-validation"))
            hide(selector = paste0("#navbar li a[data-value=", i, "]"))
        updateTabsetPanel(session, "navbar", selected = "Connection")
        hide(id = "run_crossval_sing")
        assign("crossval", NULL, .GlobalEnv)
    }

    observeEvent(input$run_analysis, {
        if (!is.null(getInfile())) {
            assign("analysis", setRGCCA(), .GlobalEnv)
            show(selector = "#tabset li a[data-value=RGCCA]")
            for (i in c("Connection", "AVE", "Samples", "Corcircle", "Fingerprint"))
                show(selector = paste0("#navbar li a[data-value=", i, "]"))
            for (i in c("navbar", "nboot_custom", "run_boot"))
                show(id = i)
            toggle(id = "run_crossval_single", condition = !is.null(rgcca_out$call$response))
            updateTabsetPanel(session, "navbar", selected = "Connection")
            save_connection(rgcca_out$call$connection)
            # for (i in c('bootstrap_save', 'fingerprint_save', 'corcircle_save',
            # 'samples_save', 'ave_save')) setToggleSaveButton(i)
            show("connection_save")
            save(rgcca_out, file = "rgcca_result.RData")
        }
    })

    save_connection <- function(connection){
        if_superblock <- grep("superblock", rownames(connection))
        if (length(if_superblock) > 0)
            connection <- connection[-if_superblock, -if_superblock]
        write.table(connection, file = "connection.txt", sep = "\t")
    }

    observeEvent(
        c(
            input$superblock,
            input$supervised,
            input$nb_comp,
            input$scheme,
            input$init,
            input$tau_opt,
            input$analysis_type
        ),
        {
            # Observe if analysis parameters are changed

            if (blocksExists()) {

                setNamesInput("x")
                setNamesInput("response")
                assign("nb_comp", input$nb_comp, .GlobalEnv)
                hide(id = "navbar")
                setAnalysis()

                for (i in c(
                    "bootstrap_save",
                    "fingerprint_save",
                    "corcircle_save",
                    "samples_save"#,
                    # "ave_save",
                    # "connection_save"
                ))
                    hide(i)
            }

            setCompUI()

            if (!is.null(input$tau_opt) && !input$tau_opt)
                setTauUI()

        },
        priority = 10
    )

    updateSuperblock <- function(id, value)
            updateSelectizeInput(
                session,
                inputId = id,
                choices = value,
                selected = value,
                server = TRUE
            )

    observeEvent(input$supervised, {
        if (input$supervised)
            updateSuperblock("superblock", FALSE)
    })

    observeEvent(input$superblock, {
        if (input$superblock)
            updateSuperblock("supervised", FALSE)
    })

    observeEvent(input$run_boot, {
        if (blocksExists())
            getBoot()
    })
    
    observeEvent(input$run_perm, {
        if (blocksExists()) {
            getPerm()
        }
    })

    observeEvent(input$run_crossval, {
        if (blocksExists() && input$supervised)
            getCrossVal()
    })

    observeEvent(input$run_crossval_single, {
        if (blocksExists() && input$supervised)
            getCrossVal2()
    })

    observeEvent(input$names_block_x, {
        isolate({
            if (blocksExists() && !is.null(input$names_block_x)) {
                if (as.integer(input$names_block_x) > round(length(blocks))) {
                    reac_var(length(blocks))
                    assign("id_block", reac_var(), .GlobalEnv)
                } else {
                    reac_var(as.integer(input$names_block_x))
                    assign("id_block", reac_var(), .GlobalEnv)
                }
            }
        })
    }, priority = 30)

    observeEvent(c(input$superblock, input$supervised), {
        reac_var(length(blocks))
        assign("id_block", reac_var(), .GlobalEnv)
        assign("id_block_y", reac_var(), .GlobalEnv)
    }, priority = 20)

    observeEvent(input$names_block_y, {
        isolate({
            if (blocksExists() && !is.null(input$names_block_y)) {
                if (as.integer(input$names_block_y) > round(length(blocks))) {
                    reac_var(length(blocks))
                    assign("id_block_y", reac_var(), .GlobalEnv)
                } else {
                    reac_var(as.integer(input$names_block_y))
                    assign("id_block_y", reac_var(), .GlobalEnv)
                }
            }
        })
    }, priority = 30)

    observeEvent(input$names_block_response, {
        # Observe if graphical parameters are changed

        if (blocksExists()) {
            if (input$supervised || input$analysis_type == "RA")
                reac_var(as.integer(input$names_block_response))
            else
                reac_var(as.integer(input$names_block_response) - 1)

            assign("id_block_resp", reac_var(), .GlobalEnv)
            assign("nb_comp", input$nb_comp, .GlobalEnv)
            setAnalysis()
        }

    })

    observeEvent(input$save_all, {
        if (blocksExists()) {
            save_plot("samples_plot.pdf", samples())
            try(save_plot("corcircle.pdf", corcircle()), silent = TRUE)
            try(save_plot("fingerprint.pdf", fingerprint(input$indexes)), silent = TRUE)
            save_plot("AVE.pdf", ave())
            if(any(NCOL(blocks) == 1))
                compy <- 1
            else 
                compy <- 2
            save_var(rgcca_out, file = "variables.txt")
            save_ind(rgcca_out, file = "individuals.txt")
            save(analysis, file = "rgcca_result.RData")
            if(!is.null(boot))
                save_plot("bootstrap.pdf", plotBoot())
            # if(!is.null(perm))
            #     save("perm.pdf", plot_permut_2D(perm))
            msgSave()
        }
    })

    msgSave <- function()
        showWarn(message(paste("Save in", getwd())), show = FALSE)

    observeEvent(c(input$text, input$compx, input$compy, input$nb_mark), {
            if (!is.null(analysis)) {
                assign("if_text", input$text, .GlobalEnv)
                assign("compx", input$compx, .GlobalEnv)
                assign("compy", input$compy, .GlobalEnv)
                if (!is.null(input$nb_mark))
                    assign("nb_mark", input$nb_mark, .GlobalEnv)
            }
        })

    observeEvent(input$response, {
        if (!is.null(input$response)) {
            assign("response_file",
                input$response$datapath,
                .GlobalEnv)
            assign("response", load_responseShiny(), .GlobalEnv)
            setUiResponse()
            showWarn(samples(), warn = TRUE)
            showWarn(message(
                paste0(input$response$name, " loaded as a group file.")
            ),
            show = FALSE)
        }

    }, priority = 10)

    ################################################ Outputs ################################################
    
    output$connectionPlot <- renderVisNetwork({
        getDynamicVariables()
        if (!is.null(analysis)) {
            observeEvent(input$connection_save, {
                save_plot(paste0("connection.", input$format), design2)
                msgSave()
            })
            design()
        }
    })

    output$AVEPlot <- renderPlot({
        
        getDynamicVariables()
        
        if (!is.null(analysis)) {
            observeEvent(input$ave_save, {
                save_plot(paste0("ave.", input$format), ave()) 
                msgSave()
            })
            ave()
        }
        
    })

    output$samplesPlot <- renderPlotly({
        getDynamicVariables()

        if (!is.null(analysis)) {
            observeEvent(input$samples_save, {
                save_plot("samples_plot.pdf", samples())
                msgSave()
            })

            save_ind(rgcca_out, file = "individuals.txt")
            p <- samples()

            if (is(p, "gg")) {
            p <- showWarn(
                modify_hovertext(
                    plot_dynamic(p, NULL, "text", TRUE, format = input$format),
                    if_text
                ), warn = FALSE)

            if (is.null(crossval) && 
                length(unique(na.omit(response))) < 2 ||
                (length(unique(response)) > 5 &&
                !is.character2(na.omit(response))))
                p <- p %>% layout(showlegend = FALSE)

            }

            p

        }
    })

    output$corcirclePlot <- renderPlotly({
        tryCatch({
            getDynamicVariables()
    
            if (!is.null(analysis)) {
                observeEvent(input$corcircle_save, {
                    save_plot("corcircle.pdf", corcircle())
                    msgSave()
                })
                
                save_var(rgcca_out, file = "variables.txt")
                p <- corcircle()

                if (is(p, "gg")) {
                    p <- modify_hovertext(plot_dynamic(p, NULL, "text", format = input$format), if_text)
                    n <- length(p$x$data)
                    (style(
                        p,
                        hoverinfo = "none",
                        traces = c(n, n - 1)
                    ))
                }
            }
        }, error = function(e) {
        })
    })

    output$fingerprintPlot <- renderPlotly({
        tryCatch({
            getDynamicVariables()
    
            if (!is.null(analysis)) {
                observeEvent(input$fingerprint_save, {
                    save_plot("fingerprint.pdf", fingerprint(input$indexes))
                    msgSave()
                })
                modify_hovertext(plot_dynamic(fingerprint(input$indexes), type = "var1D", format = input$format), hovertext = F, type = "var1D")
            }
        }, error = function(e) {
        })
    })

    output$bootstrapPlot <- renderPlotly({
        tryCatch({
            getDynamicVariables()
            refresh <- c(input$names_block_x, id_block, input$blocks_names_custom_x)
    
            if (!is.null(analysis) & !is.null(boot)) {

                assign(
                    "selected.var", 
                    get_bootstrap(boot, compx, id_block),
                    .GlobalEnv
                )

                # observeEvent(input$bootstrap_save, {
                #     save_plot("bootstrap.pdf", plotBoot())
                #     msgSave()
                # })

               modify_hovertext(plot_dynamic(plotBoot(), type = "boot1D", format = input$format), type = "boot1D", hovertext = FALSE)
            }
        }, error = function(e) {
        })

    })

    output$bootstrapTable <- DT::renderDataTable({
        tryCatch({
            getDynamicVariables()
            refresh <- c(input$names_block_x, id_block, input$blocks_names_custom_x)
            
            if (!is.null(analysis) & !is.null(boot)) {
                
                assign(
                    "selected.var", 
                    get_bootstrap(boot, compx, id_block),
                    .GlobalEnv
                )

                df <- round(get_bootstrap(boot, compx, id_block, display_order = F), 3)[, -c(1, 3)]
                colnames(df) <- c("RGCCA weight", "Lower limit", "Upper limit", "P-value", "B.H.")

                observeEvent(input$bootstrap_t_save, {
                    write.table(df, "summary_bootstrap.txt", sep = "\t")
                    msgSave()
                })

                df
            }
        }, error = function(e) {
        })
        
    }, options = list(pageLength = 10))


    output$permutationPlot <- renderPlotly({

        getDynamicVariables()

        if (!is.null(perm)) {
        #     observeEvent(input$permutation_save, {
        #         save("perm.pdf", plot_permut_2D(perm))
        #         msgSave()
        #     })
            modify_hovertext(plot_dynamic(plot_permut_2D(perm), type = "perm", format = input$format), type = "perm", hovertext = F, perm = perm)
        }

    })
    
    output$permutationTable <- renderDataTable({
        
        getDynamicVariables()
        
        if (!is.null(perm)) {

            s_perm <- summary.perm(perm)

            observeEvent(input$permutation_t_save, {
                write.table(s_perm, "summary_permutation.txt", sep = "\t", row.names = FALSE)
                msgSave()
            })

            s_perm
        }
        
    }, options = list(pageLength = 10))

    output$cvPlot <- renderPlotly({

        getDynamicVariables()

        if (!is.null(cv)) {
            # observeEvent(input$cv_save, {
            #     save("cv.pdf", plot(cv))
            #     msgSave()
            # })
            modify_hovertext(plot_dynamic(plot(cv), type = "cv", format = input$format), type = "cv", hovertext = F, perm = cv)
        }

    })
}