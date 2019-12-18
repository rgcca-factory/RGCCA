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
        # Set dynamically the maximum number of component that should be used
        # in the analysis

        sliderInput(
            inputId = "nb_comp",
            label = "Number of components",
            min = 2,
            max = getDefaultComp(),
            value = 2,
            step = 1
        )

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
        isolate(uiComp("x", 1, input$navbar != "Fingerprint"))
    })

    output$compy_custom <- renderUI({
        refresh <- refreshAnalysis()
        uiComp("y", 2)
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

    ################################################ UI function ################################################

    setTauUI <- function(superblock = NULL) {
        refresh <- c(input$superblock, input$supervised)

        if (!is.null(input$analysis_type) &&
            input$analysis_type == "SGCCA") {
            par_name <- "Sparsity"
            cond <- "input.analysis_type == SGCCA"
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
                                    1, input[[paste0("tau", i)]]),
                                step = .01
                            )
                        }))
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

    uiComp <- function(x, y, bool = TRUE) {
        label <- "Component"

        if (bool)
            label <- paste0("Component for the ", x, "-axis")

        sliderInput(
            inputId = paste0("comp", x),
            label = label,
            min = 1,
            max = input$nb_comp,
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

    output$tau_opt_custom <- renderUI({
        ui <- checkboxInput(inputId = "tau_opt",
                            label = "Use an optimal tau",
                            value = TRUE)

        if (BSPLUS)
            ui <- shinyInput_label_embed(
                ui,
                icon("question") %>%
                    bs_embed_tooltip(title = "A tau near 0 maximize the the
                    correlation whereas a tau near 1 maximize the covariance")
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
                    '1' for a connection between two blocks, or '0' otherwise.")
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

    getMinComp <- function() {
        # Get the maximum number of component allowed in an analysis based
        # on the minimum umber of column among the blocks

        if (!is.null(input$blocks)) {
            blocks <- getInfile()
            if (!is.null(blocks)) {
                min <- min(sapply(blocks, NCOL))
                if (min > 5)
                    return(5)
                else
                    return(min)
            }
        }
        return(2)

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

    getDefaultComp <- function() {
        # Set the maximum of component to the minimum
        # number of column among the blocks but not higher than 5

        min <- getMinComp()

        if (min < 5)
            return (min)
        else
            return (5)
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
            input$boot,
            input$text,
            input$names_block_response,
            input$supervised,
            input$run_analysis,
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
        isolate({
            plot_ind(
                rgcca = rgcca_out,
                resp = response,
                compx = compx,
                compy = compy,
                i_block = id_block,
                text = if_text,
                i_block_y = id_block_y,
                reponse_name = getExtension(input$response$name)
            )
        })
    }

    corcircle <- function()
        plot_var_2D(
            rgcca = rgcca_out,
            compx = compx,
            compy = compy,
            i_block = id_block,
            text = if_text,
            n_mark = nb_mark
        )

    fingerprint <- function()
        plot_var_1D(
            rgcca = rgcca_out,
            comp = compx,
            n_mark = nb_mark,
            i_block = id_block,
            type = "cor"
        )

    ave <- function()
        plot_ave(rgcca = rgcca_out)

    design <- function()
        plot_network2(rgcca_out)

    plotBoot <- function()
        plot_bootstrap(boot, compx, nb_mark, id_block)

    ################################################ Analysis ################################################

    getTau <- function() {
        tau <- integer(0)
        for (i in 1:(length(blocks_without_superb) + ifelse(input$superblock, 1, 0)))
            tau <- c(tau, input[[paste0("tau", i)]])

        return(tau)
    }

    setParRGCCA <- function(verbose = TRUE) {
        blocks <- blocks_without_superb

        if (is.null(analysis_type) | is.null(input$analysis_type))
            analysis_type <- "RGCCA"
        else
            analysis_type <- input$analysis_type

        # Tau is set to optimal by default
        if (is.null(input$tau_opt)  || (input$tau_opt && analysis_type != "SGCCA"))
            tau <- "optimal"
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
                ( analysis_type %in% c("PCA", "RGCCA", "SGCCA")) ||
                analysis_type %in% multiple_blocks_super )
            blocks <- c(blocks, Superblock = list(Reduce(cbind, blocks)))

        assign("tau", tau, .GlobalEnv)
        assign("analysis_type", analysis_type, .GlobalEnv)

        return(blocks)
    }

    setRGCCA <- function() {
        # Load the analysis

        isolate({
            if (length(grep("[SR]GCCA", analysis_type)) == 1 && !input$tau_opt)
                tau <- getTau()
        })

        if (!is.null(input$supervised) && input$supervised)
            response <- input$supervised
        else
            response <- NULL

        assign("rgcca_out",
               showWarn(
                   rgcca.analyze(
                        blocks_without_superb,
                        connection = connection,
                        response = input$names_block_response,
                        superblock = (!is.null(input$supervised) &&
                            !is.null(input$superblock) && input$superblock),
                        tau = tau,
                        ncomp = input$nb_comp,
                        scheme = input$scheme,
                        scale = FALSE,
                        init = input$init,
                        bias = TRUE,
                        type = analysis_type
                    )
                ),
                .GlobalEnv)

        # print(length(rgcca_out))
        # 
        # if (length(rgcca_out) == 1) {
        #     assign("analysis", NULL, .GlobalEnv)
        #     return(NULL)
        # }

        #getBoot()
    }

    getBoot <- function()
        assign(
            "boot",
            bootstrap(
                blocks,
                input$boot,
                connection,
                tau,
                ncomp,
                input$scheme,
                input$scale,
                input$init,
                TRUE,
                analysis_type
            ),
            .GlobalEnv
        )


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
            assign("analysis", NULL, .GlobalEnv)
            assign("boot", NULL, .GlobalEnv)
        }

    }

    setAnalysis <- function() {
        blocks <- setParRGCCA()

        if (!is.null(blocks)) {
            assign("analysis", NULL, .GlobalEnv)
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
                    && length(input$blocks$datapath) > 2
            ),
            id = id)


    setToggle2 <- function(id)
            toggle(
                condition = (input$analysis_type %in% c("RA", "RGCCA","SGCCA")),
                   id = id)


    setToggleSaveButton <- function(id)
            toggle(condition = !is.null(analysis), id = id)


    observe({
        # Event related to input$analysis_type
        toggle(
            condition = (input$analysis_type == "RGCCA"),
            id = "tau_opt")
        setToggle("tau_custom")
        setToggle("scheme")
        setToggle("superblock")
        setToggle("connection")
        setToggle2("blocks_names_response")
        setToggle("supervised")
        hide(selector = "#tabset li a[data-value=Graphic]")
        toggle(
            condition = (length(input$blocks$datapath) > 1), 
            id = "blocks_names_custom_x")
        toggle(
            condition = (length(input$blocks$datapath) > 1), 
            id = "blocks_names_custom_y")
    })


    observeEvent(c(input$navbar, input$tabset), {
        toggle(
            condition = (input$navbar == "Fingerprint"),
               id = "nb_mark_custom")
        toggle(
            condition = (input$navbar != "Fingerprint"),
               id = "text")
        toggle(
            condition = (input$navbar != "Fingerprint"),
               id = "compy_custom")
        toggle(
            condition = (input$navbar == "Samples" && length(input$blocks$datapath) > 1),
               id = "blocks_names_custom_y")
        toggle(
            condition = (input$navbar == "Samples"),
               id = "response")
        toggle(
            condition = (
                !is.null(analysis) && !input$navbar %in% c("Connection", "AVE")
            ),
            selector = "#tabset li a[data-value=Graphic]"
        )
    })


    observeEvent(input$navbar, {
        if (!is.null(analysis) && input$navbar %in% c("Connection", "AVE"))
            updateTabsetPanel(session, "tabset", selected = "RGCCA")
        else if (!is.null(analysis))
            updateTabsetPanel(session, "tabset", selected = "Graphic")
    })


    observe({
        # Initial events

        hide(selector = "#tabset li a[data-value=RGCCA]")
        hide(selector = "#navbar li a[data-value=Bootstrap]")
        hide(id = "run_boot")
        hide(id = "boot")
        hide(id = "header")
        hide(id = "init")
        hide(id = "navbar")
        hide(id = "connection_save")
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
        names <- paste(input$blocks$name, collapse = ",")

        assign("analysis", NULL, .GlobalEnv)
        hide(id = "navbar")

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
                scaling(blocks_unscaled, input$scale, TRUE),
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
            assign("analysis", NULL, .GlobalEnv)
        }
    })


    observeEvent(input$run_analysis, {
        if (!is.null(getInfile())) {
            assign("analysis", setRGCCA(), .GlobalEnv)

            show(id = "navbar")

            # for (i in c('bootstrap_save', 'fingerprint_save', 'corcircle_save',
            # 'samples_save', 'ave_save')) setToggleSaveButton(i)
        }
    })

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
                    "samples_save",
                    "ave_save",
                    "connection_save"
                ))
                    hide(i)
            }

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
            save_plot("corcircle.pdf", corcircle())
            save_plot("fingerprint.pdf", fingerprint())
            save_plot("AVE.pdf", ave())
            save_var(rgcca_out, 1, 2)
            save_ind(rgcca_out, 1, 2)
            save(analysis, file = "rgcca_result.RData")
            msgSave()
        }
    })

    msgSave <- function()
            showWarn(message(paste("Save in", getwd())), show = FALSE)

    observeEvent(c(input$text, input$compx, input$compy, input$nb_mark),
                {
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

    output$samplesPlot <- renderPlotly({
        getDynamicVariables()

        if (!is.null(analysis)) {
            observeEvent(input$samples_save, {
                save_plot("samples_plot.pdf", samples())
                msgSave()
            })

            p <- showWarn(
                modify_hovertext(
                    plot_dynamic(samples(), NULL, "text", TRUE, TRUE),
                    if_text
                ), warn = FALSE)

            if (length(unique(na.omit(response))) < 2 || (length(unique(response)) > 5 && !is.character2(na.omit(response))))
                p <- p %>% layout(showlegend = FALSE)
            p

        }
    })

    output$corcirclePlot <- renderPlotly({

        getDynamicVariables()

        if (!is.null(analysis)) {
            observeEvent(input$corcircle_save, {
                save_plot("corcircle.pdf", corcircle())
                msgSave()
            })

            p <- modify_hovertext(plot_dynamic(corcircle(), NULL, "text"), if_text)
            n <- length(p$x$data)
            (style(
                p,
                hoverinfo = "none",
                traces = c(n, n - 1)
            ))
        }

    })

    output$fingerprintPlot <- renderPlotly({

        getDynamicVariables()

        if (!is.null(analysis)) {

            observeEvent(input$fingerprint_save, {
                save_plot("fingerprint.pdf", fingerprint())
                msgSave()
            })

            p <- modify_mousehover(plot_dynamic(fingerprint(), ax2, "text"))
            n <- sapply(p$x$data, function(x) !is.null(x$orientation))

            for (i in 1:length(n[n]))
                p$x$data[[i]]$text <-
                round(as.double(
                    sub(
                        "order: .*<br />df\\[, 1\\]: (.*)<.*",
                        "\\1\\",
                        p$x$data[[i]]$text
                    )
                ), 3)
            p
        }

    })


    output$AVEPlot <- renderPlot({

        getDynamicVariables()

        if (!is.null(analysis)) {
            observeEvent(input$ave_save, {
                save_plot("AVE.pdf", ave())
                msgSave()
            })
            ave()
        }

    })

    output$connectionPlot <- renderVisNetwork({
        getDynamicVariables()
        if (!is.null(analysis)) {
            design()
        }
    })

    output$bootstrapPlot <- renderPlotly({

        getDynamicVariables()

        if (!is.null(analysis) & !is.null(boot)) {
            observeEvent(input$bootstrap_save, {
                save_plot("bootstrap.pdf", plotBoot())
                msgSave()
            })
            plot_dynamic_histogram(plotBoot())
        }

    })

}
