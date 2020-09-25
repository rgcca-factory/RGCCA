# Author: Etienne CAMENEN
# Date: 2020
# Contact: arthur.tenenhaus@l2s.centralesupelec.fr
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.)
# and produces textual and graphical outputs (e.g. variables and samples
# plots).

rm(list = ls())
options(shiny.maxRequestSize = 30 * 1024 ^ 2)

setInfo <- function(., text) {
    shinyInput_label_embed(
        icon("question") %>%
            bs_embed_tooltip(title = text))
}

# Global variables
one_block <<- c(`Principal Component Analysis` = "PCA")
two_blocks <<- c(
        `Canonical Correlation Analysis` = 'CCA',
        `Interbattery Factor Analysis` = "IFA",
        `Partial Least Squares Regression` = 'PLS',
        `Redundancy analysis` = 'RA'
    )
multiple_blocks  <<- c(
        `Regularized Generalized CCA (RGCCA)` = 'RGCCA',
        `Sparse Generalized CCA (SGCCA)` = 'SGCCA',
        `SUM of CORrelations method` = 'SUMCOR',
        `Sum of SQuared CORrelations method` = 'SSQCOR',
        `Sum of ABSolute value CORrelations method` = 'SABSCOR',
        `SUM of COVariances method` = 'SUMCOV',
        `Sum of SQuared COVariances method` = 'SSQCOV',
        `Sum of ABSolute value COVariances method` = 'SABSCOV',
        `MAXBET` = 'MAXBET',
        `MAXBETB` = 'MAXBET-B'
    )
multiple_blocks_super  <<- c(
        `Generalized CCA (GCCA)` = 'GCCA',
        `Hierarchical PCA` = 'HPCA',
        `Multiple Factor Analysis` = 'MFA'
    )
analyse_methods  <<- list(one_block, two_blocks, multiple_blocks, multiple_blocks_super)
reac_var  <<- reactiveVal()
id_block_y <<- id_block <<- id_block_resp <<- analysis <<- connection <<- perm <<- boot <<-
boot <<- analysis_type <<- crossval <<- selected.var <<- crossval <<- NULL
clickSep <<- FALSE
if_text <<- TRUE
compx <<- 1
nb_comp <<- compy <<- 2
nb_mark <<- 100
BSPLUS <<- R.Version()$minor >= 3
ax2 <<- list(linecolor = "white",
        tickfont = list(size = 10, color = "grey"))

# config for shinyapps.io
appDir <- ifelse("packrat" %in% list.files(), "", "../../R/")
# Load functions

for (f in list.files(appDir))
        source(paste0(appDir, f))

# maxdiff-b, maxdiff, maxvar-a, maxvar-b, maxvar, niles, r-maxvar,
# rcon-pca, ridge-gca, , ssqcov-1, ssqcov-2, , sum-pca, sumcov-1, sumcov-2

load_libraries(c(
    "ggplot2",
    "scales",
    "igraph",
    "plotly",
    "visNetwork",
    "shiny",
    "shinyjs",
    "MASS",
    "rlang",
    "DT"
))

if (BSPLUS) {
    load_libraries("devtools")
    if (!("bsplus" %in% installed.packages()[, "Package"]))
        devtools::install_github("ijlyttle/bsplus", upgrade = "never")
    library("bsplus", warn.conflicts = FALSE, quiet = TRUE)
}

ui <- fluidPage(
    titlePanel("R/SGCCA - The Shiny graphical interface"),
    tags$div(
        tags$p(
            "Etienne CAMENEN, Ivan MOSZER, Arthur TENENHAUS (",
            tags$a(href = "arthur.tenenhaus@l2s.centralesupelec.fr",
            "arthur.tenenhaus@l2s.centralesupelec.fr"),
            ")"
        ),
        tags$i("Multi-block data analysis concerns the analysis of several sets of variables (blocks) observed on the same group of samples. The main aims of the RGCCA package are: to study the relationships between blocks and to identify subsets of variables of each block which are active in their relationships with the other blocks."),
        tags$br(), tags$br()
    ),
    tags$a(href = "https://github.com/rgcca-factory/RGCCA/blob/release/3.0.0/inst/shiny/tutorialShiny.md", "Go to the tutorial"),
    tags$strong("|"),
    tags$a(href = "https://www.youtube.com/watch?v=QCkEBsoP-tc", "Watch a demo", target = "_blank"),
    tags$br(), tags$br(),
    tags$style(".fa-camera {color:#c7c7c7}"),
    tags$style(".fa-camera:hover {color:#7c7c7c}"),
    tags$style("#connection_save, #ave_save {border-color:white; left: 0%}"),
    tags$style("#connection_save:hover, #ave_save:hover {background-color:white}"),
    tags$style("#connection_save:focus, #ave_save:focus {outline:none; background-color:white}"),
    tags$style("#connection_save:active, #ave_save:active {box-shadow:none}"),
    tags$style(".js-plotly-plot .plotly .modebar {left: 0%}"),
    useShinyjs(),
    sidebarLayout(sidebarPanel(
        tabsetPanel(
            id = "tabset",
            tabPanel(
                "Data",
                uiOutput("file_custom"),
                uiOutput("sep_custom"),
                checkboxInput(
                    inputId = "header",
                    label = "Consider first row as header",
                    value = TRUE
                )
            ),


            # Analysis parameters

            tabPanel(
                "RGCCA",
                uiOutput("analysis_type_custom"),
                checkboxInput(
                    inputId = "each_ncomp",
                    label = "Tune the components for each block",
                    value = FALSE
                ),
                uiOutput("nb_compcustom"),
                uiOutput("scale_custom"),
                radioButtons(
                    "init",
                    label = "Mode of initialization",
                    choices = c(SVD = "svd",
                                Random = "random"),
                    selected = "svd"
                ),

                uiOutput("superblock_custom"),
                checkboxInput(
                    inputId = "supervised",
                    label = "Supervised analysis",
                    value = FALSE
                ),

                conditionalPanel(
                    condition = "input.supervised || input.analysis_type == 'RA'",
                    uiOutput("blocks_names_response")),

                uiOutput("connection_custom"),
                uiOutput("scheme_custom"),
                uiOutput("tau_opt_custom"),
                uiOutput("each_tau_custom"),
                uiOutput("tau_custom"),
                uiOutput("val_custom"),
                sliderInput(
                    inputId = "ncv",
                    label = "Number of cross-validation",
                    min = 1,
                    max = 100,
                    value = 1,
                    step = 1
                ),
                sliderInput(
                    inputId = "kfold",
                    label = "Number of folds",
                    min = 2,
                    max = 10,
                    value = 5,
                    step = 1
                ),
                actionButton(
                    inputId = "run_crossval",
                    label = "Run cross-validation"),
                uiOutput("nperm_custom"),
                actionButton(inputId = "run_perm",
                    label = "Run permutation"),
                # sliderInput(
                #     inputId = "power",
                #     label = "Power of the factorial",
                #     min = 2,
                #     max = 6,
                #     value = 2,
                #     step = 1
                # ),
                actionButton(
                    inputId = "run_analysis",
                    label = "Run analysis"),
                uiOutput("nboot_custom"),
                actionButton(inputId = "run_boot",
                    label = "Run bootstrap"),
                actionButton(
                    inputId = "run_crossval_single",
                    label = "Evaluate the model")
            ),

            # Graphical parameters

            tabPanel(
                "Graphic",
                radioButtons(
                    "format",
                    label = "Output image format",
                    choices = c(
                        `jpeg` = "jpeg",
                        `png` = "png"#,
                        #`svg` = "svg"
                        # `tiff` = "tiff",
                        # `pdf` = "pdf"
                    ),
                    selected = "png"
                ),
                checkboxInput(
                    inputId = "text",
                    label = "Display names",
                    value = TRUE
                ),
                uiOutput("blocks_names_custom_x"),
                uiOutput("blocks_names_custom_y"),
                uiOutput("compx_custom"),
                uiOutput("compy_custom"),
                uiOutput("nb_mark_custom"),
                uiOutput("response_custom"),
                checkboxInput(
                    inputId = "show_crossval",
                    label = "Display cross-validation",
                    value = TRUE
                ),
                radioButtons(
                    "indexes",
                    label = "Type of indexes",
                    choices = c(
                        Correlation = "cor",
                        Weights = "weight")
                ),
                uiOutput("b_x_custom"),
                uiOutput("b_y_custom"),
                actionButton(inputId = "save_all", label = "Save all")
            )

        )
    ),

    mainPanel(
        tabsetPanel(
            type = "tabs",
            id = "navbar",
            tabPanel(
                "Connection",
                actionButton("connection_save", "", icon = icon("camera")),
                visNetworkOutput("connectionPlot")
            ),
            tabPanel(
                "AVE",
                actionButton("ave_save", "", icon = icon("camera")),
                plotOutput("AVEPlot")
            ),
            tabPanel(
                "Samples",
                plotlyOutput("samplesPlot", height = 500),
                actionButton("samples_save", "Save")
            ),
            tabPanel(
                "Corcircle",
                plotlyOutput("corcirclePlot", height = 500),
                actionButton("corcircle_save", "Save")
            ),
            tabPanel(
                "Fingerprint",
                plotlyOutput("fingerprintPlot", height = 700),
                actionButton("fingerprint_save", "Save")
            ),
            tabPanel(
                "Bootstrap",
                plotlyOutput("bootstrapPlot", height = 700),
                actionButton("bootstrap_save", "Save")
            ),
            tabPanel(
                "Bootstrap Summary",
                DT::dataTableOutput("bootstrapTable"),
                actionButton("bootstrap_t_save", "Save")
            ),
            tabPanel(
                "Permutation",
                plotlyOutput("permutationPlot", height = 700),
                # actionButton("permutation_save", "Save")
            ),
            tabPanel(
                "Permutation Summary",
                dataTableOutput("permutationTable"),
                actionButton("permutation_t_save", "Save")
            ),
            tabPanel(
                "Cross-validation",
                plotlyOutput("cvPlot", height = 700),
                #actionButton("cv_save", "Save")
            )
        )

    ))
)