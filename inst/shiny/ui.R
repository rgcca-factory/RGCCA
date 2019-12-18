# Author: Etienne CAMENEN
# Date: 2019
# Contact: arthur.tenenhaus@l2s.centralesupelec.fr
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.)
# and produces textual and graphical outputs (e.g. variables and individuals
# plots).

rm(list = ls())
options(shiny.maxRequestSize = 30 * 1024 ^ 2)

setInfo <- function(., text) {
    shinyInput_label_embed(icon("question") %>%
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
id_block_y <<- id_block <<- id_block_resp <<- analysis <<- 
boot <<- analysis_type <<- NULL
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
    "RGCCA",
    "ggplot2",
    "scales",
    "igraph",
    "plotly",
    "visNetwork",
    "shiny",
    "shinyjs"
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
        tags$strong("Authors: "),
        tags$p(
            "Etienne CAMENEN, Ivan MOSZER, Arthur TENENHAUS (",
            tags$a(href = "arthur.tenenhaus@l2s.centralesupelec.fr",
            "arthur.tenenhaus@l2s.centralesupelec.fr"),
            ")"
        )
    ),
    tags$a(href = "https://github.com/BrainAndSpineInstitute/rgcca_Rpackage/blob/master/inst/shiny/tutorialShiny.md", "Go to the tutorial"),
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

                conditionalPanel(condition = "input.supervised || input.analysis_type == 'RA'",
                                uiOutput("blocks_names_response")),

                uiOutput("connection_custom"),
                uiOutput("tau_opt_custom"),
                uiOutput("tau_custom"),
                uiOutput("scheme_custom"),

                sliderInput(
                    inputId = "boot",
                    label = "Number of boostraps",
                    min = 5,
                    max = 100,
                    value = 10,
                    step = 5
                ),
                actionButton(inputId = "run_boot",
                            label = "Run bootstrap"),
                actionButton(inputId = "run_analysis",
                            label = "Run Analysis")
            ),

            # Graphical parameters

            tabPanel(
                "Graphic",
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

                actionButton(inputId = "save_all",
                            label = "Save all")
            )

        )
    ),

    mainPanel(
        tabsetPanel(
            type = "tabs",
            id = "navbar",
            tabPanel(
                "Connection",
                visNetworkOutput("connectionPlot"),
                actionButton("connection_save", "Save")
            ),
            tabPanel(
                "AVE",
                plotOutput("AVEPlot"),
                actionButton("ave_save", "Save")
            ),
            tabPanel(
                "Samples",
                plotlyOutput("samplesPlot", height = 700),
                actionButton("samples_save", "Save")
            ),
            tabPanel(
                "Corcircle",
                plotlyOutput("corcirclePlot"),
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
            )
        )

    ))
)
