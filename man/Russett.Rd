% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Russett.R
\docType{data}
\name{Russett}
\alias{Russett}
\title{Russett data}
\format{
A data frame with 47 rows and 12 variables.
}
\usage{
data(Russett)
}
\description{
The Russett data set (Russett, 1964) is studied in Gifi (1990). Three
blocks of variables have been defined for 47 countries. The first block
is related to "Agricultural Inequality", the second to
"Industrial Development", and the last one describes the
"Political Instability". Russett collected this data to study
relationships between Agricultural Inequality, Industrial Development and
Political Instability. Russett's hypotheses can be formulated as follows:
It is difficult for a country to escape dictatorship when its agricultural
inequality is above-average and its industrial development below-average.

\describe{

\item{X1}{Agricultural Inequality \itemize{
     \item GINI: Inequality of land distribution,
     \item FARM: Percentage of farmers that own half of the land,
     \item RENT: Percentage of farmers that rent all their land.
}}

\item{X2}{Industrial Development \itemize{
     \item GNPR: Gross national product per capita ($1955),
     \item LABO: Percentage of labor forced employed in agriculture.
}}

\item{X3}{Political Instability \itemize{
     \item INST: Instability of executive (45-61),
     \item ECKS: Number of violent internal war incidents (46-61),
     \item DEAT: Number of people killed as a result of civic group
violence (50-62),
     \item DEMOSTAB: Stable democracy,
     \item DEMOINST: Unstable democracy,
     \item DICTATOR: Dictatorship.
}}
}
}
\examples{
#Loading of the Russett dataset
data(Russett)
#Russett is partitioned into three blocks (X_agric, X_ind, X_polit)
X_agric <- Russett[, c("gini", "farm", "rent")]
X_ind <- Russett[, c("gnpr", "labo")]
X_polit <- Russett[, c("inst", "ecks", "death", "demostab",
                       "demoinst", "dictator")]
A <- list(X_agric, X_ind, X_polit)
}
\references{
Russett B.M. (1964), Inequality and Instability: The Relation of
Land Tenure to Politics, World Politics 16:3, 442-454.

Gifi, A. (1990), Nonlinear multivariate analysis,
Chichester: Wiley.
}
\keyword{datasets}
