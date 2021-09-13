\name{RGCCA-package}
\alias{RGCCA-package}
\alias{RGCCA}
\docType{package}
\title{
\packageTitle{RGCCA}
}
\description{
\packageDescription{RGCCA}
}

\author{
\packageAuthor{RGCCA}

Maintainer: \packageMaintainer{RGCCA}
}

\keyword{ package }
\seealso{\link{rgcca},\link{bootstrap},\link{rgcca_cv},\link{rgcca_permutation}
}
\examples{
# Create the dataset
data(Russett)
blocks = list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11]
)

# Tune the model to find the optimal regularization parameters
perm.out = rgcca_permutation(blocks, n_cores = 1, par_type = "tau", n_perms = 10)

print(perm.out)
plot(perm.out)

# fit the rgcca model based on these optimal parameters
fit.rgcca = rgcca(perm.out)

# Build bootstrap confidence intervals.
b = bootstrap(fit.rgcca, n_cores = 1, n_boot = 10)
plot(b, n_cores = 1)
}
