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

# Tune the model to find the best sparsity coefficients (all the blocs are connected together)
perm = rgcca_permutation(blocks, n_cores = 1, par_type = "sparsity", n_perms = 10)
print(perm)
plot(perm)

res_sgcca = rgcca(blocks, type = "sgcca", sparsity = perm$bestpenalties)
plot(res_sgcca, type = "network")
plot(res_sgcca, type = "ave")

# Select the most significant variables
b = bootstrap(res_sgcca, n_cores = 1, n_boot = 100)
plot(b, n_cores = 1)

# Tune the model for explaining the politic block (politic connected to the two other blocks)
cv = rgcca_cv(blocks, response = 3, ncomp = 2, n_cores = 1)
print(cv)
plot(cv)

res_rgcca = rgcca(blocks, response = 3, ncomp = 2, tau = cv$bestpenalties)
plot(res_rgcca, type = "both")

b = bootstrap(res_rgcca, n_cores = 1, n_boot = 10)
plot(b, n_cores = 1)
}