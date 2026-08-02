# Simulate data according to the simulation scenario presented in AC-RGCCA paper and export for the package

# Simulation scenario:
# - 5 batches
# - 10 individuals per batch, exposed to different levels of a pollutant or drug (common dose effect)
# - 2 blocks, all presenting w same common dose effect + batch specific effect (possibly confounded) + random noise
#   * block 1: batch effects are based on matrix lambda_1
#   * block 2: batch effects are based on matrices lambda_1 and lambda_2 (in which effect in not uniform across individuals)

# Prepare objects and parameters
example_confounders <- list(
  blocks = NULL,
  confounders_Z_A = NULL,
  confounders_Z_B = NULL,
  common_pattern = NULL,
  labs = NULL,
  dose_groups = NULL
)

n <- 10
B <- 5
p <- c(500, 500)

alpha <- 2.5
epsi_sd <- 0.25

X <- list(block1 = matrix(0, nrow = n * B, ncol = p[1]),
          block2 = matrix(0, nrow = n * B, ncol = p[2]))
H <- list(block1 = matrix(0, nrow = p[1], 2),
          block2 = matrix(0, nrow = p[2], 2))
gamma_list <- list(block1 = list(lambda1 = matrix(0, nrow = n * B, ncol = p[1]),
                                 lambda2 = NULL,
                                 gamma = NULL),
                   block2 = list(lambda1 = matrix(0, nrow = n * B, ncol = p[2]), 
                                 lambda2 = matrix(0, nrow = n * B, ncol = p[2]),
                                 gamma = matrix(0, nrow = n * B, ncol = p[2])))

## Create factors for batch and dose
example_confounders$labs <- factor(x = rep(1:B, each = n))
example_confounders$dose_groups <- factor(x = rep(1:n, B))

# Define fixed objects
## Shared component (effect of dose on individuals)
y1 <- scale(1:n)
ones <- matrix(1, nrow = n, ncol = 1)
sigma <- 0.25 * exp(-as.matrix(dist(y1))^2/4)
set.seed(47358) # try 11, 28, 91, 214, 222, 354, 6135, 16085, 47358
y2 <- MASS::mvrnorm(n = 1, mu = rep(0, n), Sigma = sigma)

#### Save true dose effect
example_confounders$common_pattern <- data.frame("y1" = y1, "y2" = y2)

## Create dataset
### Blocks
for (j in 1:(length(p))) {
  #### Shared component
  h <- matrix(rnorm(n = 2 * p[j]), nrow = p[j], ncol = 2)
  psi <- tcrossprod(cbind(y1, y2), h)
  
  #### Batches
  for (i in 1:B) {
    #### Batch effect (+ store matrices in gamma_list)
    if (j == 1) {
      r <- matrix(rnorm(n = p[j]), nrow = p[j], ncol = 1)
      lambda_1 <- tcrossprod(ones, r)
      
      gamma <- lambda_1
      gamma_list[[j]][[1]][(1+(i-1)*n):(i*n),] <- lambda_1
      
    } else if (j == 2) {
      r <- matrix(rnorm(n = p[j]), nrow = p[j], ncol = 1)
      lambda_1 <- tcrossprod(ones, r)
      
      s <- matrix(rnorm(n = p[j]), nrow = p[j], ncol = 1)
      B_i_sample <- sample(n, size = 3)
      B_i <- rep(0, n)
      B_i[B_i_sample] <- runif(n = 3, min = 0, max = 2)
      lambda_2 <- tcrossprod(B_i, s)
      
      gamma <- lambda_1 + lambda_2
      gamma_list[[j]][[1]][(1+(i-1)*n):(i*n),] <- lambda_1
      gamma_list[[j]][[2]][(1+(i-1)*n):(i*n),] <- lambda_2
      gamma_list[[j]][[3]][(1+(i-1)*n):(i*n),] <- gamma
    } 
    
    #### Random noise
    epsilon <- matrix(rnorm(n = n * p[j], sd = epsi_sd), nrow = n, ncol = p[j])
    
    ### Store datasets
    #### Fill matrix X
    X[[j]][(1+(i-1)*n):(i*n),] <- psi + alpha * gamma + epsilon
  }
  
  #### Fill weights matrix for effect size of dose pattern
  H[[j]] <- h
}

#### Create batch model matrix
Z_a <- model.matrix(~ -1 + example_confounders$labs)
rownames(Z_a) <- rownames(X[[1]])

Z_b <- matrix(0, nrow = n * B, ncol = B * (B-1) * n / 2)
l <- 1
for (i in 1:(B-1)) {
  for (j in (i+1):B) {
    for (o in 1:n) {
      Z_b[(i-1)*n + o, l] <- 1
      Z_b[(j-1)*n + o, l] <- -1
      l <- l + 1
    }
  }
}
rownames(Z_b) <- rownames(X[[1]])

#### Center and scale matrices
X <- lapply(X, scale, center = T, scale = T)
Z_a <- scale(x = Z_a, center = T, scale = T)
Z_b <- scale(x = Z_b, center = T, scale = T)

# Store matrices
example_confounders$blocks <- X
example_confounders$confounders_Z_A <- Z_a
example_confounders$confounders_Z_B <- Z_b

# Save dataset
usethis::use_data(example_confounders, compress = "xz", overwrite = T)
