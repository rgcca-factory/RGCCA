<!-- badges: start -->
  [![CRAN status](https://www.r-pkg.org/badges/version/RGCCA)](https://CRAN.R-project.org/package=RGCCA)
  [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

<!-- badges: end -->

# R/SGCCA

##### Version: 3.0.2

##### Authors:
Fabien GIRKA, Etienne CAMENEN,  Caroline PELTIER, Vincent GUILLEMOT, 
Arnaud GLOAGUEN, Laurent LE BRUSQUET, Arthur TENENHAUS

##### Key-words:
Regularized Generalized Canonical Correlation Analysis, multi-block data analysis

##### Contact:
arthur.tenenhaus@centralesupelec.fr

##### Short description
Performs multiblock component methods (PCA, CCA, PLS, MCOA, GCCA, CPCA, MAXVAR, R/SGCCA, etc.) and produces graphical outputs (e.g. variables and individuals plots) and statistics to assess the robustness/significance of the analysis.

---

## Contents
  - [Description](#description)
  - [Algorithm](#algorithm)
  - [Installation](#installation)
  - [Installation of a development branch from the git repository](#installation-of-a-development-branch-from-the-git-repository)
  - [References](#references)

## Descriptiont
A package for multiblock data analysis (RGCCA - Regularized Generalized Canonical Correlation Analysis) as described in [1-4]. The software produces graphical outputs and statistics to assess the robustness/significance of the analysis.

## Algorithm
We consider $J$ data matrices $\mathbf X_1 , \dots, \mathbf X_J$. Each $n \times p_j$ data matrix 
$\mathbf X_j = \left[ x_{j1}, \dots, x_{jp_j} \right]$ 
is called a block and represents a set of $p_j$ variables observed on $n$ individuals. The number and the nature of the variables may differ from one block to another, but the individuals must be the same across blocks. We assume that all variables are centered. The objective of RGCCA is to find, for each block, a weighted composite of variables (called block component) $\mathbf y_j = \mathbf X_j  \mathbf a_j, ~ j = 1 ,..., J$ (where $\mathbf a_j$ is a column-vector with $p_j$ elements) summarizing the relevant information between and within the blocks. The block components are obtained such that (i) block components explain well their own block and/or (ii) block components that are assumed to be connected are highly correlated. In addition, RGCCA integrates a variable selection procedure, called SGCCA, allowing the identification of the most relevant features.

RGCCA subsumes fifty years of multiblock component methods and is defined as the following optimization problem:
$$\underset{\mathbf a_1, \dots, \mathbf a_J}{\text{maximize}} \sum_{j, k = 1}^J c_{jk} g(\text{cov}(\mathbf X_j \mathbf a_j, \mathbf X_k \mathbf a_k)) \text{ s.t. } (1 - \tau_j)\text{var}(\mathbf X_j \mathbf a_j) + \tau_j \Vert \mathbf a_j \Vert^2 = 1, ~ j = 1, \dots, J.$$

- The **scheme function** $g$ is any continuous convex function and allows to consider different optimization criteria. Typical choices of $g$ are the identity (horst scheme, leading to maximizing the sum of covariances between block components), the absolute value (centroid scheme, yielding maximization of the sum of the absolute values of the covariances), the square function (factorial scheme, thereby maximizing the sum of squared covariances), or, more generally, for any even integer $m$, $g(x) = x^m$ ($m$-scheme, maximizing the power of $m$ of the sum of covariances). The horst scheme penalizes structural negative correlation between block components while both the centroid scheme and the $m$-scheme enable two components to be negatively correlated. According to [5], a fair model is a model where all blocks contribute equally to the solution in opposition to a model dominated by only a few of the $J$ sets. If fairness is a major objective, the user must choose $m = 1$. $m > 1$ is preferable if the user wants to discriminate between blocks. In practice, $m$ is equal to 1, 2 or 4. The higher the value of $m$ the more the method acts as block selector [5].

- The **design matrix** $\mathbf C$ is a symmetric $J \times J$ matrix of nonnegative elements describing the network of connections between blocks the user wants to take into account. Usually, $c_{jk} = 1$ for two connected blocks and 0 otherwise.

- The $\tau_j$ are called **shrinkage parameters** or **regularization parameters** ranging from 0 to 1. $\tau_j$ enables interpolate smoothly between maximizing the covariance and maximizing the correlation. Setting the $\tau_j$ to 0 will force the block components to unit variance ($\text{var}(\mathbf X_j \mathbf a_j) = 1$). In this case, the covariance criterion boils down to the correlation. The correlation criterion is better in explaining the correlated structure across datasets, thus discarding the variance within each individual dataset. Setting $\tau_j$ to 1 will normalize the block weight vectors ($\Vert \mathbf a_j \Vert = 1$), which applies the covariance criterion. A value between 0 and 1 will lead to a compromise between the two first options and correspond to the following constraint $(1 − \tau_j)  \text{var}(\mathbf X_j \mathbf a_j) + \tau_j \Vert \mathbf a_j \Vert^2 = 1$. In the RGCCA package, for each block, the determination of the shrinkage parameter can be made fully automatic by using the analytical formula proposed by (Schäfer and Strimmer 2005 [6]), by permutation or K fold cross-validation.
Moreover, we can define the choice of the shrinkage parameters by providing interpretations on the properties of the resulting block components:

    - $\tau_j = 1$ yields the maximization of a covariance-based criterion. It is recommended when the user wants a stable component (large variance) while simultaneously taking into account the correlations between blocks. The user must, however, be aware that variance dominates over correlation.

    - $\tau_j = 0$ yields the maximization of a correlation-based criterion. It is recommended when the user wants to maximize correlations between connected components. This option can yield unstable solutions in case of multi-collinearity and cannot be used when a data block is rank deficient (e.g. $n < p_j$).

    - $0 < \tau_j < 1$ is a good compromise between variance and correlation: the block components are simultaneously stable and as well correlated as possible with their connected block components. This setting can be used when the data block is rank deficient.

The quality and interpretability of the RGCCA block components $\mathbf y_j = \mathbf X_j \mathbf a_j, ~ j = 1 , \dots, J$ are likely affected by the usefulness and relevance of the variables of each block. Accordingly, it is an important issue to identify within each block a subset of significant variables which are active in the relationships between blocks. **SGCCA** extends RGCCA to address this issue of variable selection. Specifically, RGCCA with all $\tau_j$ equal to 1 is combined with an L1-penalty that gives rise to SGCCA [3]. The SGCCA optimization problem is defined with $s_j$, a user defined positive constant that determines the amount of sparsity through the additional constraint $\Vert \mathbf a_j \Vert_1 \leq s_j, ~ j = 1, \dots, J$. The smaller the $s_j$, the larger the degree of sparsity for $\mathbf a_j$. The sparsity parameter $s_j$ is usually set by cross-validation or permutation. Alternatively, values of $s_j$ can simply be chosen to result in desired amounts of sparsity.

## Installation
Required:

- Software: R (≥ 3.2.0)

- R libraries: see the [DESCRIPTION](https://github.com/rgcca-factory/RGCCA/blob/main/DESCRIPTION) file.

```
install.packages("RGCCA")
```

See the [vignette](https://rgcca-factory.github.io/RGCCA/articles/RGCCA.pdf) for an introduction to the package.


## Installation of a development branch from the git repository
Required:

- Software: R (≥ 3.2.0)

- R libraries: see the [DESCRIPTION](https://github.com/rgcca-factory/RGCCA/blob/main/DESCRIPTION) file.

- The R library `devtools`.

```
remove.packages("RGCCA")
devtools::install_github(repo="https://github.com/rgcca-factory/RGCCA.git", ref = "main")
```

## References
1. Tenenhaus, M., Tenenhaus, A., & Groenen, P. J. (2017). Regularized generalized canonical correlation analysis: a framework for sequential multiblock component methods. Psychometrika, 82(3), 737-777.
2. Tenenhaus, A., Philippe, C., & Frouin, V. (2015). Kernel generalized canonical correlation analysis. Computational Statistics & Data Analysis, 90, 114-131.
3. Tenenhaus, A., Philippe, C., Guillemot, V., Le Cao, K. A., Grill, J., & Frouin, V. (2014). Variable selection for generalized canonical correlation analysis. Biostatistics, 15(3), 569-583.
4. Tenenhaus, A., & Tenenhaus, M. (2011). Regularized generalized canonical correlation analysis. Psychometrika, 76(2), 257.
5. Van de Geer, J. P. (1984). Linear relations among K sets of variables. Psychometrika, 49(1), 79-94.
6. Schäfer, J., & Strimmer, K. (2005). A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statistical applications in genetics and molecular biology, 4(1).
7. Tenenhaus, A., & Tenenhaus, M. (2014). Regularized generalized canonical correlation analysis for multiblock or multigroup data analysis. European Journal of operational research, 238(2), 391-403.
