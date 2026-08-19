<!-- badges: start -->
  [![CRAN status](https://www.r-pkg.org/badges/version/RGCCA)](https://CRAN.R-project.org/package=RGCCA)
  [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

<!-- badges: end -->

# AC-RGCCA

##### Version: 3.0.2

##### Authors of RGCCA package:
Fabien GIRKA, Etienne CAMENEN,  Caroline PELTIER, Vincent GUILLEMOT, Arnaud GLOAGUEN, Laurent LE BRUSQUET, Arthur TENENHAUS

##### AC-RGCCA Authors:
Elen GOUJON, Arthur TENENHAUS, Laurent LE BRUSQUET, Sandrine FRELON, Olivier ARMANT, Imène GARALI

##### Key-words:
Regularized Generalized Canonical Correlation Analysis, multi-block data analysis, confounding variation

##### Contact:
arthur.tenenhaus@centralesupelec.fr

##### Short description
The RGCCA package performs multiblock component methods (PCA, CCA, PLS, MCOA, GCCA, CPCA, MAXVAR, R/SGCCA, etc.) and produces graphical outputs (e.g. variables and individuals plots) and statistics to assess the robustness/significance of the analysis. Additionally, AC-RGCCA extends RGCCA's framework to allow accounting for confounding variation within multiblock data analysis.

---

## Contents
  - [Description](#description)
  - [Algorithm for RGCCA](#algorithm-for-rgcca)
  - [AC-RGCCA](#ac-rgcca)
  - [Installation](#installation)
  - [Installation of the development branch AC-RGCCA](#installation-of-the-development-branch-ac-rgcca)
  - [References](#references)

## Description
A package for multiblock data analysis (RGCCA - Regularized Generalized Canonical Correlation Analysis) as described in [1-4]. The software produces graphical outputs and statistics to assess the robustness/significance of the analysis. The AC-RGCCA method, presented in [8], allows handling confounding variation directly within multiblock data analysis. This method is implemented in the RGCCA package and accessible on the AC-RGCCA development branch. 

## Algorithm for RGCCA
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

## AC-RGCCA
AC-RGCCA extends the RGCCA framework to adjust for unwanted variation. This new method aims to extract components that summarize the information present within the blocks and shared across connected blocks, while at the same time limiting the amount of confounding variation captured by components. Building on the existing framework for RGCCA presented in the previous section, we now introduce the mathematical objects used by AC-RGCCA and the associated optimization problem.

For each block $j=1, \dots, J$, we assume that a set of $q_j$ covariates pollutes the study of the information contained in block $\mathbf X_j$. These confounding variables observed on the $n$ individuals are represented in a column-centered $n \times q_j$ matrix $\mathbf Z_j = [\mathbf z_j^{(1)}, \ldots, \mathbf z_j^{(q_j)}]$. In presence of confounders, we wish to limit the amount of unwanted variation captured by the components. Based on this philosophy, AC-RGCCA thus implements the following optimization problem:
$$\underset{\mathbf a_1, \dots, \mathbf a_J}{\text{maximize}} \sum_{j, k = 1}^J c_{jk} g(\text{cov}(\mathbf X_j \mathbf a_j, \mathbf X_k \mathbf a_k)) - \sum_{l=1}^J \frac{\gamma_l}{n} \mathbf a_l^\top \mathbf X_l^\top \mathbf K_l \mathbf X_l \mathbf a_l$$
$$\text{ s.t. } (1 - \tau_j)\text{var}(\mathbf X_j \mathbf a_j) + \tau_j \Vert \mathbf a_j \Vert^2 = 1, ~ j = 1, \dots, J.$$
For each block $\mathbf X_j$, we associate a penalty parameter $\gamma_j \geq 0$ used to control the regularization imposed, and a $n \times n$ confounders kernel matrix $\mathbf K_j$, with elements $[K_j]_{ii'}$ measuring a similarity between observations $i$ and $i'$ in $\mathbf Z_j$. Thanks to RGCCA-specific parameters $\mathbf C$, $g$, $\{\tau_j\}_{j = 1, \dots, J}$, and the chosen deflation strategy, this new formulation allows controlling for confounding variation in various multiblock scenarios, including correlation-based or supervised models with a response block. 

A natural choice for the confounders kernel matrix is the linear kernel $\mathbf K_j = \mathbf Z_j \mathbf Z_j^\top$. With this kernel, the quantity $\mathbf a_j^\top \mathbf X_j^\top \mathbf Z_j \mathbf Z_j^\top \mathbf X_j \mathbf a_j$ is proportional to $\sum_{k=1}^{q_j} {\text{cov}}^2 (\mathbf z_j^{(k)}, \mathbf X_j \mathbf a_j)$, and the corresponding penalty therefore encourages the component $\mathbf X_j \mathbf a_j$ to be orthogonal to the confounders in $\mathbf Z_j$. Importantly, we note that more sophisticated kernel functions can also be used.

Two arguments are introduced within the `rgcca` function for users to customize the AC-RGCCA methods to their analysis:
- `gamma_confounders`: The **penalty tuning parameters** $\gamma_1, \dots, \gamma_J \geq 0$ control the strength of the penalty that will be applied on each block $j=1, \dots, J$. Choosing $\gamma_1 = \dots = \gamma_J = 0$ imposes no penalization and leads to performing classical RGCCA. On the other hand, a large enough value will shift the solutions towards being orthogonal to the confounding variables (if $\mathbf K_j = \mathbf Z_j \mathbf Z_j^\top$). Tuning $\boldsymbol \gamma = (\gamma_1, \dots, \gamma_J)$ can be achieved by studying the ratio $R(\boldsymbol \gamma)$:
  $$R(\boldsymbol \gamma) = \left( \sum_{j=1}^J \mathbf a_j^\top \mathbf X_j^\top \mathbf K_j \mathbf X_j \mathbf a_j \right)/\left( \sum_{j,k=1}^J c_{jk}^\star \mathbf a_j^\top \mathbf X_j^\top \mathbf X_k \mathbf a_k \right).$$
  In this ratio, $\mathbf C^\star$ is the connection matrix used for the model with all its diagonal elements $c_{jj}^\star, j=1, \dots, J$ changed to $1$; this allows to include the variance explained by each block in the denominator, in addition to the covariance terms. This function measures the relative importance of the unwanted variation captured by the components in relation to the total variation extracted (with terms $j = k$ of the denominator sum giving the component variance, and terms $j \neq k$ giving the covariance between components of connected blocks). We suggest selecting the smallest value $\gamma$ such that $R(\gamma) < 0.05 R(\gamma = 0)$ in the components studied, with the possibility of users tailoring the threshold level of $0.05$ to their application.
- `confounders`: The **confounders matrices** incode the unwanted information, whose influence on the estimated components is to be minimized by AC-RGCCA. It can be given as either $\mathbf Z_1, \dots, \mathbf Z_J$, or directly $\mathbf K_1, \dots, \mathbf K_J$. If users input $\mathbf Z_1, \dots, \mathbf Z_J$, a linear kernel will be applied to compute the confounders kernel matrices. Otherwise, if users wish to adapt AC-RGCCA's optimization problem to their biological question, they are free to create and input their own confounders kernel matrices $\mathbf K_1, \dots, \mathbf K_J$.

See the documentation of `rgcca` in the AC-RGCCA branch for examples on how to use the method, and further information about the function's arguments.

AC-RGCCA is implemented on the namesake development branch within the RGCCA GitHub repository. Please see below for its installation. 

## Installation
Required:

- Software: R (≥ 3.2.0)

- R libraries: see the [DESCRIPTION](https://github.com/rgcca-factory/RGCCA/blob/main/DESCRIPTION) file.

```
install.packages("RGCCA")
```

See the [vignette](https://rgcca-factory.github.io/RGCCA/articles/RGCCA.pdf) for an introduction to the package.


## Installation of the development branch AC-RGCCA
Required:

- Software: R (≥ 3.2.0)

- R libraries: see the [DESCRIPTION](https://github.com/rgcca-factory/RGCCA/blob/AC-RGCCA/DESCRIPTION) file.

- The R library `devtools`.

```
remove.packages("RGCCA")
devtools::install_github(repo = "rgcca-factory/RGCCA", ref = "AC-RGCCA")
```

## References
1. Tenenhaus, M., Tenenhaus, A., & Groenen, P. J. (2017). Regularized generalized canonical correlation analysis: a framework for sequential multiblock component methods. Psychometrika, 82(3), 737-777.
2. Tenenhaus, A., Philippe, C., & Frouin, V. (2015). Kernel generalized canonical correlation analysis. Computational Statistics & Data Analysis, 90, 114-131.
3. Tenenhaus, A., Philippe, C., Guillemot, V., Le Cao, K. A., Grill, J., & Frouin, V. (2014). Variable selection for generalized canonical correlation analysis. Biostatistics, 15(3), 569-583.
4. Tenenhaus, A., & Tenenhaus, M. (2011). Regularized generalized canonical correlation analysis. Psychometrika, 76(2), 257.
5. Van de Geer, J. P. (1984). Linear relations among K sets of variables. Psychometrika, 49(1), 79-94.
6. Schäfer, J., & Strimmer, K. (2005). A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statistical applications in genetics and molecular biology, 4(1).
7. Tenenhaus, A., & Tenenhaus, M. (2014). Regularized generalized canonical correlation analysis for multiblock or multigroup data analysis. European Journal of operational research, 238(2), 391-403.
8. Goujon, E., Tenenhaus, A., Le Brusquet, L., Frelon, S., Armant, O., & Garali, I. (To be published) AC-RGCCA: Adjusting for Confounding Variation Within Multiblock Data Analysis.
