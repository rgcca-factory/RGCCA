#' Regularized Generalized Canonical Correlation Analysis (RGCCA)
#'
#' RGCCA is a general statistical framework for multiblock data analysis.
#' The rgcca() function implements this framework and is the main entry point
#' of the package.
#'
#' @details
#' Given \eqn{J}{J} data matrices
#' \eqn{\mathbf X_1, \mathbf X_2, \dots, \mathbf X_J}{X1, X2, ..., XJ}
#' that represent \eqn{J}{J} sets of variables
#' observed on the same set of \eqn{n}{n} individuals. These matrices
#' \eqn{\mathbf X_1, \mathbf X_2, \dots, \mathbf X_J}{X1, X2, ..., XJ},
#' called blocks must have the same number of rows, but may (and usually will)
#' have different numbers of columns.
#'
#' RGCCA aims to study the relationships between these \eqn{J}{J} blocks.
#' It constitutes a general framework for many multi-block component methods
#' (see Tenenhaus and Tenenhaus, 2011 ; Tenenhaus et al. 2017). It combines the
#' power of multi-block data analysis  methods (maximization of well identified
#' criteria) and the flexibility of PLS path modeling (the researcher decides
#' which blocks are connected and which are not). Hence, the use of RGCCA
#' requires the construction (user specified) of a design matrix
#' \eqn{\mathbf C}{C} that
#' characterizes the connections between blocks. Elements of the (symmetric)
#' design matrix \eqn{\mathbf C = (c_{jk})}{C = (c_{jk})}
#' are positive (and usually equal to 1 if blocks \eqn{j}{j}
#' and \eqn{k}{k} are connected, and 0 otherwise).
#' The rgcca() function implements
#' a monotone global convergent algorithm: the bounded criteria to be
#' maximized increases at each step of the iterative procedure and hits, at
#' convergence, a stationary point of the RGCCA optimization problem.
#'
#' Moreover,
#' when the tau argument is used, depending on the dimensionality of each block
#' \eqn{\mathbf X_j, j = 1, \ldots, J}{Xj, j = 1, \ldots, J},
#' the primal algorithm (when \eqn{n \geq p_j}{n >= p_j}) or the dual algorithm
#' (when \eqn{n < p_j}{n < p_j}) is used (see Tenenhaus et al. 2015).
#'
#' When sparsity is
#' specified, SGCCA extends RGCCA to address the issue of variable selection
#' (Tenenhaus et al, 2014). Specifically, RGCCA is combined with an L1-penalty
#' that gives rise to Sparse GCCA (SGCCA). The SGCCA algorithm is very similar
#' to the RGCCA algorithm and keeps the same convergence properties (i.e. the
#' bounded criteria to be maximized increases at each step of the iterative
#' procedure and hits at convergence a stationary point).
#'
#' At last, a deflation strategy can be used to compute several block
#' components
#' (specified by ncomp) per block. Within each block, components or weight
#' vectors are guaranteed to be orthogonal. It should be noted that the numbers
#' of components per block can differ from one block to another.
#'
#' The rgcca() function handles missing values (punctual or blockwise missing
#' structure) using the algorithm described in (Tenenhaus et al, 2005).
#'
#' Guidelines describing how to use RGCCA in practice are provided in
#' (Garali et al., 2018).
#' 
#' Additionally, AC-RGCCA extends RGCCA to enable simultaneous handling of 
#' confounding variation (Goujon et al, to be published). When the arguments 
#' confounders and gamma_confounders are used, the function automatically 
#' applies AC-RGCCA. The confounders matrix (or matrices) represents the 
#' unwanted variation whose influence on the components the user wishes to 
#' limit. This is done through the addition of a penalty in RGCCA's 
#' optimization problem. The penalty tuning parameter gamma_confounders allows 
#' to control the strength of the regularization imposed by the penalty. 
#' 
#' As for classical RGCCA, a dual strategy facilitates the application of 
#' AC-RGCCA to "wide" block matrices (\eqn{n \ll p_j}{n << pj}). Also, 
#' components of higher orders are retrieved using the same deflation approach 
#' implemented for RGCCA. AC-RGCCA extends naturally to the supervised case 
#' through the existing framework for RGCCA. In supervised AC-RGCCA, we assume 
#' that the variables from the response block are not affected by the 
#' confounding effects. 
#' 
#' @param blocks A list that contains the \eqn{J} blocks of variables
#' \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}{X1, X2, ..., XJ}.
#' Block \eqn{\mathbf{X}_j}{Xj} is a matrix of dimension
#' \eqn{n \times p_j}{n x p_j} where \eqn{n} is the number of
#' observations and \eqn{p_j} the number of variables. The blocks argument can
#' be also a fitted cval, rgcca or permutation object.
#' @param method A string specifying which multiblock component
#' method to consider. Possible values are found using
#' \link{available_methods}.
#' @param scale A logical value indicating if variables are standardized.
#' @param scale_block A logical value or a string indicating if each block is
#' scaled.
#'
#' If TRUE or "inertia", each block is divided by the squared root of the sum 
#' of eigenvalues of its empirical covariance matrix.
#'
#' If "lambda1", each block is divided by
#' the square root of the highest eigenvalue of its empirical covariance matrix.
#'
#' If standardization is applied (scale = TRUE), the block scaling applies on
#' the standardized blocks.
#' @param connection  A (\eqn{J \times J}{J x J}) symmetric matrix describing
#' the network of connections between blocks (default value: 1-diag(J)).
#' @param scheme A string or a function specifying the scheme function applied
#' to
#' covariance maximization among "horst" (the identity function), "factorial"
#'  (the square function - default value), "centroid" (the absolute value
#'  function). The scheme function can be any continuously differentiable convex
#'  function and it is possible to design explicitly the scheme function
#'  (e.g. function(x) x^4) as argument of the function.  See (Tenenhaus et al,
#'  2017) for details.
#' @param ncomp A numerical value or a vector of length \eqn{J} indicating
#' the number of components per block. If a single value is provided,
#' the same number of components is extracted for every block.
#' @param tau Either a numerical value, a numeric vector of size
#' \eqn{J}{J}, or a
#' numeric matrix of dimension
#' \eqn{\mathrm{max}(\textrm{ncomp}) \times J}{max(ncomp) x J}
#' containing the values of the regularization parameters
#' (default: tau = 1, for each
#' block and each dimension), or a string equal to "optimal".
#' The regularization parameters varies from 0 (maximizing the correlation) to
#' 1 (maximizing the covariance).
#'
#' If tau is a numerical
#' value, tau is identical across all constraints applied to all
#' block weight vectors.
#'
#' If tau is a vector, tau[j] is used for the constraints applied to
#' all the block weight vectors associated to block \eqn{\mathbf X_j}{Xj}.
#'
#' If tau is a matrix, tau[k, j] is associated with the constraints
#' applied to the kth block weight vector corresponding to block
#' \eqn{\mathbf X_j}{Xj}.
#'
#' If tau = "optimal" the regularization
#' parameters are estimated for each block and each dimension using the Schafer
#' and Strimmer (2005) analytical formula. The tau parameters can also be
#' estimated using
#' \link{rgcca_permutation} or \link{rgcca_cv}.
#' @param sparsity Either a numerical value, a numeric vector of
#' size \eqn{J} or a numeric matrix
#' of dimension \eqn{\textrm{max}(\textrm{ncomp}) \times J} encoding the L1
#' constraints applied to the
#' block weight vectors. For block \eqn{j}, the amount of
#' sparsity varies between
#' \eqn{1/\textrm{sqrt}(p_j)} and 1 (larger values of sparsity
#' correspond to less penalization).
#'
#' If sparsity is a numerical value, then sparsity is identical across
#' all constraints applied to all block weight vectors.
#'
#' If sparsity is a vector, sparsity[j] is identical across the constraints
#' applied to the block weight vectors associated to block
#' \eqn{\mathbf X_j}{Xj}:
#' \deqn{
#'    \forall k, \Vert a_{j,k} \Vert_{1} \le \textrm{sparsity}[j] \sqrt{p_j}.
#' }{for all k, ||ajk||1 <= sparsity(j) sqrt(pj).}
#'
#' If sparsity is a matrix, sparsity[k, j] is associated with the constraints
#' applied to the kth block weight vector corresponding to block
#' \eqn{\mathbf X_j}{Xj}:
#' \deqn{
#'    \Vert a_{j,k}\Vert_{1} \le \textrm{sparsity}[k,j] \sqrt{p_j}.
#' }{||ajk||1 <= sparsity(k, j) sqrt(pj).}
#'
#' The sparsity parameter can be estimated by using \link{rgcca_permutation} or
#' \link{rgcca_cv}.
#' @param init A string giving the type of initialization to use in
#' the RGCCA algorithm. It could be either by
#' Singular Value Decompostion ("svd")
#' or by random initialization ("random") (default: "svd").
#' @param bias A logical value for biased (\eqn{1/n}) or unbiased
#' (\eqn{1/(n-1)}) estimator of the variance/covariance
#' (default: bias = TRUE).
#' @param tol The stopping value for the convergence of the algorithm
#' (default: tol = 1e-08).
#' @param response A numerical value giving the position of the response block.
#' When the response argument is filled, the supervised mode is automatically
#' activated.
#' @param superblock A logical value indicating if the
#' superblock option is used.
#' @param confounders Either a matrix with \eqn{n} observations or a list of 
#' \eqn{J} such matrices encoding the unwanted information, whose influence on 
#' the estimated components is to be minimized by AC-RGCCA. If only one matrix 
#' is given, all blocks are assumed to suffer from the same unwanted variation. 
#' On the other hand, users can specify a list of \eqn{J} matrices for specific 
#' confounders for each block. 
#' 
#' For each block \eqn{j=1,\ldots,J}{j=1,...,J}, a confounders kernel matrix 
#' \eqn{\mathbf K_j}{Kj} of dimensions \eqn{n \times n}{n x n} is defined. Its 
#' elements \eqn{[K_j]_{ii'}}{Kj[i,i']} measure a similarity between 
#' observations \eqn{i} and \eqn{i'} in \eqn{\mathbf Z_j}{Zj}. By default, the 
#' linear kernel \eqn{\mathbf K_j = \mathbf Z_j \mathbf Z_j^\top}{Kj=Zj Zj'} is 
#' computed but alternative kernel matrices can be computed by users to 
#' implement more sophisticated penalties. 
#' 
#' Indeed, users can choose to input the confounders matrices with either 
#' \eqn{\mathbf Z_j}{Zj} (\eqn{n \times q_j}{n x qj} matrix with \eqn{q_j}{qj} 
#' confounders observed on the \eqn{n} individuals) or \eqn{\mathbf K_j}{Kj} 
#' (kernel matrix of size \eqn{n \times n}{n x n}). Internal checking 
#' procedures then determines if the input matrices correspond to 
#' \eqn{\mathbf Z_j}{Zj} or \eqn{\mathbf K_j}{Kj} for 
#' \eqn{j=1,\ldots,J}{j=1,...,J} based on the matrix characteristics. 
#' 
#' If \eqn{\mathbf Z_1, \ldots, \mathbf Z_J}{Z1,...,ZJ} are given, the matrices 
#' are centered by columns and optionally scaled depending on the argument 
#' scale, thus imposing the same scaling as for the blocks 
#' \eqn{\mathbf X_1, \ldots, \mathbf X_J}{X1,...,XJ}. A linear kernel is then 
#' applied.
#' 
#' Alternatively, if users wish to implement a specific kernel and input 
#' \eqn{\mathbf K_1, \ldots, \mathbf K_J}{K1,...,KJ}, the following 
#' transformation is applied to center data in the feature space:
#' \deqn{\mathbf K_j^C = \left(\mathbf I_n - \frac{1}{n} \mathbf 1_n 
#' \mathbf 1_n^\top \right) \mathbf K_j \left(\mathbf I_n - \frac{1}{n} 
#' \mathbf 1_n \mathbf 1_n^\top \right),}{KjC = (In - (1/n)1n 1n') Kj (In - 
#' (1/n)1n 1n'),}
#' where \eqn{\mathbf K_j^C}{KjC} is the kernel matrix of the centered data 
#' (Shawe-Taylor and Cristianini, 2004). Also, if scale is TRUE, the 
#' normalization of data in the feature space is performed using the 
#' transformation:
#' \deqn{\mathbf K_j^N = \mathbf D^{-1/2} \mathbf K_j^C \mathbf D^{-1/2},}{KjN =
#' D^{-1/2} KjC D^{-1/2},}
#' with \eqn{\mathbf D}{D} being the diagonal matrix with the inverse of the 
#' norm of centered feature vectors as diagonal elements: \eqn{d_k = 
#' [K_j^C]_{kk}, k=1, \dots, n}{dk = KjC[kk], k=1,...,n} and \eqn{\mathbf K_j^N
#' }{KjN} giving the kernel matrix obtained on centered and normalized data.
#' 
#' If users wish to define confounders matrices 
#' only for certain blocks, \code{NULL} elements can be given in the list at 
#' positions of blocks for which confounder adjustment is not necessary 
#' (equivalent to setting the corresponding \eqn{\gamma_j}{gammaj} elements to 
#' 0 in the vector given to gamma_confounders). 
#' 
#' When superblock = TRUE, the matrix \eqn{\mathbf K_j}{Kj} and the parameter 
#' \eqn{\gamma_j}{gammaj} are assumed to be identical for all blocks including 
#' the superblock. 
#' @param gamma_confounders Either a numerical value, or a numeric vector of 
#' size \eqn{J} used to define the coefficients 
#' \eqn{\gamma_j \ge 0}{gammaj >= 0} for \eqn{j=1, \ldots, J}{j=1,...,J} which 
#' control the strength of the penalization imposed on the confounding 
#' variation. If one scalar is given, the penalty coefficient is assumed to be 
#' the same for all blocks. Otherwise, users can give a vector of length 
#' \eqn{J} with values chosen for each block. 
#' @param NA_method  A string indicating the method used for
#' handling missing values ("na.ignore", "na.omit"). (default: "na.ignore").
#' \itemize{
#' \item "na.omit" corresponds to perform RGCCA on the fully observed
#' observations (observations from which missing values have been removed).
#' \item "na.ignore" corresponds to perform RGCCA algorithm on available
#' data (See Tenenhaus et al, 2005).}
#' @param verbose A logical value indicating if the progress of the
#' algorithm is reported while computing.
#' @param quiet A logical value indicating if some diagnostic messages
#' are reported.
#' @param n_iter_max Integer giving the algorithm's maximum number of
#' iterations.
#' @param comp_orth A logical value indicating if the deflation should lead to
#' orthogonal block components or orthogonal block weight vectors.
#' @param A Deprecated argument, please use blocks instead.
#' @param C Deprecated argument, please use connection instead.
#' @return A fitted rgcca object.
#' @return \item{Y}{A list of \eqn{J} elements. The jth element
#' of the list \eqn{Y}
#' is a matrix that contains the block components for block j.}
#' @return \item{a}{A list of \eqn{J} elements. The jth element
#' of the list \eqn{a}
#' is a matrix that contains the block weight vectors for block j.}
#' @return \item{astar}{A list of \eqn{J} elements. Each column of astar[[j]] is
#' a vector such that Y[[j]] = blocks[[j]] \%*\% astar[[j]].}
#' @return \item{crit}{A list of vector of length max(ncomp). Each vector of
#' the list is related to one specific deflation stage and reports the values
#' of the criterion for this stage across iterations.}
#' @return \item{primal_dual}{A vector of length J. Element \eqn{j}{j}
#' is either "primal" or "dual", depending on whether the primal or dual
#' RGCCA algorithm was used for block \eqn{j}{j}.}
#' @return \item{AVE}{A list of numerical values giving the indicators of model
#' quality based on the Average Variance Explained (AVE): AVE(for each block),
#' AVE(outer model), AVE(inner model).}
#' @return \item{optimal}{A logical value indicating if the Schaffer and
#' Strimmer formula was applied for estimating the optimal tau parameters.}
#' @return \item{opt}{A list containing some options of
#' the fitted RGCCA object.}
#' @return \item{call}{Call of the function.}
#' @return \item{blocks}{A list that contains the \eqn{J}
#' blocks of variables
#' \eqn{\mathbf X_1, \mathbf X_2, \dots, \mathbf X_J}{X1, X2, ..., XJ}.
#' Block \eqn{\mathbf X_j}{Xj} is a matrix of dimension
#' \eqn{n \times p_j}{n x pj} where \eqn{p_j}{pj} is the number of
#' variables in \eqn{\mathbf X_j}{Xj}. These blocks are preprocessed
#' according to the values of
#' scale/scale_block/NA_method.}
#' @references Garali I, Adanyeguh IM, Ichou F, Perlbarg V, Seyer A, Colsch B,
#' Moszer I, Guillemot V, Durr A, Mochel F, Tenenhaus A. (2018) A strategy for
#' multimodal data integration: application to biomarkers identification
#' in spinocerebellar ataxia. Briefings in Bioinformatics. 19(6):1356-1369.
#' @references Tenenhaus M., Tenenhaus A. and Groenen P. J. (2017). Regularized
#' generalized canonical correlation analysis: a framework for sequential
#' multiblock component methods. Psychometrika, 82(3), 737-777.
#' @references Tenenhaus A., Philippe C. and Frouin, V. (2015). Kernel
#' generalized canonical correlation analysis. Computational Statistics and
#' Data Analysis, 90, 114-131.
#' @references Tenenhaus A., Philippe C., Guillemot V., Le Cao K. A., Grill J.
#' and Frouin, V. (2014), Variable selection for generalized canonical
#' correlation analysis, Biostatistics, 15(3), pp. 569-583.
#' @references Tenenhaus A. and Tenenhaus M., (2011). Regularized Generalized
#' Canonical Correlation Analysis, Psychometrika, 76(2), pp 257-284.
#' @references Tenenhaus, M., Vinzi, V. E., Chatelin, Y. M., & Lauro, C. (2005).
#' PLS path modeling. Computational statistics & data analysis, 48(1), 159-205.
#' @references Schafer J. and Strimmer K. (2005). A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics. Statistical Applications in Genetics and Molecular Biology 4:32.
#' @references Arnaud Gloaguen, Vincent Guillemot, Arthur Tenenhaus.
#' An efficient algorithm to satisfy l1 and l2 constraints.
#' 49emes Journees de Statistique, May 2017, Avignon, France. (hal-01630744)
#' @references Goujon E. et al. AC-RGCCA: Adjusting for Confounding Variation 
#' Within Multiblock Data Analysis. To be published.
#' @references Shawe-Taylor J. and Cristianini N, (2004). Kernel Methods for 
#' Pattern Analysis. Cambridge University Press.
#' @examples
#' ####################
#' # Example 1: RGCCA #
#' ####################
#' # Create the dataset
#' data(Russett)
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#'
#' politic <- as.factor(apply(Russett[, 9:11], 1, which.max))
#'
#' # RGCCA with default values : Blocks are fully connected, factorial scheme
#' # tau = 1 for all blocks, one component per block.
#' fit_rgcca <- rgcca(blocks = blocks)
#'
#' print(fit_rgcca)
#'
#' plot(fit_rgcca, type = "weight", block = 1:3)
#'
#' plot(fit_rgcca,
#'   type = "sample", block = 1:2,
#'   comp = rep(1, 2), resp = politic
#' )
#'
#'
#' ############################################
#' # Example 2: RGCCA and multiple components #
#' ############################################
#' #  By default rgcca() returns orthogonal block components.
#' fit_rgcca <- rgcca(blocks,
#'   method = "rgcca",
#'   connection = 1 - diag(3),
#'   superblock = FALSE,
#'   tau = rep(1, 3),
#'   ncomp = c(2, 2, 2),
#'   scheme = "factorial",
#'   comp_orth = TRUE,
#'   verbose = TRUE
#' )
#'
#' print(fit_rgcca)
#'
#' plot(fit_rgcca,
#'   type = "sample", block = 1,
#'   comp = 1:2, resp = politic
#' )
#'
#' plot(fit_rgcca, type = "weight",
#'      block = 1:3, display_order = FALSE)
#'
#' ##############################
#' # Example 3: MCOA with RGCCA #
#' ##############################
#'
#' fit_rgcca <- rgcca(blocks, method = "mcoa", ncomp = 2)
#' print(fit_rgcca)
#'
#' # biplot representation
#' plot(fit_rgcca, type = "biplot", block = 4, resp = politic)
#'
#' \dontrun{
#'   ####################################
#'   # Example 4: RGCCA and permutation #
#'   ####################################
#'
#'   # Tune the model to find the best set of tau parameters.
#'   # By default, blocks are fully connected.
#'
#'   set.seed(27) #favorite number
#'   perm_out <- rgcca_permutation(blocks,
#'     n_cores = 1,
#'     par_type = "tau",
#'     n_perms = 50
#'   )
#'
#'   print(perm_out)
#'   plot(perm_out)
#'
#'   # all the parameters were imported from a fitted permutation object
#'   fit_rgcca <- rgcca(perm_out)
#'   print(fit_rgcca)
#'
#'
#'   #######################################
#'   # Example 5: RGCCA and dual algorithm #
#'   #######################################
#'   # Download the dataset's package at http://biodev.cea.fr/sgcca/ and install
#'   # it from the package archive file.
#'   # You can do it with the following R commands:
#'   if (!("gliomaData" %in% rownames(installed.packages()))) {
#'     destfile <- tempfile()
#'     download.file(
#'       "http://biodev.cea.fr/sgcca/gliomaData_0.4.tar.gz", destfile
#'      )
#'     install.packages(destfile, repos = NULL, type = "source")
#'   }
#'
#'   data("ge_cgh_locIGR", package = "gliomaData")
#'
#'   blocks <- ge_cgh_locIGR$multiblocks
#'   Loc <- factor(ge_cgh_locIGR$y)
#'   levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#'   blocks[[3]] <- Loc
#'   sapply(blocks, NCOL)
#'
#'   # rgcca algorithm using the dual formulation for X1 and X2
#'   # and the dual formulation for X3. X3 is the group coding matrix associated
#'   # with the qualitative variable Loc. This block is considered
#'   # as response block and specified using the argument response.
#'
#'   fit_rgcca <- rgcca(
#'     blocks = blocks,
#'     response = 3,
#'     method = "rgcca",
#'     tau = c(1, 1, 0),
#'     ncomp = 1,
#'     scheme = function(x) x^2, #factorial scheme,
#'     verbose = TRUE,
#'   )
#'
#'   fit_rgcca$primal_dual
#'   print(fit_rgcca)
#'
#'   ###########################################
#'   # Example 6: RGCCA and variable selection #
#'   ###########################################
#'
#'   # Variable selection and RGCCA : the sgcca algorithm
#'   fit_sgcca <- rgcca(
#'     blocks = blocks,
#'     method = "sgcca",
#'     response = 3,
#'     sparsity = c(.071, .2, 1), ncomp = 1,
#'     scheme = "factorial", verbose = TRUE,
#'   )
#'
#'   print(fit_sgcca)
#'
#'
#'   ############################################
#'   #  Example 7: RGCCA, multiple components   #
#'   #  and different penalties per component   #
#'   ############################################
#'
#'   # S/RGCCA algorithm with multiple components and different
#'   # penalties for each components (-> sparsity is a matrix)
#'
#'   fit_rgcca <- rgcca(blocks, response = 3,
#'     tau = matrix(c(.5, .5, 0, 1, 1, 0), nrow = 2, byrow = TRUE),
#'     ncomp = c(2, 2, 1), scheme = "factorial")
#'
#'   print(fit_rgcca)
#'
#'
#'   # the same applies for SGCCA
#'   fit_sgcca <- rgcca(blocks, response = 3,
#'     sparsity = matrix(c(.071, 0.2,  1,
#'                         0.06, 0.15, 1), nrow = 2, byrow = TRUE),
#'     ncomp = c(2, 2, 1), scheme = "factorial")
#'
#'   print(fit_sgcca)
#'
#'   ##################################################
#'   # Example 8: Supervised mode en cross validation #
#'   ##################################################
#'   # Prediction of the location from GE and CGH
#'
#'   # Tune sparsity values based on the cross-validated accuracy.
#'   set.seed(27) #favorite number
#'   cv_out <- rgcca_cv(blocks, response = 3,
#'                      par_type = "sparsity",
#'                      par_length = 10,
#'                      ncomp = 1,
#'                      prediction_model = "lda",
#'                      metric = "Accuracy",
#'                      k = 3, n_run = 5,
#'                      n_cores = 2)
#'   print(cv_out)
#'   plot(cv_out, display_order = TRUE)
#'
#'   # all the parameters were imported from the fitted cval object.
#'   fit_rgcca <- rgcca(cv_out)
#'   print(fit_rgcca)
#' 
#'   #######################################################
#'   # Example 9: Accounting for confounders with AC-RGCCA #
#'   #######################################################
#'   # Load the dataset
#'   data(example_confounders)
#'   
#'   # Apply methods
#'   ## Classical RGCCA
#'   fit_RGCCA <- rgcca(blocks = example_confounders$blocks, tau = 1, 
#'                      ncomp = 2, superblock = T, scale = T)
#'   print(fit_RGCCA)
#'
#'   ## AC-RGCCA scheme A
#'   fit_AC_RGCCA_A <- rgcca(blocks = example_confounders$blocks, tau = 1, 
#'                           ncomp = 2, superblock = T, scale = T, 
#'                           confounders = example_confounders$confounders_Z_A,
#'                           gamma_confounders = 5)
#'   print(fit_AC_RGCCA_A)
#'   
#'   ## AC-RGCCA scheme B
#'   fit_AC_RGCCA_B <- rgcca(blocks = example_confounders$blocks, tau = 1, 
#'                           ncomp = 2, superblock = T, scale = T, 
#'                           confounders = example_confounders$confounders_Z_B,
#'                           gamma_confounders = 5)
#'   print(fit_AC_RGCCA_B)
#'   
#'   # Compare results on superblock components
#'   plot(fit_RGCCA, type = "samples", block = 3, 
#'        response = example_confounders$dose_groups)
#'   plot(fit_AC_RGCCA_A, type = "samples", block = 3, 
#'        response = example_confounders$dose_groups)
#'   plot(fit_AC_RGCCA_B, type = "samples", block = 3, 
#'        response = example_confounders$dose_groups)
#' }
#'
#' @export
#' @importFrom graphics plot
#' @importFrom stats cor quantile sd p.adjust rnorm pnorm qnorm
#' @importFrom stats model.matrix
#' @importFrom methods is
#' @seealso \code{\link[RGCCA]{plot.rgcca}}, \code{\link[RGCCA]{summary.rgcca}},
#' \code{\link[RGCCA]{rgcca_cv}},
#' \code{\link[RGCCA]{rgcca_permutation}}
#' \code{\link[RGCCA]{rgcca_predict}}
rgcca <- function(blocks, connection = NULL, tau = 1, ncomp = 1,
                  scheme = "factorial", scale = TRUE, init = "svd",
                  bias = TRUE, tol = 1e-08, verbose = FALSE,
                  scale_block = "inertia", method = "rgcca",
                  sparsity = 1, response = NULL,
                  superblock = FALSE, 
                  confounders = NULL, gamma_confounders = 0,
                  NA_method = "na.ignore", quiet = TRUE,
                  n_iter_max = 1000, comp_orth = TRUE,
                  A = NULL, C = NULL) {
  # Check for deprecated arguments
  if (!missing(A)) {
    warning("Argument A is deprecated, use blocks instead.")
    blocks <- A
  }
  if (!missing(C)) {
    warning("Argument C is deprecated, use connection instead.")
    connection <- C
  }
  
  rgcca_args <- as.list(environment())
  ### If specific objects are given for blocks, parameters are imported from
  #   these objects.
  tmp <- get_rgcca_args(blocks, rgcca_args)
  opt <- tmp$opt
  rgcca_args <- tmp$rgcca_args
  rgcca_args$quiet <- quiet
  rgcca_args$verbose <- verbose
  
  blocks <- rgcca_args$blocks
  
  if (opt$disjunction) {
    blocks[[rgcca_args$response]] <- as_disjunctive(
      blocks[[rgcca_args$response]]
    )
  }
  
  ### Apply strategy to deal with NA, scale and prepare superblock
  tmp <- handle_NA(blocks, NA_method = rgcca_args$NA_method)
  na.rm <- tmp$na.rm
  blocks <- scaling(tmp$blocks,
                    scale = rgcca_args$scale,
                    bias = rgcca_args$bias,
                    scale_block = rgcca_args$scale_block,
                    na.rm = na.rm
  )
  if (rgcca_args$superblock) {
    blocks[["superblock"]] <- Reduce(cbind, blocks)
    colnames(blocks[["superblock"]]) <- paste0(
      "s-", colnames(blocks[["superblock"]])
    )
  }
  
  ### Call the gcca function
  gcca_args <- rgcca_args[c(
    "connection", "ncomp", "scheme", "init", "bias", "tol",
    "verbose", "superblock", "response", "n_iter_max", "comp_orth",
    "confounders", "gamma_confounders"
  )]
  gcca_args[["na.rm"]] <- na.rm
  gcca_args[["blocks"]] <- blocks
  gcca_args[["disjunction"]] <- opt$disjunction
  gcca_args[[opt$param]] <- rgcca_args[[opt$param]]
  func_out <- do.call(rgcca_outer_loop, gcca_args)
  
  ### Format the output
  func_out <- format_output(func_out, rgcca_args, opt, blocks)
  
  class(func_out) <- "rgcca"
  invisible(func_out)
}
