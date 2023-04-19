# RGCCA 3.0.0

* Added a `NEWS.md` file to track changes to the package.

Many differences have been introduced since the last version published on CRAN.
We list below the most important ones:
* The function `rgcca` is now the main entry point of the package. 
* Many well-known methods of the multiblock literature are now directly
available by setting the `method` argument in the `rgcca` function.
The list of methods can be found using the `available_methods` function.
* Arguments `A` and `C` of the `rgcca` function are now deprecated. `blocks`
and `connection` should be used instead.
* Small utility functions are no longer exported. 
* Functions `rgcca_transform` and `rgcca_predict` have been added. They allow
applying fitted RGCCA models to project new subjects onto the learned subspaces
and make predictions of response blocks using `caret` predictive models.
* Functions `rgcca_cv` and `rgcca_permutation` have been added. They allow for
finding the best parameters of the models based on a cross-validation or
permutation criterion.
* The function `rgcca_bootstrap` has been added to evaluate the weights of a
model using a bootstrap procedure.
* The function `rgcca_stability` has been added to evaluate the stability of
the variable selection performed by SGCCA.
* Print and plot functions have been added for the outputs of functions
`rgcca`, `rgcca_cv`, `rgcca_permutation`, `rgcca_bootstrap`,
and `rgcca_stability`.
