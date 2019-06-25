update_dgx <- function(scheme, Y, dg, n, L, j) { 
	if (mode(scheme) == "function"){
	    assign(formalArgs(scheme), cov2(Y[, j], Y))
	    dgx_tmp  = as.vector(attr(eval(dg), "grad"))
	    dgx      = matrix(rep(dgx_tmp, n), n, L, byrow = TRUE)
    }else{
        switch(scheme,
           ############
           ## HORST  ##
           ############
           "horst"     = {return(matrix(rep(1, n), n, L, byrow = TRUE))},
           ################
           ## FACTORIAL  ##
           ################
           "factorial" = {return(matrix(rep(cov2(Y[, j], Y), n), n, L, byrow = TRUE))},
           ###############
           ## CENTROID  ##
           ###############
           "centroid"  = {return(sign(matrix(rep(cov2(Y[, j], Y), n), n, L, byrow = TRUE)))}
        )
    }
}