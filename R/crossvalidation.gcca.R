crossvalidation.gcca <- function(
    rgcca,
    A,
    bloc_to_pred,
    validation = c("test", "kfold", "loo"),
    type = "regression",
    fit = "lm",
    new_scaled = TRUE,
    k = 5,
    rep = 10,
    nb_cores = parallel::detectCores() - 1) {

    f <- quote(
        function(){

            Atrain <- lapply(bigA, function(x) x[-inds, ])

            if (class(rgcca) == "sgcca")
                tau <- rgcca$c1
            else
                tau <- rgcca$tau

            rgcca_k <- rgcca.analyze(
                Atrain,
                rgcca$C,
                tau = tau,
                ncomp = rgcca$ncomp,
                scheme = rgcca$scheme,
                scale = FALSE,
                type = class(rgcca),
                verbose = FALSE
            )

            for (k in seq_len(length(bigA))) {
                for (j in seq_len(ncol(rgcca$a[[k]]))) {
                    if (cor(rgcca$a[[k]][, j], rgcca_k$a[[k]][, j]) < 0) 
                        rgcca_k$a[[k]][, j] <- -rgcca_k$a[[k]][, j]
                }
            }
            
            predict.gcca(
                rgcca_k,
                A = Atrain,
                newA = lapply(bigA, function(x) x[inds, ]),
                type = type,
                fit = fit,
                bloc_to_pred = bloc_to_pred,
                bigA = bigA,
                new_scaled = TRUE
            )
        }
    )

    bigA <- A
    
    if(validation == "loo")
        v_inds <- 1:nrow(A[[1]])
    if(validation == "kfold"){
        v_inds <- list()
        for (i in 1:rep){
            inds <- sample(nrow(A[[1]]))
            inds <- split(inds, sort(inds %% k))
            v_inds <- c(v_inds, inds)
        }
    }

    if(validation == "test"){
        inds <- sample(nrow(A[[1]]), size = nrow(A[[1]]) * 0.3)
        scores <- list(eval(f)())
        preds <- scores$res
    }else{
        scores <- parallel::mclapply(
            v_inds, 
            function(i){
                inds <- i
                eval(f)()
            }, mc.cores = nb_cores
        )
    }
    

    if(validation %in% c("loo", "kfold")){

        preds <- lapply(
            1:length(A),
            function (x) Reduce(
                rbind, 
                lapply(
                    scores,
                    function(y) y$pred[[x]]
                    )
                )
            )

        if(validation == "kfold"){
            preds <- lapply(
                1:length(A),
               function (x) Reduce(
                    rbind,
                    lapply(
                        1:nrow(A[[1]]),
                        function(y) apply( preds[[x]][ row.names(preds[[x]]) == rownames(A[[1]])[y], ], 2, mean)
                    )
                )
             )
        }

    for(x in 1:length(preds))
        row.names(preds[[x]]) <- row.names(A[[1]])


    }
    
    scores <- mean(unlist(lapply(scores, function(x) x$score)))

    return(list(scores, preds))
}
