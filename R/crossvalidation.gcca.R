crossvalidation.gcca <- function(
    rgcca,
    A,
    validation = c("testset", "kfold", "leave-one-out"),
    type = "regression",
    fit = "lm",
    bloc_to_pred = "clinic",
    new_scaled = TRUE,
    k = 5,
    rep = 20,
    nb_cores = parallel::detectCores() - 1) {

    f <- quote(
        predict.gcca(
            rgcca,
            A = lapply(bigA, function(x) x[-inds, ]),
            newA = lapply(bigA, function(x) x[inds, ]),
            type = type,
            fit = fit,
            bloc_to_pred = bloc_to_pred,
            bigA = bigA,
            new_scaled = TRUE
        )
    )

    bigA <- A
    
    if(validation == "leave-one-out")
        v_inds <- 1:nrow(A[[1]])
    if(validation == "kfold"){
        v_inds <- list()
        for (i in 1:rep){
            inds <- sample(nrow(A[[1]]))
            inds <- split(inds, sort(inds %% k))
            v_inds <- c(v_inds, inds)
        }
        
    }

    if(validation == "testset"){
        inds <- sample(nrow(A[[1]]), size = nrow(A[[1]]) * 0.3)
        scores <- list(eval(f))
    }else{
        scores <- parallel::mclapply(v_inds, function(i){
            inds <- i
            eval(f)
        }, mc.cores = nb_cores
        )
    }

    mean(unlist(lapply(scores, function(x) x$score)))

    # # TODO: case for different number of comps between blocks
    # preds <- lapply(1:length(A), 
    #     function(x) sapply(1:max(rgcca$ncomp), 
    #         function(y) mean(abs(sapply(scores, 
    #             function(z) z$pred[[x]][, y])))
    #     )
    # )
    # 
    # return(list(scores, preds))
}
