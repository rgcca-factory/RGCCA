
predict.mgcca = function (object, newA, type = c("regression", "classification","coefficient"), fit = c("lm", "cor", "lda", "logistic"), 
                         bloc_to_pred = NULL, y.train = NULL, y.test = NULL, new_scaled = T, scale_size_bloc = T, cutoff = 0.5){
  
  type  = match.arg(type)
  astar = object$astar
  A     = object$A
  p     = sapply(A, ncol)
  L     = length(A)
  
  if (type != "coefficient") {
    if (missing(newA)) {
      pred = object$Y
    }
    else {
      #Arguments checking
      if (missing(fit) || missing(type) || ( missing(bloc_to_pred) && ( missing(y.train) || missing(y.test) ) )){
        stop("Please, define type, fit and bloc_to_pred or y.train/y.test")
      }
      if (type == "classification" && (fit == "cor" || fit == "lm")) stop("Please, classification prediction only works with LDA and LOGISTIC")
      if (type == "regression" && (fit == "lda" || fit == "logistic")) stop("Please, regression prediction only works with LM and COR")
      
      #Compute test parameters
      if (is.null(dim(newA[[1]]))) newp  = sapply(newA, length) else newp  = sapply(newA, ncol)
      
      newL  = length(newA)
      
      #Check similarity between TRAIN and TEST set
      if (is.null(names(A)) ||  is.null(names(newA))) stop("Please, blocs do not have names")
      if (length(names(A)) != L &&  length(names(newA)) != newL) stop("Please, number of blocs and number of blocs named are not the same")
      MATCH = match(names(newA), names(A))
      if (sum(is.na(MATCH)) != 0) stop("Please, blocs in new data did not exist in old data")
      if (!identical(newp, p[MATCH])) stop("The dimension of test dataset is inappropriate!")
      
      #y definition
      if (missing(bloc_to_pred)){
        if (names(y.train) != names(y.test)) stop("Please, train and test sets do not have the same name")
        bloc_to_pred = names(y.train)
      }
      bloc_y    = match(bloc_to_pred, names(A))
      newbloc_y = match(bloc_to_pred, names(newA))
      if(missing(y.train)){
        if (is.na(bloc_y)) stop("Please, bloc_to_pred is absent from training set")
        y.train = A[[bloc_y]]
      }else{
        if (is.na(bloc_y)) warning("Please, bloc_to_pred is absent from training set")
      }
      if(missing(y.test)){
        if (is.na(newbloc_y)) stop("Please, bloc_to_pred is absent from testing set")
        y.test = newA[[newbloc_y]]
      }else{
        if (is.na(newbloc_y)) warning("Please, bloc_to_pred is absent from testing set")
      }
      
      
      # #Scaling
      # if (!new_scaled){
      #   scl_fun = function(data, scaled) {
      #     if (is.null(dim(data))){
      #       if (scale_size_bloc){
      #         scale(t(as.matrix(data)), center = attr(scaled, "scaled:center"), scale = attr(scaled, "scaled:scale")) / sqrt(length(data))
      #       }else{
      #         scale(t(as.matrix(data)), center = attr(scaled, "scaled:center"), scale = attr(scaled, "scaled:scale"))
      #       }
      #     }else{
      #       if (scale_size_bloc){
      #         scale(data, center = attr(scaled, "scaled:center"), scale = attr(scaled, "scaled:scale")) / sqrt(NCOL(data))
      #       }else{
      #         scale(data, center = attr(scaled, "scaled:center"), scale = attr(scaled, "scaled:scale"))
      #       }
      #     }
      #   }
      #   newA = mapply(scl_fun, newA, A, SIMPLIFY=FALSE)
      # }
      
      #Dimension Reduction
      pred = lapply( MATCH, function(u) newA[[u]] %*% astar[[u]])
      
      #Prediction
      if (!is.na(bloc_y) && fit != "cor"){
        comp.Train = Reduce("cbind", object$Y[-bloc_y])
        NAMES      = paste(unlist(mapply(function(name, times){rep(name, times)}, names(A)[-bloc_y], object$ncomp[-bloc_y])), colnames(comp.Train), sep = "_")
      }else{
        comp.Train = Reduce("cbind", object$Y)
        NAMES      = paste(unlist(mapply(function(name, times){rep(name, times)}, names(A), object$ncomp)), colnames(comp.Train), sep = "_")
      } 
      if (!is.na(newbloc_y) && fit != "cor"){
        comp.Test  = Reduce("cbind", pred[-newbloc_y])
        NAMES      = paste(unlist(mapply(function(name, times){rep(name, times)}, names(A)[-bloc_y], object$ncomp[-newbloc_y])), colnames(comp.Train), sep = "_")
      }else{
        comp.Test  = Reduce("cbind", pred)
        NAMES      = paste(unlist(mapply(function(name, times){rep(name, times)}, names(A), object$ncomp)), colnames(comp.Train), sep = "_")
      } 
      comp.Train           = as.data.frame(comp.Train)
      comp.Test            = as.data.frame(comp.Test)
      colnames(comp.Train) = colnames(comp.Test) = NAMES
      if (type == "regression"){
        score = switch(
          fit, 
          "lm"  = {
            reslm  = lm(y.train ~ ., data = comp.Train)
            ychapo = predict(reslm, comp.Test)
            if (any(is.na(ychapo))) warning("NA in predictions.")
            mean((y.test - ychapo)**2)
          },
          "cor" = {
            if (is.null(newA[[1]])) {
              comp.Test
            }else{
              sum((cor(comp.Test)*object$C)[upper.tri(object$C)]**2)
            }
          }
        )
        class.fit = NULL
      } else if (type == "classification"){
        ngroups   = nlevels(as.factor(y.train))
        class.fit = switch(
          fit, 
          "lda"      = {
            reslda     = lda(x = comp.Train, grouping = y.train)
            class.fit  = predict(reslda, comp.Test)$class
          },
          "logistic" = {
            if (ngroups > 2) {
              reslog      = nnet::multinom(y.train ~ ., data = comp.Train, trace = FALSE)
              class.fit   = predict(reslog, newdata = comp.Test)
            }else if (ngroups == 2) {
              reslog      = glm(y.train ~ ., data = comp.Train, family = binomial)
              class.fit   = predict(reslog, type = "response", newdata = comp.Test)
              if (fit.type == "class") {
                class.fit.class = class.fit > cutoff
                class.fit       = factor(as.numeric(class.fit.class))
              }
            }
          }
        )
        score = sum(class.fit == y.test)/length(y.test)
      }
    }
  }
  if (type == "coefficient") {
    pred = astar
  }
  invisible(list(pred = pred, class.fit = class.fit, score = score))
}
