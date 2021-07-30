shave.matlist <- function(mat_list, nb_cols)
  mapply(function(m, nbcomp) m[, 1:nbcomp, drop = FALSE],
         mat_list, nb_cols,
         SIMPLIFY = FALSE
  )

shave.veclist <- function(vec_list, nb_elts)
  mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY = FALSE)
