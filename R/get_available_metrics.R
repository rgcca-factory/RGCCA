#' Get the metrics computed by the prediction model based on the task.
#' @noRd
get_available_metrics <- function(classification) {
  if (classification) {
    return(c(
      "Accuracy", "Kappa", "F1", "Sensitivity", "Specificity",
      "Pos_Pred_Value", "Neg_Pred_Value", "Precision",
      "Recall", "Detection_Rate", "Balanced_Accuracy"
    ))
  } else {
    return(c("RMSE", "MAE"))
  }
}
