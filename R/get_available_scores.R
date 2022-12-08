#' Get the scores computed by the prediction model based on the task.
#' @noRd
get_available_scores <- function(classification) {
  if (classification) {
    return(c(
      "Accuracy", "Kappa", "Mean_F1", "Mean_Sensitivity", "Mean_Specificity",
      "Mean_Pos_Pred_Value", "Mean_Neg_Pred_Value", "Mean_Precision",
      "Mean_Recall", "Mean_Detection_Rate", "Mean_Balanced_Accuracy"
    ))
  } else {
    return(c("RMSE", "MAE"))
  }
}
