#' Function to identify the most promising treatment-biomarker-linked subgroup
#'
#' @description
#' This function is used to estimate the effect size of each subgroup and to
#' select the most promising subgroup.
#' @param formula a formula object, with the combinations of treatment and
#' biomarker term, e.g., formula = "T1:B1+T1:B2+T2:B1+T2:B2"
#' @param surv survival time
#' @param event the status indicator, 0=alive, 1=dead
#' @param data a data.frame in which to interpret the variables named in the formula
#'
#' @return conduct.phase2() select the most effective subgroup and returns the
#' estimated hazard ratio.
#'
#' @export
#' @examples
#' conduct.phase2(formula = "T1:B1+T1:B2+T2:B1+T2:B2", surv = "surv",
#' event = "death", data = "example.1")
#'
#'
conduct.phase2 = function(formula, surv, event, data){

  fit   = parse(text=paste0("coxph(Surv(",surv, ",", event,")~",formula, ",",data,")"))
  itest = summary(eval(fit))

  hr = as.vector(itest$coefficients[,2])
  # rank the hazard ratio and select the smallest one
  indx = which(hr == min(hr))

  output = c(rownames(itest$coefficients)[indx], round(c(itest$coefficients[indx,c(2,5)]),digits = 2))
  names(output) = c("Selected Subgroup","hazard ratio", "p-value")

  #output$coef = itest$coefficients[indx]

  return(output)
}
