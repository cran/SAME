#' Function to estimate the hazard ratios and other statistics of the selected
#' subgroup
#'
#' @description
#' This function is used to estimate the effect size of the selected
#' subgroup.
#'
#' @param data a data.frame in which to interpret the variables named in the formula
#' @param eta a cutoff probability for the strength of evidence for decision-making
#' @param theta a clinically meaningful treatment effect size defined by clinicians
#'
#' @return conduct.phase3()
#'
#' @export
#' @examples
#' \donttest{conduct.phase3(example.2,eta=0.8, theta=0.95)}
#'
#'
conduct.phase3 = function(data,eta,theta){

  # aim to make final decision - eta
  post = main.model.ph3(data)
  h.rate0 = exp(post$result$statistics[1,1])
  h.rate1 = exp(post$result$statistics[1,1] + post$result$statistics[2,1])
  h.ratio = h.rate1/h.rate0
  h.ratio.sd = post$result$statistics[3,'SD']
  prp = ngt(theta,log(h.ratio),h.ratio.sd) + nlt(-theta,log(h.ratio),h.ratio.sd)

  out = list()
  out$post.h.rate = post$post.hrate
  out$final.dec = ifelse(prp > eta,"success","fail")
  out$summary.stat = summary(out$post.h.rate)

  return(out)
  # output: 1.final decision (fail/success)
  #         2.estimate effect size/95% CI (bayes mean/median/)
  #         3.raw effect size/posterior parameter (hr)

}
