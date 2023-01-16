#' Function to calibrate the cutoff points under null hypothesis
#'
#' @description
#' This function is used to calibrate the cutoff points under null hypothesis
#' using a multi-arm multi-stage biomarker-enriched design with time-to-event
#' endpoints.
#' @param median.c The median survival time for control group
#' @param K Number of biomarkers
#' @param L Information fraction in terms of the accumulative events in
#' phase II stage, e.g., K = c(1/4,1/2,1)
#' @param lfu Follow-up time
#' @param alpha One-sided familywise error rate
#' @param power Power
#' @param accrate Accrual rate
#' @param theta A clinically meaningful treatment effect size defined by clinicians
#' @param bio.preva Prevalence of biomarker(s)
#' @param FAtime.phase3 the study ending time of phase III
#' @param N.iter Number of iterations
#'
#' @return find.cutoffs() returns the calibrated cutoff points that can control the
#'     type I error rate.
#' @export
#'
#' @examples
#' \donttest{
#' find.cutoffs(median.c=12,K=2,L=c(1/4,1/2,1),lfu=0,alpha=0.05,power=0.9,
#'              accrate=15,theta=log(1.25),bio.preva=c(0.4,0.6),FAtime.phase3=48,
#'              N.iter=3)}


find.cutoffs = function(median.c,K,L,lfu,alpha,power,accrate,theta,bio.preva,
                        FAtime.phase3,N.iter){


  #1st stage parameters calibration
  eta.grid = seq(0.94, 0.96, 0.01)
  j = 0; prn.null = 1;

  while (prn.null >= alpha & j <= length(eta.grid)){

    j = j + 1

    results = mams(median.c,hr=c(1,1,1,1),K,L,lfu,alpha,power,accrate,
                   futility = 0,superiority = 1,theta,bio.preva,eta = eta.grid[j],
                   FAtime.phase3, N.iter = 5)

    prn.null = results$prn
  }

  eta = eta.grid[j]
  out.null.1 = c(eta, prn.null)
  names(out.null.1) = c("eta", "alpha")

  #2nd stage parameters calibration
  fu.grid = 0.2#seq(0.03,0.04,0.01)
  sup.grid = 0.999#seq(0.95,0.95,0.01)

  out.null.2 = c()

  for (i in 1:length(fu.grid)){
    for (j in 1:length(sup.grid)){

      results = mams(median.c,hr=c(1,1,1,1),K,L,lfu,alpha,power,accrate,
                     futility = fu.grid[i],superiority = sup.grid[j],theta,
                     bio.preva,eta = out.null.1[1],FAtime.phase3,N.iter=10)

      out.null.2 = rbind(out.null.2, c(fu.grid[i], sup.grid[j], results$prn))
    }
  }

  # identify the cutoffs
  colnames(out.null.2) = c("futility", "superiority", "alpha")
  out.null.2 = as.data.frame(out.null.2)
  index = min(which(abs(out.null.2$alpha-alpha)==min(abs(out.null.2$alpha-alpha))))

  cutpoints = c(eta, out.null.2$futility, out.null.2$superiority)
  names(cutpoints) = c("eta", "futility", "superiority")

  out = list()

  out$alpha = out.null.2$alpha
  out$cutoff.points = cutpoints

  return(out)

}
