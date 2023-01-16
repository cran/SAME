#' Function to simulate Bayesian seamless multi-arm biomarker-enriched phase
#' II/III designs
#'
#' @description
#'     This function finds the required number of events
#'     using a multi-arm multi-stage biomarker-enriched design with time-to-event
#'     endpoints.
#' @param median.c The median survival time for control group
#' @param hr Alternative hazard ratio
#' @param K Number of biomarkers
#' @param L Information fraction in terms of the accumulative events in phase II stage, e.g., K = c(1/4,1/2,1)
#' @param lfu Follow-up time
#' @param alpha One-sided familywise error rate
#' @param power Power
#' @param accrate Accrual rate
#' @param theta A clinically meaningful treatment effect size defined by clinicians
#' @param bio.preva Prevalence of biomarker(s)
#' @param FAtime.phase3 the study ending time of phase III
#' @param N.iter Number of iterations
#'
#'
#'
#' @return sim_trial() returns the nominal type I error rate and calibrated cutoff
#' points, nominal power under user-defined hypothesis, empirical power under
#' user-defined number of simulations, the duration of trial(time), the number
#' of events (num_evs), the number of patients (num_pts) from different stages.
#' The function can also display the number of events and patients under the
#' selected subgroup, the distribution of decision zones and the estimated hazard
#' ratio for the final analysis.
#'
#'
#'
#' @export
#'
#' @import coda
#' @import extraDistr
#' @import ggplot2
#' @import survival
#' @import rjags
#' @import boot
#' @import expint
#' @importFrom utils tail
#' @importFrom stats dunif integrate pnorm qnorm rbinom rexp rmultinom runif update var
#'
#' @examples
#' \donttest{
#' sim.trial(median.c=12,hr=c(1,1,1,0.6),K=2,L=c(1/4,1/2,1),lfu=0,
#'           alpha=0.05,power=0.9,accrate=15,theta=log(1.25),
#'           bio.preva=c(0.4,0.6),FAtime.phase3=48,N.iter=5)
#'}
#'


sim.trial = function(median.c,hr,K,L,lfu,alpha,power,accrate,theta,bio.preva,
                     FAtime.phase3,N.iter){


  #1st stage parameters calibration
  eta.grid = 0.8 #seq(0.2, 0.5, 0.1)
  j = 0; prn.null = 1;

  while (prn.null >= alpha & j <= length(eta.grid)){

    j = j + 1

    results = mams(median.c,hr=c(1,1,1,1),K,L,lfu,alpha,power,accrate,
                   futility = 0,superiority = 1,theta,bio.preva,eta = eta.grid[j],
                   FAtime.phase3, N.iter)

    prn.null = results$prn
  }

  eta = eta.grid[j]
  out.null.1 = c(eta, prn.null)
  names(out.null.1) = c("eta", "alpha")

  #2nd stage parameters calibration
  fu.grid = 0.2#seq(0.03,0.04,0.01)
  sup.grid = 0.8#seq(0.995,0.999,0.001)

  out.null.2 = c()

  for (i in 1:length(fu.grid)){
    for (j in 1:length(sup.grid)){

      results = mams(median.c,hr=c(1,1,1,1),K,L,lfu,alpha,power,accrate,
                     futility = fu.grid[i],superiority = sup.grid[j],theta,
                     bio.preva,eta = out.null.1[1],FAtime.phase3,N.iter)

      out.null.2 = rbind(out.null.2, c(fu.grid[i], sup.grid[j], results$prn))
    }
  }

  # identify the cutoffs
  colnames(out.null.2) = c("futility", "superiority", "alpha")
  out.null.2 = as.data.frame(out.null.2)
  index = min(which(abs(out.null.2$alpha-alpha)==min(abs(out.null.2$alpha-alpha))))

  # using the identified cutoffs to run simulation

  out.alter = mams(median.c,hr,K,L,lfu,alpha,power,accrate,
                   futility = out.null.2[index,1], superiority = out.null.2[index,2], theta,
                   bio.preva,eta = out.null.1[1],FAtime.phase3,N.iter)

  cutpoints = c(eta, out.null.2[index,1], out.null.2[index,2])
  names(cutpoints) = c("eta", "futility", "superiority")

  out = list()
  out$power = out.alter$prn
  out$alpha = out.null.2$alpha[index]
  out$cutoff.points = cutpoints
  out$design = out.alter$design
  out$subgroup.design = out.alter$subgroup.phaseII
  out$zone = out.alter$zone
  out$hazard.ratio = out.alter$hazard.ratio
  out$subgroup  = out.alter$subgroup


  return(out)

}
