#' Function to simulate Bayesian seamless multi-arm biomarker-enriched phase
#' II/III designs with user-defined cutoff points
#'
#' @description
#' This function finds the required number of events using a multi-arm multi-stage
#' biomarker-enriched design with time-to-event endpoints with the user-defined
#' cutoff points.
#' @param median.c The median survival time for control group
#' @param hr Alternative hazard ratio
#' @param K Number of biomarkers
#' @param L Information fraction in terms of the accumulative events in phase II stage, e.g., K = c(1/4,1/2,1)
#' @param lfu Follow-up time
#' @param alpha One-sided family-wise error rate
#' @param power Power
#' @param accrate Accrual rate
#' @param theta A clinically meaningful treatment effect size defined by clinicians
#' @param bio.preva Prevalence of biomarker(s)
#' @param eta A cutoff probability for the strength of evidence for decision-making
#' @param FAtime.phase3 the study ending time of phase III
#' @param eta A cutoff probability for the strength of evidence for decision-making
#' and defined by user.
#' @param futility cutoff point for futility termination
#' @param superiority cutoff point for superiority termination
#' @param N.iter Number of iterations
#'
#' @export
#'
#' @return sim.trial.2() returns the nominal type I error rate,  nominal power under
#'     user-defined hypothesis, empirical power under user-defined number of
#'     simulations, the duration of trial(time), the number of events (num_evs),
#'     the number of patients (num_pts) from different stages. The function can
#'     also display the number of events and patients under the selected subgroup,
#'     the distribution of decision zones and the estimated hazard ratio for the
#'     final analysis.
#'
#' @examples
#' \donttest{
#'  sim.trial.2(median.c=12,hr=c(1,1,1,0.6),K=2,L=c(1/4,1/2,1),lfu=0,alpha=0.05,
#'              power=0.9,accrate=15,theta=log(1.25),bio.preva=c(0.4,0.6),
#'              FAtime.phase3=48,eta=0.2,futility=0.1,superiority=0.9,
#'              N.iter=3)
#'              }


sim.trial.2 = function(median.c,hr,K,L,lfu,alpha,power,accrate,theta,bio.preva,
                     FAtime.phase3,eta,futility,superiority,N.iter){

  # using the use-defined cutoffs to run simulation to obtain type I error rate

  out.null = mams(median.c,hr=c(1,1,1,1),K,L,lfu,alpha,power,accrate,
                 futility,superiority,theta,bio.preva,eta,
                 FAtime.phase3, N.iter)


  # using the use-defined cutoffs to run simulation

  out.alter = mams(median.c,hr,K,L,lfu,alpha,power,accrate,futility, superiority,
                   theta,bio.preva,eta,FAtime.phase3,N.iter)


  out = list()
  out$power = out.alter$prn
  out$alpha = out.null$prn
  out$design = out.alter$design
  out$selected.subgroup = out.alter$subgroup.phaseII
  out$zone = out.alter$zone
  out$hazard.ratio = out.alter$hazard.ratio
  out$subgroup  = out.alter$subgroup

  return(out)

}
