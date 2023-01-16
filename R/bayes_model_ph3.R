#Bayesian cox model for phase III

main.model.ph3 = function(data){

  cens = data$survtime
  data$survtime[data$death == 0] = NA # Censored
  cens[data$death == 1] = 0 # Uncensored
  is.censored = as.numeric(is.na(data$survtime))

  modelstring = "
  model{
    for(i in 1:N){

      is.censored[i] ~ dinterval(time[i], cens[i])

      time[i] ~ dweib(rho,lbd[i])
      log(lbd[i]) <- coef[1] + coef[2]*treatment[i]
      #+ coef[3]*biomarker[i]+ coef[4]*treatment[i]*biomarker[i]

    }

     for (i in 1:Ncoefs){
       coef[i] ~ dnorm(0,1)
     }

     rho ~ dunif(0,10)
     alpha <- exp(coef[1])
     hr <- exp(coef[2])

  }
 "

  data.jags = list(N = nrow(data), time=data$survtime, cens = cens,
                   #biomarker = cbind(data$B1,data$B2),
                   #treatment = cbind(data$T1,data$T2),
                   #biomarker = data$profiles,
                   treatment = data$treatments,
                   is.censored = is.censored,
                   Ncoefs=2)

  initial.jags = function(){
    list(coef=rep(0,2),
         rho=dunif(1,0,10))
  }

  para.jags = c("coef", "hr", "rho")

  ph.model = jags.model(textConnection(modelstring), data = data.jags,
                        inits = initial.jags, n.chains = 3)
  update(ph.model,2000)
  res = coda.samples(ph.model, variable.names=para.jags, n.iter=10000,n.thin=10)
  result = summary(as.mcmc(do.call(rbind, res)))

  out = list()
  out$result = result
  out$post.hrate = res[,3][[1]]
  return(out)
}

