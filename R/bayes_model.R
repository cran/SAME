#Bayesian cox model
main.model = function(data){

  cens = data$survtime
  data$survtime[data$death == 0] = NA # Censored
  cens[data$death == 1] = 0 # Uncensored
  is.censored = as.numeric(is.na(data$survtime))

  modelstring = "
  model{
    for(i in 1:N){

      is.censored[i] ~ dinterval(time[i], cens[i])

      #time[i] ~ dweib(rho, lbd[i])
      time[i] ~  dexp(lbd[i])

      a[i] <- coef[1] + coef[2]*treatment[i,1] + coef[3]*treatment[i,2]+ coef[4]*biomarker[i,1] + coef[5]*biomarker[i,2]
      b[i] <- coef[6]*treatment[i,1]*biomarker[i,1] + coef[8]*treatment[i,2]*biomarker[i,1]
      c[i] <- coef[7]*treatment[i,1]*biomarker[i,2] + coef[9]*treatment[i,2]*biomarker[i,2]

      log(lbd[i]) <- (a[i]+b[i]+c[i])
    }

     for (i in 1:Ncoefs){
       coef[i] ~ dnorm(0,1)
     }

      rho ~ dunif(0,10)
     alpha <- exp(coef[1])


    for(i in 1:4){
      postprob[i,1]<-step(coef[2]+profiles[i,1]*coef[6]+profiles[i,2]*coef[7])
      postprob[i,2]<-step(coef[5]+profiles[i,1]*coef[8]+profiles[i,2]*coef[9])
    }
  }
 "

  data.jags = list(N = nrow(data), time=data$survtime, cens = cens,
                   biomarker = cbind(data$B1, data$B2),
                   treatment = cbind(data$T1, data$T2),
                   profiles  = structure(.Data=c(0,1,0,1,0,0,1,1),.Dim=c(4,2)),
                   is.censored = is.censored,
                   Ncoefs=9)

  initial.jags = function(){
    list(coef=c(0,0,0,0,0,0,0,0,0),
         rho=dunif(1,0,10))
  }

  para.jags = c("coef", "postprob", "rho")

  ph.model = jags.model(textConnection(modelstring), data = data.jags,
                        inits = initial.jags, n.chains = 3)
  update(ph.model,2000)
  res = coda.samples(ph.model, variable.names=para.jags, n.iter=10000,n.thin=10)
  result = as.mcmc(do.call(rbind, res))
  return(summary(result))
}

