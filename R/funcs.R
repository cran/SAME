rowMax=function(mat){
  return(apply(mat,1,max))
}

nlt = function(x, mean, sd){
  pnorm(x, mean, sd, lower.tail=TRUE)
}
ngt = function(x, mean, sd){
  pnorm(x, mean, sd, lower.tail=FALSE)
}

#convert vector of allocation to matrix:
dummy.convert = function(vect,K){
  mat=matrix(0,length(vect),K)
  for(i in 1:length(vect)){
    if(vect[i]==0){
      mat[i,]=c(0,0)}else if(vect[i]==1){
        mat[i,]=c(1,0)}else if(vect[i]==2){
          mat[i,]=c(0,1)}
  }
  return(mat)
}


###data generation function
sim.data = function(treatments, treatments.dummy, profiles, lambda, lfu,
                    accrate, rho){

  N = length(treatments)
  y = rep(NA, N)
  v = runif(N)

  # weibull
  for(i in 1:N){
    #y[i] = rexp(1, rate=lambda[i]) #actual survival time
    y[i] = (-log(v[i])/lambda[i])^(1/rho)
  }

  #nsite = 1
  ## accrual
  arate=accrate#*nsite
  accrtime=N/arate  ##accrual duration
  a=runif(N)*accrtime  ##uniform accrual
  death=rep(1,N)

  if (lfu!=0){
    followup=rexp(N,rate=lfu) ##follow-up time
    death[followup < y]=0  ##censored due to lost of follow-up
    survtime=pmin(followup,y)  ##survival time observed
  } else{
    survtime=y
  }

  asurv=a+survtime
  data=data.frame(y,a,survtime,death,asurv,treatments, treatments.dummy,
                  profiles)

  data=data[order(data$asurv),]
  return(data)
}

######Predictive Probability#####
post_cum_m = function(d,V,x,sigma){

  # A function to get the cumulative posterior distribution of median survival time
  # d:number of events
  n = 1-d
  Z = log(2)*V*(12*x)^(-1/sigma)
  integrand = function(s){tmp = s^(n-2)*exp(-Z/s); res=rep(0,length(tmp)); idx=which(is.finite(tmp)); res[idx]=tmp[idx]; res}

  a1 = integrate(integrand, lower =10e-4, upper = 1)
  a  = a1$value

  output = 12^(-d/sigma)*a*(V*log(2))^d/gamma(d+0.0001)*x^(-d/sigma)

  return(output)
}

close = function(x, value, tol=NULL){
  if(!is.null(tol)){
    fx = x[abs(x-value) <= tol]
  } else {
    fx = x[order(abs(x-value))][1:2]
    sort(match(fx,x),decreasing = FALSE)
  }
}

ints.median = function(N.sampling,d,V,sigma){

  # A function to conduct inversion transform sampling and randomly sample median survival time
  u = runif(N.sampling,0,1)

  # Interpolation
  t  = seq(0.1/12,30/12, by=0.01)

  fx = c()

  for (i in 1:length(t)){
    try(fx[i] <- post_cum_m(d,V,t[i],sigma), silent = TRUE)
  }
  #plot(fx)
  rv.m = c()

  for (i in 1:length(u)){
    indx = close(fx, value=u[i]) #find the index of two values closed to u[i]
    #c = (u[i]-fx[indx[1]])/(fx[indx[2]]-u[i]) #need to deal with if denominator is 0 and sign
    c = abs((u[i]-fx[indx[1]])/(fx[indx[2]]-u[i])) #constant
    rv.m[i] = (c*t[indx[2]] + t[indx[1]])/(c+1)*12
  }

  return(rv.m)
}


PrP = function(N.sampling,d,V,hr0,hr1,sigma,subgroup,N.phase3,lfu,accrate,
               eta,theta,FAtime.phase3){

  # d, V: information obtained from observed data
  # sample median survival time
  # N.sampling: number of sampling
  # N.phase3: sample size for phase III
  # thrsh_sup: threshold for treatment efficacy
  # m: number of events need for phase III

  # metrics
  m0=m1=mean.ratio=sd=post.hr.s2=rep(NA,N.sampling);

  for (i in 1:N.sampling){

    if (N.phase3 == 0){
      m0 = ints.median(N.sampling,d[1],V[1],sigma) #control
      m1 = ints.median(N.sampling,d[2],V[2],sigma) #treatment

      mean0 = mean(m0); mean1 = mean(m1)

      #hazard ratio
      mean.ratio[i] = log(mean0/mean1)
      var.ratio  = var(m0/m1)

      sd[i] = sqrt(log(var.ratio/(mean0/mean1)^2+1))

      post.hr.s2[i] = ngt(theta, mean.ratio[i],sd[i]) + nlt(-theta,mean.ratio[i],sd[i])

    }else{
      lambda.prp = rep(c(hr0,hr1),each = N.phase3)
      allo.ph3.prp.x = rep(c(0,1), each = N.phase3)

      allo.ph3.prp = rep(c(0,subgroup[1]), each = N.phase3)
      allo.ph3.prp.dummy = dummy.convert(allo.ph3.prp,K=2)

      prof.ph3.prp = rep(subgroup[2], 2*N.phase3)
      prof.ph3.prp.dummy = dummy.convert(prof.ph3.prp,K=2)

      data.PrP = sim.data(allo.ph3.prp.x, allo.ph3.prp.dummy, prof.ph3.prp.dummy, lambda.prp,
                          lfu, accrate, rho=1)

      colnames(data.PrP)[c(7:10)] = c("T1","T2","B1","B2")

      #decide the interim analysis time: which is the predefined event size
      edata = data.PrP[data.PrP$death == 1,]

      #FAtime = edata$asurv[m] #time point when observing enough events

      data.PrP = data.PrP[(data.PrP$a < FAtime.phase3),] # only keep obs enrolled before IA
      data.PrP$death.final = data.PrP$death # event at IA
      data.PrP$death.final[data.PrP$asurv > FAtime.phase3] = 0 # censored at IA
      data.PrP$surv = data.PrP$survtime # for subjects have event

      # censoring time truncated for IA
      idx = which(data.PrP$death.final==0)
      data.PrP$surv[idx] = FAtime.phase3 - data.PrP$a[idx]

      V_tilde = d_tilde = V.star = d.star = c()

      for (j in 1:2){
        V_tilde[j] = sum(exp(data.PrP[which(data.PrP$treatments==j-1),]$surv/(sigma*12)))
        d_tilde[j] = sum(data.PrP$death.final[data.PrP$treatments==j-1])

        V.star[j] = V_tilde[j] + V[j]
        d.star[j] = d_tilde[j] + d[j]
      }

      # sample median survival time again with new data set
      m0.updated = ints.median(N.sampling,d.star[1],V.star[1],sigma) #control
      m1.updated = ints.median(N.sampling,d.star[2],V.star[2],sigma) #treatment

      mean0 = mean(m0.updated); mean1 = mean(m1.updated)

      #mean0;mean1
      mean.ratio[i] = log(mean0/mean1)
      var.ratio  = var(m0.updated/m1.updated)

      sd[i] = sqrt(log(var.ratio/(mean0/mean1)^2+1))

      post.hr.s2[i] = ngt(theta, mean.ratio[i], sd[i]) + nlt(-theta, mean.ratio[i], sd[i])
    }
  }

  prp = mean(post.hr.s2 > eta, na.rm = TRUE)

  results = list()
  results$output = cbind(mean.ratio, sd)
  results$prp = ifelse(is.na(prp),0,prp)

  return(results)
}

mams = function(median.c,hr,K,L,lfu,alpha,power,accrate,futility,superiority,theta,
                bio.preva,eta,FAtime.phase3,N.iter){

  #default
  a = b = 1
  N.sampling = 100
  N.phase3.min = 10
  N.phase3.step = 10

  # use hr to get the parameters (beta0/beta/delta)
  beta0 = log(2)/median.c
  beta  = c(0,0)
  gamma = c(0,0)
  delta = log(hr)

  # compute an overall event size based on hazard ratios
  if (all(hr==1)){
    Z_alpha = qnorm(alpha/2,lower.tail = FALSE)
    Z_beta  = qnorm(1-power,lower.tail = FALSE)

    A = (Z_alpha + Z_beta)^2
    B = (log(0.6))^2*0.5^2

    nevents = A/B

    e.cum = ceiling(nevents*L)

    n.cum = 4*e.cum;
    maxN = tail(n.cum,1)

    N.phase3.max = maxN

  }else{
    Z_alpha = qnorm(alpha/2,lower.tail = FALSE)
    Z_beta  = qnorm(1-power,lower.tail = FALSE)

    A = (Z_alpha + Z_beta)^2
    B = (log(min(hr)))^2*0.5^2

    nevents = A/B

    e.cum = ceiling(nevents*L)

    n.cum = 4*e.cum;
    maxN = tail(n.cum,1)

    N.phase3.max = maxN

  }
  #hr = exp(min(delta))


  J = length(L) # number of stages for phase II
  # metrics
  ## phase II
  IAtime = matrix(NA, nrow = N.iter, J); hr.2 = matrix(NA, nrow = N.iter, 4);
  subgrp = matrix(NA,N.iter,2); num_pts = num_evs = matrix(NA, nrow = N.iter, J+1);
  num_evs_sbp = num_pts_sbp = rep(NA, N.iter);
  FAtime = zone = final.dec = rep(NA, N.iter); prp.s2 = rep(0, N.iter)

  ## phase III
  prp.FA = rep(0,N.iter);
  #prp_list = array(NA,dim = c((N.phase3.max-N.phase3.min)/N.phase3.step, 2, N.iter));
  #output.phase2 = array(NA,dim = c(N.sampling,2,N.iter))
  #output.phase3 = array(NA,dim = c(N.sampling*((N.phase3.max-N.phase3.min)/N.phase3.step),3,N.iter))
  N.phase3.prp = matrix(NA, nrow = N.iter, ncol = 2)
  mean.hratio = sd.hr.FA = rep(NA,N.iter)

  #find number of patients to be recruited between each analysis
  n.tilde = c(n.cum[1],n.cum[-1]-n.cum[-length(n.cum)])
  e.tilde = c(e.cum[1],e.cum[-1]-e.cum[-length(e.cum)])

  for (i in 1:N.iter){

    # treatments profile
    profiles=NULL
    for(i1 in 0:1){
      for(i2 in 0:1){
        profiles=rbind(profiles,c(i2,i1))
      }
    }

    #start by determining all patients biomarker profiles:
    all.profiles=matrix(0,maxN,K)
    for(k in 1:maxN){
      all.profiles[k,]=rbinom(K,1,bio.preva)
    }

    #convert to binary score:
    all.score=rep(1,maxN)
    for(k in 1:K){
      all.score = all.score+all.profiles[,k]*2^{k-1}
    }

    #get patient numbers recruited at each stage:
    patientranges=matrix(0,J,2)
    patientranges[1,]=c(1,n.cum[1])
    for(k in 2:J){
      patientranges[k,]=c(n.cum[k-1]+1,n.cum[k])
    }

    #set initial allocation probabilities for linked BAR approach:
    allo.prob = matrix(0,2^K,K+1)
    allo.prob[1,] = rep(1/(K+1),K+1)
    allo.prob[2,] = c(0.5,0.5,0)
    allo.prob[3,] = c(0.5,0,0.5)
    allo.prob[4,] = rep(1/(K+1),K+1)
    pi = allo.prob[,-1] #posterior mean that treatment is better than control


    ##set initial allocation probabilities for regular BAR approach:
    #if(informativeprior==0){
    #  allo.prob=matrix(1/(K+1),2^K,K+1)
    #  pi = allo.prob[,-1]
    #}

    allocation = rep(0,maxN)
    response   = rep(0,maxN)
    survtime   = rep(0,maxN)
    allo.subgroup=matrix(0,2^K,K+1) #generate a matrix to store the randomization result
    diff = rep(0,2^K)

    for(j in 1:J){#for each stage

      #get allocation and response vectors for specific stage j:
      temp.allo = rep(0,n.tilde[j])
      temp.resp = rep(0,n.tilde[j])

      #get binary score for each patient for specific stage j:
      temp.score = all.score[patientranges[j,1]:patientranges[j,2]]
      temp.profiles = all.profiles[patientranges[j,1]:patientranges[j,2],]

      if(j==1){
        for(k in 1:length(temp.score)){
          #based on score to assign group
          temp.allo[k]=which(as.double(rmultinom(1,1,allo.prob[temp.score[k],]))==1)-1
          # count the number of patients in each subgroup
          allo.subgroup[temp.score[k],(temp.allo[k]+1)] = allo.subgroup[temp.score[k],(temp.allo[k]+1)]+1
          diff=rowMax(allo.subgroup[,(2:3)])-allo.subgroup[,1] ### will be used in next stage
        }
      }

      if(j>1){
        for(k in 1:length(temp.score)){
          #allocation probability for experimental groups
          temp.allo.prob=cbind(rep(0,length(profiles[,1])),
                               pi^(a*((patientranges[j,1]+k)/maxN)^b))
          #allocation probability for control group
          temp.allo.prob[,1]=apply(temp.allo.prob,1,max)*(exp(diff)^((1/(K+1))*(patientranges[j,1]+k)/maxN))
          #normalize them
          allo.prob=temp.allo.prob/rowSums(temp.allo.prob)

          #assign patient based on the computed allocation probability
          temp.allo[k]= which(as.double(rmultinom(1,1,allo.prob[temp.score[k],]))==1)-1
          # count number of patients in each subgroup
          allo.subgroup[temp.score[k],(temp.allo[k]+1)]=allo.subgroup[temp.score[k],(temp.allo[k]+1)]+1
          diff=rowMax(allo.subgroup[,(2:3)])-allo.subgroup[,1]
        }
      }

      #get hazard rate for each patient
      temp.allo.dummy = dummy.convert(temp.allo, K=2)

      int = matrix(NA, nrow = length(temp.allo), ncol = 4)

      for (k in 1:length(temp.allo)){
        int[k,] = temp.profiles[k,]%*%t(temp.allo.dummy[k,])
      }

      # lambda: hazard rate for each person
      # beta0 : baseline hazard rate (scale parameter)

      lambda = beta0 * exp(temp.allo.dummy %*% beta + temp.profiles %*% gamma
                           + int%*%delta)

      #lambda = exp(sapply(1:length(temp.allo),function(x){return(beta0 +
      #             ifelse(temp.allo[x]>1,beta[temp.allo[x]-1] +
      #             as.double(delta[temp.allo[x]-1,]%*%temp.profiles[x,]),0) +
      #             temp.profiles[x,]%*%gamma)}))

      #simulate time to event data
      data = sim.data(temp.allo, temp.allo.dummy, temp.profiles, lambda,
                      lfu, accrate, rho=1)

      colnames(data)[c(7:10)] = c("T1","T2","B1","B2")

      ##Fit JAGs model to get new allocation probabilities:
      #first get data for patients who are assessed by time of interim analysis:

      #interim
      #decide the interim analysis time: which is the predefined event size
      data$a = data$a + ifelse(j==1, 0, IAtime[i,j-1])
      data$asurv = data$asurv + ifelse(j==1, 0, IAtime[i,j-1])
      edata.2 = data[data$death == 1,]

      IAtime_pool = c()
      for (k in 1:K){
        IAtime_pool[k] = edata.2[edata.2$treatments==0|edata.2$treatments==k,]$asurv[e.tilde[j]]
      }

      IAtime[i,j] = max(IAtime_pool,na.rm = TRUE) # need to add previous IA time

      data.int = data[(data$a < IAtime[i,j]),] # only keep obs enrolled before IA
      data.int$death.int = data.int$death # event at IA
      data.int$death.int[data.int$asurv > IAtime[i,j]] = 0 # censored at IA
      data.int$surv = data.int$survtime # for subjects have event

      # censoring time truncated for IA
      data.int$surv[data.int$death.int == 0] = IAtime[i,j] - data.int$a[data.int$death.int == 0]

      num_pts[i,j] = dim(data.int)[1]
      num_evs[i,j] = dim(data.int[which(data.int$death.int==1),])[1]

      #record patients info before time truncated
      assign(paste('data.s',j,sep = ''), data.int)

      if (j < J){
        #Bayesian model to compute posterior probability that treatment is better than control
        posterior = main.model(data.int)

        # a matrix with each row corresponding to a row of all profile
        pi = t(matrix(posterior$statistics[10:17,1],K,dim(profiles)[1], byrow = T)) + 0.00001
        #add on a tiny amount to avoid problems with all 0s
      }
    }

    #store information for final analysis
    data.p2 = c()
    for(j in 1:J){
      temp = get(paste('data.s',j,sep = ''))
      data.p2 = rbind(data.p2,temp)
    }


    data.p2$biomarkers = as.factor(ifelse(data.p2$B1 == 0 & data.p2$B2 == 0, 0,
                                          ifelse(data.p2$B1 == 1 & data.p2$B2 == 0, 1,
                                                 ifelse(data.p2$B1 == 0 & data.p2$B2 == 1, 2,
                                                        ifelse(data.p2$B1 == 1 & data.p2$B2 == 1, 3,0)))))


    ###final analysis - for all patients at phase II together
    possibleError = tryCatch(
      coxph(Surv(surv, death.int)~T1:B1+T1:B2+T2:B1+T2:B2, data.p2),
      error = function(e) {
        message(paste("An error occurred for item", i, "coxph() function didnt work:\n"), e)
      }
    )

    if(inherits(possibleError, "error")) next

    itest = coxph(Surv(surv, death.int)~T1:B1+T1:B2+T2:B1+T2:B2, data.p2)

    ms = survfit(Surv(surv,death.int)~treatments,data.p2)

    hr.2[i,]= as.vector(summary(itest)$coef[,'exp(coef)'])

    # rank the hazard ratio and select the smallest one
    subgrp.ind = which(hr.2[i,] == min(hr.2[i,]))
    delta.int = delta[subgrp.ind]
    candidates = matrix(c(1,1,1,2,2,1,2,2),ncol = 2, byrow = TRUE)
    subgrp[i,] = candidates[subgrp.ind,] #treatment and biomaker number

    ## only keep phase II data with selected treatment and biomarker
    data.phase.2 = data.p2[which(data.p2$treatments ==0|
                                   data.p2$treatments == subgrp[i,1]),]
    data.phase.2$treatments = ifelse(data.phase.2$treatments==0,0,1)

    if(subgrp[i,2]==1){
      data.phase.2 = data.phase.2[which(data.phase.2$B1 == 1&data.phase.2$B2 == 0),]
      }else if(subgrp[i,2]==2){
      data.phase.2 = data.phase.2[which(data.phase.2$B1 == 0& data.phase.2$B2 == 1),]
      }

    data.phase.2$profiles = rep(1, dim(data.phase.2)[1])
    #data.phase.2$survtime = data.phase.2$survtime#/12

    #calculate the median survival time
    ms = survfit(Surv(surv,death.int)~treatments,data.phase.2)
    mst0 = summary(ms)$table['treatments=0','median']
    mst1 = summary(ms)$table[paste0('treatments=1'),'median']

    hr0 = log(2)/mst0
    hr1 = log(2)/mst1

    sigma = 1
    V=d=c()
    for (j in 0:1){
      d[j+1] = sum(data.phase.2[which(data.phase.2$treatments==j),]$death.int,na.rm = TRUE)
      V[j+1] = sum(exp(data.phase.2[which(data.phase.2$treatments==j),]$surv/(sigma*12)))
    }


    result = PrP(N.sampling,d,V,hr0,hr1,sigma,subgrp[i,],0,lfu,accrate,eta,
                 theta,FAtime.phase3)
    prp.s2[i] = ifelse(is.na(result$prp),0,result$prp)
    #output.phase2[,,i] = result$output

    tryCatch({
      ### Decision Rule
      if (prp.s2[i] < futility){
        zone[i] = 'futility'
      }else if (prp.s2[i] > superiority){
        zone[i] = 'superiority'
      }else if (prp.s2[i] >= futility & prp.s2[i] <= superiority){
        zone[i] = 'promising'
      }
    }, error = function(e){})

    if(zone[i] == 'futility'| zone[i] == 'superiority'){

      # terminate the trial
      FAtime[i] = IAtime[i,J] #final time is the phase II time
      num_pts[i,4] = 0#dim(data.p2)[1]
      num_evs[i,4] = 0#dim(data.p2[which(data.p2$death.int==1),])[1]

      num_evs_sbp[i] = dim(data.phase.2[which(data.phase.2$death.int==1),])[1]
      num_pts_sbp[i] = dim(data.phase.2)[1]

    }else if (zone[i] == 'promising'){

      #out = info = c()
      prp = prp.s2[i]
      N.phase3 = N.phase3.min

      while (prp < superiority & N.phase3 < N.phase3.max){
        N.phase3 = N.phase3 + N.phase3.step
        result = PrP(N.sampling,d,V,hr0,hr1,sigma,subgrp[i,],N.phase3,lfu,accrate,
                     eta,theta,FAtime.phase3)
        prp = result$prp
        #info = rbind(info, cbind(rep(N.phase3,N.sampling),result$output))
        #colnames(info) = c('sample size', 'hazard ratio', 'sd.hazard ratio')
      }

      #output.phase3[,,i] = info
      N.phase3.prp[i,] = c(N.phase3,prp)

      allo.ph3 = rep(c(0,subgrp[i,1]),each = N.phase3)
      prof.ph3 = rep(1,2*N.phase3)
      allo.ph3.x = rep(c(0,1),each = N.phase3)

      # dummy variable subgrp[i,1]
      allo.ph3.dummy = dummy.convert(allo.ph3, K=2)
      prof.ph3.dummy = dummy.convert(rep(subgrp[i,2],2*N.phase3), K=2)

      lambda.ph3 = beta0 * exp(allo.ph3.x*beta[subgrp[i,1]]  + prof.ph3*gamma[subgrp[i,2]] +
                                 allo.ph3.x*delta.int)

      # continue to enroll patients with estimated event size
      data.3 = sim.data(allo.ph3.x, allo.ph3.dummy, prof.ph3.dummy, lambda.ph3,
                        lfu, accrate, rho=1)
      data.3$a = data.3$a+IAtime[i,J]
      data.3$asurv = data.3$asurv+IAtime[i,J]

      #decide the interim analysis time: which is the predefined event size
      edata.3 = data.3[data.3$death == 1,]

      FAtime[i] = IAtime[i,J] + FAtime.phase3
      #FAtime[i] = edata.3$asurv[m.plan] #time point when observing enough events

      data.p3 = edata.3[(edata.3$a < FAtime[i]),] # only keep obs enrolled before IA
      data.p3$death.int = data.p3$death # event at IA
      data.p3$death.int[data.p3$asurv > FAtime[i]] = 0 # censored at IA
      data.p3$surv = data.p3$survtime # for subjects have event

      # censoring time truncated for IA
      idx = which(data.p3$death.final==0)
      data.p3$surv[idx] =  FAtime.phase3 - data.p3$a[idx]

    }

    if(zone[i] != 'futility'){
      if (zone[i] == 'promising'){
        # re-organize data
        colnames(data.p3)[c(7:10)] = c("T1","T2","B1","B2")
        data.phase.3 = data.p3

        #combine all data from phase II and III to perform final efficacy test
        num_pts[i,4] = sum(num_pts[i,c(1:3)]) + dim(data.phase.3)[1]
        num_evs[i,4] = sum(num_evs[i,c(1:3)]) + dim(data.phase.3[which(data.phase.3$death.int==1),])[1]

        data.phase.3$profiles = rep(1,dim(data.phase.3)[1])
        data.phase.3$biomarkers = data.phase.3$profiles

        data.total = rbind(data.phase.2,data.phase.3)

        num_evs_sbp[i] = dim(data.total[which(data.total$death.int==1),])[1]
        num_pts_sbp[i] = dim(data.total)[1]

        #final efficacy test
        posterior.s3 = main.model.ph3(data.total)
        hr0 = exp(posterior.s3$result$statistics[1,1])
        hr1 = exp(posterior.s3$result$statistics[1,1] + posterior.s3$result$statistics[2,1])
        mean.hratio[i] = hr1/hr0
        sd.hr.FA[i] = posterior.s3$result$statistics[3,'SD']

        prp.FA[i] = ngt(theta,log(mean.hratio[i]),sd.hr.FA[i]) + nlt(-theta,log(mean.hratio[i]),sd.hr.FA[i])

      }else if (zone[i] == 'superiority'){

        #final efficacy test - Bayesian
        data.phase.2$profiles = rep(1,dim(data.phase.2)[1])

        posterior.s3 = main.model.ph3(data.phase.2)
        hr0 = exp(posterior.s3$result$statistics[1,1])
        hr1 = exp(posterior.s3$result$statistics[1,1] + posterior.s3$result$statistics[2,1])
        mean.hratio[i] = hr1/hr0
        sd.hr.FA[i] = posterior.s3$result$statistics[3,'SD']
        prp.FA[i] = ngt(theta,log(mean.hratio[i]),sd.hr.FA[i]) + nlt(-theta,log(mean.hratio[i]),sd.hr.FA[i])
      }
    }
  }

  # Overall summary statistics
  time= round(c(apply(IAtime, FUN = mean,2),mean(FAtime)),digits = 0)
  num.evs = round(apply(num_evs, FUN = mean,2,na.rm=TRUE),digits = 0)
  num.pts =  round(apply(num_pts, FUN = mean,2,na.rm=TRUE),digits = 0)

  design = array(c(time, num.evs, num.pts), dim = c(1,4,3),
                 dimnames = list("Design",
                               c("IA1.ph2","IA2.ph2","FA.ph2","ph3"),
                               c("time","num_evs","num_pts")))

  subgroup.phaseII = c(mean(num_pts_sbp, na.rm = TRUE),
                       mean(num_evs_sbp, na.rm = TRUE))

  names(subgroup.phaseII) = c("num_pts", "num_evs")

  zone.t = round(table(zone)/N.iter*100,digits = 2)

  # performance of arm/biomarker selection
  if (subgrp[i,1]==1 & subgrp[i,2]==1){
    final.dec[i] = 1
  }else if(subgrp[i,1]==1 & subgrp[i,2]==2){
    final.dec[i] = 2
  }else if(subgrp[i,1]==2 & subgrp[i,2]==1){
    final.dec[i] = 3
  }else{
    final.dec[i] = 4
  }
  #arm_count = as.data.frame(table(factor(subgrp[,1], levels = 1:2)))
  #arm_prob  = as.data.frame(round(table(factor(subgrp[,1], levels = 1:2))/N.iter*100.00,digits = 2))

  #tbl.arm = cbind(arm_count, arm_prob[,2])
  #bio_count = as.data.frame(table(factor(subgrp[,2], levels = 1:2)))
  #bio_prob  = as.data.frame(round(table(factor(subgrp[,2], levels = 1:2))/N.iter*100.00,digits = 2))

  #tbl.bio = cbind(bio_count, bio_prob[,2])
  #colnames(tbl.bio)= c("Biomarker","Count","%")

  final.subgroup = as.data.frame(round(table(factor(final.dec, levels = 1:4,labels = c("T1-B1","T1-B2","T2-B1","T2-B2")))/length(final.dec[!is.na(final.dec)])*100.00,digits = 2))
  colnames(final.subgroup)= c("Subgroup","%")


  output = list()
  output$prn = round(mean(prp.FA > eta,na.rm = TRUE),digits = 2)
  output$subgroup = final.subgroup
  output$design = design
  #output$treatment = tbl.arm
  #output$biomarker = tbl.bio
  output$zone = zone.t
  output$hazard.ratio = round(mean(mean.hratio,na.rm = TRUE),digits = 2)
  output$Sel.subgroup = round(subgroup.phaseII, digits = 0)

  return(output)

}
