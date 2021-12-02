##########################
#compositeBinary.R
##########################
# INPUT:
###########################

# index: Indices for stratification. individuals in the same strata should be assigned the same unique number. For example, index = c(1,3,2,2,3,3,1) would mean individuals {1,7}, {2,5,6}, and {3,4} are in the same strata.
# R: Binary outcome vector
# Z: Binary vector of indicators of treatment if assuming strong ignorability, or a binary vector of the instrumental variable if an IV study
# D: In an IV study, this is a binary vector containing the actual treatment of interest
# parameter: Which causal parameter for binary outcomes are we interested in? Options are "risk.difference", "risk.ratio", and "effect.ratio"
# null: Null value for the chosen parameter. Defaults are 0 for the effect ratio and the risk difference, and 1 for the risk ratio
# Gamma: At what levels of the unmeasured confounder do we want to test the null? Tests are always performed at Gamma=1, and can be tested at other particular values of Gamma by passing a vector of values for Gamma > 1 here
# effectDirection: Is there a known direction of effect. Default is "both", meaning that a priori we believe the effect could lie in either direction from person to person. Other options are "nonnegative" (rTij >= rCij for all i,j), and "nonpositive" (rTij <= rCij for all i,j)
# monotonicity: Only affects an IV study, and asks if we can assume monotonicity of D (dTij >= dCij for all i,j). This is also called "no defiers". Default is TRUE
# exclusionRestriction: Only affects an IV study, and asks if we can assume the exclusion restriction holds (dTij =dCij => rTij = rCij forall i,j). Default is TRUE
# alternative: Alternative hypothesis for the hypothesis test. Default is "two.sided". Can also be "greater" or "less"
# alpha: Size of the test. Default is 0.05
# CI : should we construct confidence intervals? Default is true. Note that confidence intervals are constructed using the continuous relaxation for the risk ratio and the effect ratio. This has virtiually no impact on the lengths, but can stop the solver from getting stuck because a given hypothesized value isn't feasible (which can happen as both quantities are ratios that still only take on finite numbers of values)
# confidenceLevel: Confidence Level for the intervals
# sensitivityAnalysis: Should we conduct a sensitivity Analysis? If yes, we find the value of Gamma such that we go from rejecting to failing to reject. Only proceeds if we reject at Gamma = 1
#continuousRelax: Should we use the continuous relaxation for the sensitivity analysis? Default is FALSE, but if the solver gets stuck on a given value of Gamma this will speed things up without sacrificing much in terms of accuracy.

####Note: Gurobi is required for using this script. It is freely available for academic use. Their website has information on downloading it, and for installing the R package.

####Note 2: If your p-value is 0.5 while doing a sensitivity analysis, this means that your p-value is AT LEAST 0.5


#######################
#OUTPUT
#######################
# parameter: which parameter was tested
# estimate: estimated value for that parameter
# null: null hypothesis tested
# alternative: direction for the alternative
# pval: Pvalue for the test at Gamma=1
# confint: confidence interval; null if not created
# confidenceLevel: confidence Level
# Pvalues: Matrix with Gamma values requested, and corresponding p-values at those values of Gamma
# maxGamma: Gamma at which we go from rejecting to failing to reject
# effect Direction, monotonicity, exclusionRestriction: What values did we assume for these?




compositeBinary = function(index, R, Z, D=NULL, parameter = "risk.difference", null,Gamma = 1, effectDirection= "both", monotonicity = T, exclusionRestriction = T, alternative = "two.sided", alpha = 0.05,  CI = T, confidenceLevel = 1-alpha, sensitivityAnalysis = T, continuousRelax = F)
{
  require(Matrix)
  if(missing(null))
  {
    null = 0
    if(parameter == "risk.ratio")
    {
      null = 1
    }
  }
  sort.new = function(x)
  { 
    temp = sort(unique(x))
    new = 1:length(temp)
    ret = rep(0,length(x))
    for(i in new)
    {
      ret[x == temp[i]] = i
    }
    ret
  }
  treatment = Z
  outcome = R
  dose = D
  dosage = dose
  index = sort.new(index)
  if(missing(null))
    outcome = 1*outcome
  if(any(outcome!=0 & outcome!= 1) == T)
  {
    stop("Outcomes (R) must be binary")
  }
  
  if(alternative != "two.sided" & alternative != "greater" & alternative != "less")
  {
    stop("Alternative options are two.sided, greater, or less")
    
  }
  if(parameter != "effect.ratio" & parameter != "risk.difference" & parameter != "risk.ratio")
  {
    stop("Causal parameter options are risk.difference, risk.ratio, and effect.ratio")
  }
  
  if(is.null(D) & parameter == "effect.ratio")
  {
    stop("Without an instrumental variable, the causal parameter cannot be effect.ratio")
  }
  
  param = parameter
  
  ns = table(index)
  ms = table(index[treatment*1==1]) 
  
  
  if(any(ms!=1 & ns-ms!=1))
  {
    stop("Strata must have either one treated and the rest controls, or one control and the rest treateds")
  }
  
  
  treatment = (1*treatment == 1)
  
  if(any(treatment != 1 & treatment != 0) & (any(D!=Z)))
  {
    stop("Instrumental Variable (Z) Must be Binary")
  }
  treatment = 1*treatment
  if(any(treatment != 1 & treatment != 0))
  {
    stop("Treatment Vector (Z) Must be Binary")
  }
  dosage = 1*dosage
  if(any(dosage != 1 & dosage != 0))
  {
    stop("D Must be Binary")
  }
  
  
  ER = exclusionRestriction
  MO = monotonicity
  DE = effectDirection
  if(!is.logical(exclusionRestriction))
  {
    stop("exclusionRestriction must be TRUE or FALSE")
  }
  if(!is.logical(monotonicity))
  {
    stop("monotonicity must be TRUE or FALSE")
  }
  
  
  if(DE != "both" & DE != "nonnegative" & DE!= "nonpositive")
  {
    stop("effectDirection can be both, nonnegative, or nonpositive")
  }
  
  gur = suppressWarnings(require(gurobi))
  if(gur == F)
  {
    stop("Gurobi Optimization Suite required to use this script. Gurobi is freely available for academic use.")
  }
  
  alpha.temp = (1 - confidenceLevel)/2
  Gamma.vec = Gamma
  pvalsens = rep(0, length(Gamma[Gamma>1]))
  estimate = 0
  pval = 0
  Gammachange = NULL
  confint = NULL
  if(parameter == "risk.difference")
  {
    res = CRDbinary(index, treatment, outcome, null, DE = DE, alternative = alternative, conf.int = CI, conf.level = confidenceLevel, continuous.relax = continuousRelax)
    deltahat = res$ATE.est
    estimate = deltahat
    pval = res$pval
    confint = res$CI
    
    Gammachange = NULL
    if(sensitivityAnalysis == T)
    {
      
      Gammachange = 1
      if(pval < alpha)
      {
        Gammachange = uniroot(sensCRD2, c(1,1.05), index=index, treatment=treatment, outcome=outcome, null=null, DE=DE, alternative = alternative, alpha = alpha, continuous.relax = T, extendInt = "downX")$root
        if(continuousRelax == F)
        {
          Gammachange = uniroot(sensCRD2, c(Gammachange,Gammachange+.05), index=index, treatment=treatment, outcome=outcome, null=null, DE=DE, alternative = alternative, alpha = alpha, continuous.relax = F, extendInt = "downX")$root
        }
      }
    }
    Gamma = Gamma.vec	
    Gamma = Gamma[Gamma>1]
    if(length(Gamma) > 0)
    {
      for(i in 1:length(Gamma))
      {
        pvalsens[i] = sensCRD(index, treatment, outcome, null, DE, alternative, alpha, Gamma.vec = Gamma[i], calculate.pval = T, continuous.relax = continuousRelax)$pval
      }
    }		
  }
  
  
  if(parameter == "risk.ratio")
  {
    res = CRRbinary(index, treatment, outcome, null, DE = DE, alternative = alternative, continuous.relax = continuousRelax)
    phihat = res$phihat
    estimate = phihat
    pval = res$pval
    
    if(CI == T)
    {
      
      alpha.temp = (1-confidenceLevel)/2
      
      SE.wald =  CRRbinary(index, treatment, outcome, phihat, DE = DE, alternative = alternative, continuous.relax = T)$SE
      w1.per.strat = (ns)*(tapply((treatment)*outcome, index, sum)/ms)
      w2.per.strat =  (ns)*tapply((1-treatment)*null*outcome, index, sum)/(ns-ms)
      
      
      
      ub = round((sum(w1.per.strat) - qnorm(alpha.temp)*SE.wald)/sum(w2.per.strat), 2)
      lb = round((sum(w1.per.strat) - qnorm(1-alpha.temp)*SE.wald)/sum(w2.per.strat), 2)
      if(DE == "nonpositive")
      {
        ub = min(1, ub)
      }
      if(DE == "nonnegative")
      {
        lb = max(1, lb)
      }
      upval = CRRbinary(index, treatment, outcome, ub, DE = DE, alternative = "less", continuous.relax = T)$pval
      RU = (upval < alpha.temp)
      diff = 10
      DU = (-0.01)*RU + 0.01*(!RU)
      ex = (!(RU))*.01
      if(DE == "nonpositive" & RU == F & ub == 1)
      {
        ex = 0
        diff = 1
      }
      while(diff != 1)
      {
        ub = ub + DU
        upval = CRRbinary(index, treatment, outcome, ub, DE = DE, alternative = "less", continuous.relax = T)$pval
        RU1 = (upval < alpha.temp)
        diff = abs(RU1 - RU)	
        RU = RU1  
        if(DE == "nonpositive" & ub == 1 & RU == F)
        {
          ex = 0
          diff = 1
        }
      }
      ub = ub - ex

      loval = CRRbinary(index, treatment, outcome, lb, DE = DE, alternative = "greater", continuous.relax = T)$pval
      
      RU = (loval < alpha.temp)
      diff = 10
      DU = (.01)*RU + (-.01)*(!RU)
      ex = (!(RU))*.01
      if(DE == "nonnegative" & RU == F&lb == 1)
      {
        ex = 0
        diff = 1
      }
      while(diff != 1)
      {
        lb = lb + DU
        loval = CRRbinary(index, treatment, outcome, lb, DE = DE, alternative = "greater", continuous.relax = T)$pval
        RU1 = (loval < alpha.temp)
        diff = abs(RU1 - RU)	
        RU = RU1  	
        if(DE == "nonnegative" & lb==1 & RU == F)
        {
          ex = 0
          diff = 1
        }
      }
      lb = lb + ex
      
      confint = c(lb, ub) 				
    }
    
    Gammachange = NULL
    if(sensitivityAnalysis == T)
    {
      Gammachange = 1
      if(pval < alpha)
      {
        Gammachange = uniroot(sensCRR2, c(1,1.05), index=index, treatment=treatment, outcome=outcome, null=null, DE=DE, alternative = alternative, alpha = alpha, continuous.relax = T, extendInt = "downX")$root
        if(continuousRelax == F)
        {
          Gammachange = uniroot(sensCRR2, c(Gammachange,Gammachange+.05), index=index, treatment=treatment, outcome=outcome, null=null, DE=DE, alternative = alternative, alpha = alpha, continuous.relax = F, extendInt = "downX")$root
        }
      }
    }
    Gamma= Gamma.vec
    Gamma = Gamma[Gamma>1]
    if(length(Gamma) > 0)
    {
      for(i in 1:length(Gamma))
      {
        pvalsens[i] = sensCRR(index, treatment, outcome, null, DE, alternative, alpha, Gamma.vec = Gamma[i], calculate.pval = T, continuous.relax = continuousRelax)$pval
      }
    }
  }
  
  
  
  if(parameter == "effect.ratio")
  {
    res = ERbinary(index, treatment, outcome, dose = dose, null, DE = DE, MO = MO, ER = ER, alternative = alternative, continuous.relax = continuousRelax)
    lambdahat = res$lambdahat
    estimate = lambdahat
    pval = res$pval
    
    if(CI == T)
    {
      
      alpha.temp = (1-confidenceLevel)/2
      
      SE.wald =  ERbinary(index, treatment, outcome, dose, lambdahat, DE = DE, MO = MO, ER = ER, alternative = alternative, continuous.relax = T)$SE
      wATE.per.strat = (ns)*(tapply((treatment)*outcome, index, sum)/ms - tapply((1-treatment)*outcome, index, sum)/(ns-ms))
      ATE.est = sum(wATE.per.strat)
      wATE2.per.strat = (ns)*(tapply((treatment)*dose, index, sum)/ms - tapply((1-treatment)*dose, index, sum)/(ns-ms))
      ATE2.est = sum(wATE2.per.strat)
      
      ub= (ATE.est + qnorm(1-alpha.temp/2)*SE.wald)/ATE2.est
      lb = (ATE.est + qnorm(alpha.temp/2)*SE.wald)/ATE2.est  
      if(MO == T & ER == T)
      {
        ub = min(c(ub, 1))
        lb = max(c(lb, -1))
      }
      if(DE == "nonpositive")
      {
        ub = min(0, ub)
      }
      if(DE == "nonnegative")
      {
        lb = max(0, lb)
      }
      upval = ERbinary(index, treatment, outcome, dose, ub, DE = DE, MO = MO, ER = ER, alternative = "less", continuous.relax = T)$pval
      RU = (upval < alpha.temp)
      diff = 10
      DU = (-0.01)*RU + 0.01*(!RU)
      ex = (!(RU))*.01
      if(MO == T & ER == T & RU == F & ub == 1)
      {
        ex = 0
        diff = 1
      }
      if(DE == "nonpositive" & RU == F&lb == 0)
      {
        ex = 0
        diff = 1
      }
      while(diff != 1)
      {
        ub = ub + DU
        upval = ERbinary(index, treatment, outcome, dose, ub, DE = DE, MO = MO, ER = ER, alternative = "less", continuous.relax = T)$pval		
        RU1 = (upval < alpha.temp)
        diff = abs(RU1 - RU)	
        RU = RU1
        if(MO == T & ER == T  &RU == F& ub == 1)
        {
          ex = 0
          diff = 1
        }
        if(DE == "nonpositive" & RU == F&lb == 0)
        {
          ex = 0
          diff = 1
        }
      }
      ub = ub - ex

      loval = ERbinary(index, treatment, outcome, dose, lb, DE = DE, MO = MO, ER = ER, alternative = "greater", continuous.relax = T)$pval
      
      RU = (loval < alpha.temp)
      diff = 10
      DU = (.01)*RU + (-.01)*(!RU)
      ex = (!(RU))*.01
      if(MO == T & ER == T & RU == F & lb == -1)
      {
        ex = 0
        diff = 1
      }
      if(DE == "nonnegative" & RU == F&lb == 0)
      {
        ex = 0
        diff = 1
      }
      while(diff != 1)
      {
        lb = lb + DU
        loval = ERbinary(index, treatment, outcome, dose, lb, DE = DE, MO = MO, ER = ER, alternative = "greater", continuous.relax = T)$pval
        RU1 = (loval < alpha.temp)
        diff = abs(RU1 - RU)	
        RU = RU1  	
        if(MO == T & ER == T & RU == F & lb == -1)
        {
          ex = 0
          diff = 1
        }
        if(DE == "nonnegative" & RU == F&lb == 0)
        {
          ex = 0
          diff = 1
        }
      }
      lb = lb + ex
      
      confint = c(lb, ub) 				
    }
    
    Gammachange = NULL
    if(sensitivityAnalysis == T)
    {
      Gammachange = 1
      if(pval < alpha)
      {
        Gammachange = uniroot(sensEffectRatio2, c(1,1.1), index=index, treatment=treatment, outcome=outcome, dose = dose, null=null, DE=DE, MO = MO, ER = ER, alternative = alternative, alpha = alpha, continuous.relax = T, extendInt = "downX")$root
        
        if(continuousRelax == F)
        {
          Gammachange = uniroot(sensEffectRatio2, c(Gammachange,Gammachange+0.05), index=index, treatment=treatment, outcome=outcome, dose = dose, null=null, DE=DE, MO = MO, ER = ER, alternative = alternative, alpha = alpha, continuous.relax = F, extendInt = "downX")$root
        }
      }
    }
    
    Gamma = Gamma.vec
    Gamma = Gamma[Gamma>1]
    pvalsens = rep(0, length(Gamma))
    if(length(Gamma) > 0)
    {
      for(i in 1:length(Gamma))
      {
        pvalsens[i] = sensEffectRatio(index, treatment, outcome, dose, null, DE, MO, ER, alternative, alpha, Gamma.vec = Gamma[i], calculate.pval = T, continuous.relax = continuousRelax)$pval
      }
    }
  }
  
  Pvalues = cbind(c(1, Gamma.vec[Gamma.vec>1]), c(pval, pvalsens))
  colnames(Pvalues) = c("Gamma", "P-value")
  if(parameter != "effect.ratio")
  {
    return(list(parameter = parameter, estimate = estimate, null = null, alternative = alternative, pval = pval, confint = confint, confidencelevel = confidenceLevel, Pvalues = Pvalues, maxGamma = Gammachange, effectDirection = DE))
  }
  else
  {
    return(list(parameter = parameter, estimate = estimate, null = null, alternative = alternative, pval = pval, confint = confint, confidencelevel = confidenceLevel, Pvalues = Pvalues, maxGamma = Gammachange, effectDirection = DE, monotonicity = MO, exclusionRestriction = ER))
  }
  
}



##########################
# CRDbinary
################
CRDbinary = function(index, treatment, outcome, null = 0, DE = "both", alternative = "two.sided", conf.int = T, conf.level = .95, continuous.relax = F)
{
  
  require(Matrix)
  
  index = sort.new(index)
  outcome = 1*outcome
  if(any(outcome!=0 & outcome!= 1) == T)
  {
    stop("Outcomes must be binary")
  }
  
  if(alternative != "two.sided" & alternative != "greater" & alternative != "less")
  {
    stop("Alternative options are two.sided, greater, or less")
  }
  
  
  ns = table(index)
  ms = table(index[treatment*1==1]) 
  N.total = sum(ns)
  
  nostratum = length(unique(index))
  treatment = 1*treatment
  if(any(treatment != 1 & treatment != 0))
  {
    stop("Treatment Vector Must be Binary")
  }
  treatment = (1*treatment == 1)
  
  max.ATE = (sum(outcome[treatment]) + sum(!treatment) - sum(outcome[!treatment]))/length(outcome)
  min.ATE = (sum(outcome[treatment])  - sum(outcome[!treatment]) - sum(treatment))/length(outcome)
  if(null < min.ATE || null > max.ATE)
  {
    stop("Null is infeasible given observed data. No allocation of unobserved potential outcomes could satisfy it.")
  }
  
  gur = suppressWarnings(require(gurobi))
  
  
  wATE.per.strat = (ns)/sum(ns)*(tapply((treatment)*outcome, index, sum)/ms - tapply((1-treatment)*outcome, index, sum)/(ns-ms))
  ATE.est = sum(wATE.per.strat)
  
  ntnc = cbind(ms, ns-ms)
  mins = apply(ntnc, 1, min)
  if(!all(mins == 1))
  {
    stop("Each stratum must have either one treated and many controls or one control and many treateds")
  }
  
  if(abs(null) < 1)
  {
    null = round(null*N.total)
  }
  
  
  max.e = (N.total*ATE.est > null)
  
  vec.11 = tapply(outcome*treatment, index, sum)
  vec.01 = tapply((!treatment)*outcome, index, sum)
  vec.10 = tapply((treatment)*(!outcome), index, sum)
  vec.00 = tapply((!treatment)*(!outcome), index, sum)
  V = cbind(vec.11, vec.10, vec.00, vec.01)
  id = apply(V, 1, paste, collapse = "-")
  
  num.id =  sort.new(xtfrm(id))
  nosymm = length(unique(num.id))
  cc = table(num.id)
  bb = tapply(ns, num.id, mean)
  N.sr = sum(cc-1)
  m.01 = vec.01+1
  m.00 = vec.00+1
  m.10 = vec.10 + 1
  m.11 = vec.11 + 1
  
  mult.00 = tapply(m.00, num.id, mean)
  mult.01 = tapply(m.01, num.id, mean)
  mult.11 = tapply(m.11, num.id, mean)
  mult.10 = tapply(m.10, num.id, mean)
  mult.all = mult.11*mult.10*mult.00*mult.01
  if(DE == "nonpositive")
  {
    mult.all = mult.10*mult.01
  }
  if(DE=="nonnegative")
  {
    mult.all = mult.11*mult.00
  }
  ns.type = tapply(ns, num.id, mean)
  N.vars = sum(mult.all)
  index.symm = rep(1:nosymm, mult.all)
  n.types = (mult.all)
  Diff = rep(0, N.vars)
  #V.list = vector("list", N.vars)
  
  PV = rep(0, N.vars)
  n.per = cc
  row.ind = rep(0, 2*N.vars)
  col.ind = row.ind
  values = row.ind
  b = rep(0, nosymm + 1)
  for(kk in 1:nosymm)
  {
    row.ind[which(index.symm==kk)] = rep(kk, n.types[kk])
    col.ind[which(index.symm==kk)] = which(index.symm==kk)  
    values[which(index.symm==kk)] = rep(1, (n.types[kk]))
    b[kk] = n.per[kk]
  }
  row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
  col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
  #row.ind[(2*N.vars+1):(3*N.vars)] = rep(nosymm + 2, N.vars+1)
  
  Gamma.sens = 1
  
  for(kk in 1:nosymm)
  {
    i = which(num.id==kk)[1]
    symmgroup = which(index.symm == kk)
    ind = which(index==i)
    treatstrat = treatment[ind]
    outstrat = outcome[ind]
    outsymm = c(sort(outstrat[treatstrat==F]), sort(outstrat[treatstrat==T]))
    PO = matrix(0, n.types[kk], ns[i])
    count = 1
    mt.01 = mult.01[kk]
    mt.00 = mult.00[kk]
    mt.10 = mult.10[kk]
    mt.11 = mult.11[kk]
    T.00 = matrix(1, mt.00-1, mt.00-1)
    T.00[lower.tri(T.00)] = 0
    T.00 = rbind(T.00, c(rep(0, mt.00-1)) )
    T.01 = matrix(1, mt.01-1, mt.01-1)
    T.01[lower.tri(T.01)] = 0
    T.01 = rbind(T.01, c(rep(0, mt.01-1)) )
    
    T.10 = matrix(1, mt.10-1, mt.10-1)
    T.10[lower.tri(T.10)] = 0
    T.10 = rbind(T.10, c(rep(0, mt.10-1)) )
    T.11 = matrix(1, mt.11-1, mt.11-1)
    T.11[lower.tri(T.11)] = 0
    T.11 = rbind(T.11, c(rep(0, mt.11-1)) )
    c.00 = mt.00
    c.01 = mt.01
    c.10 = mt.10
    c.11 = mt.11
    if(DE == "nonnegative")
    {
      T.10 = matrix(0, 1, mt.10-1)
      c.10 = 1
      T.01 = matrix(1, 1, mt.01-1)
      c.01 = 1
    }
    if(DE == "nonpositive")
    {
      T.11 = matrix(1, 1, mt.11-1)
      c.11 = 1
      T.00 = matrix(0, 1, mt.00-1)
      c.00 = 1
    }
    
    for(ll in 1:c.00)
    {
      for(mm in 1:c.01)
      {
        for(jj in 1:c.10)
        {
          for(oo in 1:c.11)
          {
            
            tempvec = c(T.00[ll,], T.01[mm,], T.10[jj,], T.11[oo,])
            PO[count,] = tempvec[!is.na(tempvec)]
            count = count+1
          }
        }
        
      }
    }
    for(jj in 1:n.types[kk])
    {
      ind.jj = (((jj-1)*(ns[i]-1))+1):(jj*(ns[i]-1))
      po.symm = PO[jj,]
      treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
      outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
      outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm)           
      sum.cont = sum(outcontrol)/(ns[i]-1)
      Q = (outtreat + outcontrol/(ns[i]-1) - sum.cont)*ns[i]
      if(sum(treatstrat)>1)
      {
        sum.cont = sum(outtreat)/(ns[i]-1)
        Q = -(outtreat/(ns[i]-1) + outcontrol - sum.cont)*ns[i]  
      }
      
      qi = Q*max.e - Q*(!max.e)
      ord = order(qi)
      qi.sort = sort(qi)
      
      
      Gamma.sens = 1
      
      
      
      mu = mean(qi.sort)
      sigma2 = mean(qi.sort^2) - mu^2
      
      mu[abs(mu) < 1e-8] = 0
      sigma2[sigma2 < 1e-8] = 0
      
      PV[symmgroup[jj]] = (sigma2)
      Diff[symmgroup[jj]] = sum(outtreat - outcontrol)  
    }
  }
  
  values[(N.vars+1):(2*N.vars)] = Diff
  b[nosymm+1] = 0
  const.dir = c(rep("=", nosymm+1))
  A = sparseMatrix(row.ind, col.ind, x=values)
  
  res = ATEtest(ATE.est, N.total, b, A, const.dir,PV, Diff, null, alternative)
  pval = res$pval
  tstat = res$tstat
  SE = res$SE
  TE.wald = round(ATE.est*N.total)
  if(conf.int == T)
  {
    alpha.temp = (1-conf.level)/2
    
    SE.wald =  ATEtest(ATE.est, N.total, b,A, const.dir,PV, Diff, TE.wald, alternative)$SE
    ub = round(TE.wald - qnorm(alpha.temp)*SE.wald)
    lb = round(TE.wald + qnorm(alpha.temp)*SE.wald)
    if(DE == "nonpositive")
    {
      ub = min(0, ub)
    }
    if(DE == "nonnegative")
    {
      lb = max(0, lb)
    }
    
    upval = ATEtest(ATE.est, N.total, b, A, const.dir,PV, Diff, ub, "less")$pval
    RU = (upval < alpha.temp)
    diff = 10
    DU = (-1)*RU + 1*(!RU)
    ex = !(RU)
    if(DE == "nonpositive" & RU == F & ub == 0)
    {
      diff = 1
      ex = 0
    }
    while(diff != 1)
    {
      ub = ub + DU
      upval = ATEtest(ATE.est, N.total, b, A, const.dir,PV, Diff, ub, "less")$pval
      RU1 = (upval < alpha.temp)
      diff = abs(RU1 - RU)	
      RU = RU1  		
      if(DE == "nonpositive" & ub==0 & RU == F)
      {
        diff=1
        ex = 0
      }
    }
    ub = ub - ex
    
    
    loval = ATEtest(ATE.est, N.total, b, A, const.dir,PV, Diff, lb, "greater")$pval
    RU = (loval < alpha.temp)
    diff = 10
    DU = (1)*RU + -1*(!RU)
    ex = !(RU)
    if(DE == "nonnegative" & RU == F & lb == 0)
    {
      diff = 1
      ex = 0
    }
    while(diff != 1)
    {
      lb = lb + DU
      loval = ATEtest(ATE.est, N.total, b, A, const.dir, PV, Diff, lb, "greater")$pval
      RU1 = (loval < alpha.temp)
      diff = abs(RU1 - RU)
      if(DE == "nonnegative" & lb==0 & RU == F)
      {
        diff=1
        ex = 0
      }
      RU = RU1  		
    }
    lb = lb + ex
    
    CI = c(lb, ub) 	
    return(list(ATE.est = ATE.est, SE = SE/N.total, null = null/N.total, pval = pval, alternative = alternative, CI = CI/N.total, conf.level = conf.level))
  }
  else
  {
    
    return(list(ATE.est = ATE.est, SE = SE/N.total, null = null/N.total, pval = pval, alternative = alternative, CI = NULL, conf.level = NULL))
  }
} 



ATEtest = function(ATE.est, N.total, b, A, const.dir,PV, Diff, null, alternative, continuous.relax = F)
{
  nosymm = length(b)-1
  N.vars = length(PV)
  b[nosymm+1] = null
  
  model = list()
  
  model$A = A  	
  model$obj = c(PV)
  model$sense = const.dir
  model$rhs = b
  model$vtype = c(rep("I", N.vars))
  if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}
  model$modelsense = "max"
  
  solm = gurobi(model, params = list(OutputFlag = 0))
  SE = sqrt(solm$objval)
  
  tstat = (N.total*ATE.est - null)/SE
  
  if(alternative == "two.sided")
  {
    pval = 2*pnorm(-abs(tstat))
  }
  if(alternative == "greater")
  {
    pval = 1 - pnorm((tstat))
  }
  if(alternative == "less")
  {
    pval = pnorm((tstat))
  }
  return(list(tstat = tstat, SE = SE, pval = pval))
}

sort.new = function(x)
{
  temp = sort(unique(x))
  new = 1:length(temp)
  ret = rep(0,length(x))
  for(i in new)
  {
    ret[x == temp[i]] = i
  }
  ret
}	

########
#sensCRD
#########

sensCRD = function(index, treatment, outcome, null = 0, DE = "both", alternative = "two.sided", alpha = 0.05, Gamma.vec = 1, calculate.pval = T, continuous.relax = F)
{
  
  PVAL = calculate.pval
  ns = table(index)
  ms = table(index[treatment*1==1]) 
  N.total = sum(ns)
  nostratum = length(unique(index))
  treatment = 1*treatment
  
  treatment = (1*treatment == 1)
  
  wATE.per.strat = (ns)/sum(ns)*(tapply((treatment)*outcome, index, sum)/ms - tapply((1-treatment)*outcome, index, sum)/(ns-ms))
  ATE.est = sum(wATE.per.strat)
  if(abs(null) < 1)
  {
    null = round(null*N.total)
  }
  max.e = (N.total*ATE.est > null)
  
  vec.11 = tapply(outcome*treatment, index, sum)
  vec.01 = tapply((!treatment)*outcome, index, sum)
  vec.10 = tapply((treatment)*(!outcome), index, sum)
  vec.00 = tapply((!treatment)*(!outcome), index, sum)
  V = cbind(vec.11, vec.10, vec.00, vec.01)
  id = apply(V, 1, paste, collapse = "-")
  
  num.id =  sort.new(xtfrm(id))
  nosymm = length(unique(num.id))
  cc = table(num.id)
  bb = tapply(ns, num.id, mean)
  N.sr = sum(cc-1)
  m.01 = vec.01+1
  m.00 = vec.00+1
  m.10 = vec.10 + 1
  m.11 = vec.11 + 1
  
  mult.00 = tapply(m.00, num.id, mean)
  mult.01 = tapply(m.01, num.id, mean)
  mult.11 = tapply(m.11, num.id, mean)
  mult.10 = tapply(m.10, num.id, mean)
  mult.all = mult.11*mult.10*mult.00*mult.01
  if(DE == "nonpositive")
  {
    mult.all = mult.10*mult.01
  }
  if(DE=="nonnegative")
  {
    mult.all = mult.11*mult.00
  }
  ns.type = tapply(ns, num.id, mean)
  N.vars = sum(mult.all*(ns.type-1))
  index.symm = rep(1:nosymm, mult.all*(ns.type-1))
  n.types = (mult.all)
  Diff = rep(0, N.vars)
  #V.list = vector("list", N.vars)
  
  PM = rep(0, N.vars)
  PV = rep(0, N.vars)
  n.per = cc
  row.ind = rep(0, 3*N.vars+1)
  col.ind = row.ind
  values = row.ind
  b = rep(0, nosymm + 2)
  for(kk in 1:nosymm)
  {
    row.ind[which(index.symm==kk)] = rep(kk, (ns.type[kk]-1)*n.types[kk])
    col.ind[which(index.symm==kk)] = which(index.symm==kk)  
    values[which(index.symm==kk)] = rep(1, (ns.type[kk]-1)*(n.types[kk]))
    b[kk] = n.per[kk]
  }
  row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
  col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
  row.ind[(2*N.vars+1):(3*N.vars+1)] = rep(nosymm + 2, N.vars+1)
  col.ind[(2*N.vars+1):(3*N.vars+1)] = 1:(N.vars+1)
  zscore = rep(0, length(Gamma.vec))
  Rejectvec = zscore
  pvalvec = zscore
  kappavec = zscore
  for(ee in 1:length(Gamma.vec))
  {
    Gamma.sens = Gamma.vec[ee]
    
    for(kk in 1:nosymm)
    {
      i = which(num.id==kk)[1]
      symmgroup = which(index.symm == kk)
      ind = which(index==i)
      treatstrat = treatment[ind]
      outstrat = outcome[ind]
      outsymm = c(sort(outstrat[treatstrat==F]), sort(outstrat[treatstrat==T]))
      PO = matrix(0, n.types[kk], ns[i])
      count = 1
      mt.01 = mult.01[kk]
      mt.00 = mult.00[kk]
      mt.10 = mult.10[kk]
      mt.11 = mult.11[kk]
      T.00 = matrix(1, mt.00-1, mt.00-1)
      T.00[lower.tri(T.00)] = 0
      T.00 = rbind(T.00, c(rep(0, mt.00-1)) )
      T.01 = matrix(1, mt.01-1, mt.01-1)
      T.01[lower.tri(T.01)] = 0
      T.01 = rbind(T.01, c(rep(0, mt.01-1)) )
      
      T.10 = matrix(1, mt.10-1, mt.10-1)
      T.10[lower.tri(T.10)] = 0
      T.10 = rbind(T.10, c(rep(0, mt.10-1)) )
      T.11 = matrix(1, mt.11-1, mt.11-1)
      T.11[lower.tri(T.11)] = 0
      T.11 = rbind(T.11, c(rep(0, mt.11-1)) )
      c.00 = mt.00
      c.01 = mt.01
      c.10 = mt.10
      c.11 = mt.11
      if(DE == "nonnegative")
      {
        T.10 = matrix(0, 1, mt.10-1)
        c.10 = 1
        T.01 = matrix(1, 1, mt.01-1)
        c.01 = 1
      }
      if(DE == "nonpositive")
      {
        T.11 = matrix(1, 1, mt.11-1)
        c.11 = 1
        T.00 = matrix(0, 1, mt.00-1)
        c.00 = 1
      }
      
      
      for(ll in 1:c.00)
      {
        for(mm in 1:c.01)
        {
          for(jj in 1:c.10)
          {
            for(oo in 1:c.11)
            {
              
              tempvec = c(T.00[ll,], T.01[mm,], T.10[jj,], T.11[oo,])
              PO[count,] = tempvec[!is.na(tempvec)]
              count = count+1
            }
          }
          
        }
      }
      for(jj in 1:n.types[kk])
      {
        ind.jj = (((jj-1)*(ns[i]-1))+1):(jj*(ns[i]-1))
        po.symm = PO[jj,]
        treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
        outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
        outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm) 
        sum.cont = sum(outcontrol)/(ns[i]-1)
        Q = (outtreat + outcontrol/(ns[i]-1) - sum.cont)*ns[i]
        if(sum(treatstrat)>1)
        {
          sum.cont = sum(outtreat)/(ns[i]-1)
          Q = -(outtreat/(ns[i]-1) + outcontrol - sum.cont)*ns[i]  
        }
        qi = Q*max.e - Q*(!max.e)
        ord = order(qi)
        qi.sort = sort(qi)
        
        
        mu = rep(0, length(ind)-1)
        sigma2 = rep(0, length(ind)-1)
        
        
        for(j in 1:(length(ind)-1))
        {
          mu[j] = (sum(qi.sort[1:(j)]) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]))/((j) + Gamma.sens*(ns[i]-(j)))
          sigma2[j] = (sum(qi.sort[1:(j)]^2) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]^2))/((j) + Gamma.sens*(ns[i]-(j))) - mu[j]^2
        }
        mu[abs(mu) < 1e-8] = 0
        sigma2[sigma2 < 1e-8] = 0
        
        
        PM[symmgroup[ind.jj]] = mu*(max.e) - mu*(!max.e)
        PV[symmgroup[ind.jj]] = (sigma2)
        Diff[symmgroup[ind.jj]] = sum(outtreat - outcontrol)
        
        
      }
    }
    
    values[(N.vars+1):(2*N.vars)] = Diff
    values[(2*N.vars+1):(3*N.vars+1)] = c(-PM, 1)
    b[nosymm+1] = null
    b[nosymm+2] = 0
    
    const.dir = c(rep("=", nosymm+2))
    model = list()
    alpha.opt = alpha
    if(alternative != "two.sided")
    {
      alpha.opt = 2*alpha
    }
    if(Gamma.sens==1)
    {
      model$A = sparseMatrix(row.ind[1:(2*N.vars)], col.ind[1:(2*N.vars)], x=values[1:(2*N.vars)])
      model$obj = c(PV)
      model$sense = const.dir[1:(nosymm+1)]
      model$rhs = b[1:(nosymm+1)]
      model$vtype = c(rep("I", N.vars))
      if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}
      model$modelsense = "max"
      
      
      solm = gurobi(model, params = list(OutputFlag = 0))
      zed = ((N.total*ATE.est - null)/sqrt(solm$objval))
      tstat = zed
      kappa = (N.total*ATE.est - null)^2 - qchisq(1-alpha.opt, 1)*solm$objval
      kappavec[ee] = kappa
      pval = 0
      if(alternative == "two.sided")
      {
        pval = 2*pnorm(-abs(tstat))
      }
      if(alternative == "greater")
      {
        pval = 1 - pnorm((tstat))
      }
      if(alternative == "less")
      {
        pval = pnorm((tstat))
      }
      Reject = (pval < alpha)
    }
    if(Gamma.sens != 1)
    {
      diff = 200
      kappa = qchisq(1-alpha.opt, 1)
      count=0
      while(diff > 1e-8)
      {
        Plin = -2*N.total*ATE.est*PM - kappa*PV 
        rowind.q =  1:(N.vars+1)
        colind.q = 1:(N.vars+1)
        values.q = c(rep(0, N.vars),1)
        Q = sparseMatrix(rowind.q, colind.q, x=values.q)
        model$A = sparseMatrix(row.ind, col.ind, x=values)
        model$obj = c(Plin,0)
        model$Q = Q
        model$sense = const.dir
        model$rhs = b
        model$vtype = c(rep("I", N.vars), "C")
        if(continuous.relax == T){model$vtype = c(rep("C", N.vars+1))}
        model$lb = c(rep(0, N.vars), -Inf)
        model$modelsense = "min"
        solm = gurobi(model, params = list(OutputFlag = 0))
        x = solm$x[1:N.vars]
        kappa.new = (N.total*ATE.est - sum(PM*x))^2/sum(PV*x)
        diff = abs(kappa.new - kappa)
        pval = 0
        if(PVAL == F)
        {
          diff = 0
          Reject = (kappa.new > kappa)
        }
        kappa = kappa.new
        kappavec[ee] = (N.total*ATE.est - sum(PM*x))^2 - qchisq(1-alpha.opt, 1)*sum(PV*x)
      }
      
      zed = sqrt((N.total*ATE.est - sum(PM*x))^2/sum(PV*x))
      
      if(alternative == "less")
      {
        zed = -zed
      }
      zscore[ee] = zed
      tstat = zed
      
      if(PVAL == T)
      {
        if(alternative == "two.sided")
        {
          pval = 2*pnorm(-abs(tstat))
        }
        if(alternative == "greater")
        {
          pval = 1 - pnorm((tstat))
        }
        if(alternative == "less")
        {
          pval = pnorm((tstat))
        }
        Reject = (pval < alpha)
      }
      
      
      if(sign(-sum(PM*x) + N.total*ATE.est)!=sign(N.total*ATE.est - null))
      {
        Reject = F
        pval = 0.5
        if(alternative == "two.sided")
        {
          pval = 1
        }
      }
      
      if(alternative == "greater" & sum(PM*x) < null)
      {
        pval = .5
      }
      if(alternative == "less" & sum(PM*x) > null)
      {
        pval = .5
      }
      
      
      
    }
    pvalvec[ee] = pval
    Rejectvec[ee] =  Reject
  }
  if(PVAL == F)
  {
    return(list(Gamma.vec = Gamma.vec, Reject = Rejectvec, kappa = kappavec))
  }
  if(PVAL == T)
  {
    return(list(Gamma.vec = Gamma.vec, pval = pvalvec))
  }
}  

###########
#CRRbinary

#######################
CRRbinary = function(index, treatment, outcome, null = 1, DE = "both", alternative = "two.sided", continuous.relax = F)
{
  
  require(Matrix)
  
  
  ns = table(index)
  ms = table(index[treatment*1==1]) 
  N.total = sum(ns)
  
  nostratum = length(unique(index))
  treatment = 1*treatment
  
  treatment = (1*treatment == 1)
  
  max.ATE = (sum(outcome[treatment]) + sum(!treatment) - sum(outcome[!treatment]))/length(outcome)
  min.ATE = (sum(outcome[treatment])  - sum(outcome[!treatment]) - sum(treatment))/length(outcome)
  
  
  gur = suppressWarnings(require(gurobi))
  
  phihat = sum((ns)*(tapply((treatment)*outcome, index, sum)/ms))/sum(ns*tapply((1-treatment)*outcome, index, sum)/(ns-ms))
  phihat = round(phihat, 2)
  null = round(null, 2)
  wCRR.per.strat = (ns)*(tapply((treatment)*outcome, index, sum)/ms - tapply((1-treatment)*null*outcome, index, sum)/(ns-ms))
  CRR.est = sum(wCRR.per.strat)
  
  
  
  
  
  max.e = (CRR.est > 0)
  
  vec.11 = tapply(outcome*treatment, index, sum)
  vec.01 = tapply((!treatment)*outcome, index, sum)
  vec.10 = tapply((treatment)*(!outcome), index, sum)
  vec.00 = tapply((!treatment)*(!outcome), index, sum)
  V = cbind(vec.11, vec.10, vec.00, vec.01)
  id = apply(V, 1, paste, collapse = "-")
  
  num.id =  sort.new(xtfrm(id))
  nosymm = length(unique(num.id))
  cc = table(num.id)
  bb = tapply(ns, num.id, mean)
  N.sr = sum(cc-1)
  m.01 = vec.01+1
  m.00 = vec.00+1
  m.10 = vec.10 + 1
  m.11 = vec.11 + 1
  
  mult.00 = tapply(m.00, num.id, mean)
  mult.01 = tapply(m.01, num.id, mean)
  mult.11 = tapply(m.11, num.id, mean)
  mult.10 = tapply(m.10, num.id, mean)
 
  mult.all = mult.11*mult.10*mult.00*mult.01
  
  if(DE == "nonpositive")
  {
    mult.all = mult.10*mult.01
  }
  if(DE=="nonnegative")
  {
    mult.all = mult.11*mult.00
  }
  ns.type = tapply(ns, num.id, mean)
  N.vars = sum(mult.all)
  index.symm = rep(1:nosymm, mult.all)
  n.types = (mult.all)
  Diff = rep(0, N.vars)
  #V.list = vector("list", N.vars)
  
  PV = rep(0, N.vars)
  n.per = cc
  row.ind = rep(0, 2*N.vars)
  col.ind = row.ind
  values = row.ind
  b = rep(0, nosymm + 1)
  for(kk in 1:nosymm)
  {
    row.ind[which(index.symm==kk)] = rep(kk, n.types[kk])
    col.ind[which(index.symm==kk)] = which(index.symm==kk)  
    values[which(index.symm==kk)] = rep(1, (n.types[kk]))
    b[kk] = n.per[kk]
  }
  row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
  col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
  #row.ind[(2*N.vars+1):(3*N.vars)] = rep(nosymm + 2, N.vars+1)
  
  Gamma.sens = 1
  
  for(kk in 1:nosymm)
  {
    i = which(num.id==kk)[1]
    symmgroup = which(index.symm == kk)
    ind = which(index==i)
    treatstrat = treatment[ind]
    outstrat = outcome[ind]
    outsymm = c(sort(outstrat[treatstrat==F]), sort(outstrat[treatstrat==T]))
    PO = matrix(0, n.types[kk], ns[i])
    count = 1
    mt.01 = mult.01[kk]
    mt.00 = mult.00[kk]
    mt.10 = mult.10[kk]
    mt.11 = mult.11[kk]
    T.00 = matrix(1, mt.00-1, mt.00-1)
    T.00[lower.tri(T.00)] = 0
    T.00 = rbind(T.00, c(rep(0, mt.00-1)) )
    T.01 = matrix(1, mt.01-1, mt.01-1)
    T.01[lower.tri(T.01)] = 0
    T.01 = rbind(T.01, c(rep(0, mt.01-1)) )
    T.10 = matrix(1, mt.10-1, mt.10-1)
    T.10[lower.tri(T.10)] = 0
    T.10 = rbind(T.10, c(rep(0, mt.10-1)) )
    T.11 = matrix(1, mt.11-1, mt.11-1)
    T.11[lower.tri(T.11)] = 0
    T.11 = rbind(T.11, c(rep(0, mt.11-1)) )
    c.00 = mt.00
    c.01 = mt.01
    c.10 = mt.10
    c.11 = mt.11
    if(DE == "nonnegative")
    {
      T.10 = matrix(0, 1, mt.10-1)
      c.10 = 1
      T.01 = matrix(1, 1, mt.01-1)
      c.01 = 1
    }
    if(DE == "nonpositive")
    {
      T.11 = matrix(1, 1, mt.11-1)
      c.11 = 1
      T.00 = matrix(0, 1, mt.00-1)
      c.00 = 1
    }
    
    
    for(ll in 1:c.00)
    {
      for(mm in 1:c.01)
      {
        for(jj in 1:c.10)
        {
          for(oo in 1:c.11)
          {          
            tempvec = c(T.00[ll,], T.01[mm,], T.10[jj,], T.11[oo,])
            PO[count,] = tempvec[!is.na(tempvec)]
            count = count+1
          }
        }
        
      }
    }
    for(jj in 1:n.types[kk])
    {
      ind.jj = (((jj-1)*(ns[i]-1))+1):(jj*(ns[i]-1))
      po.symm = PO[jj,]
      treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
      outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
      outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm)     
      sum.cont = sum(null*outcontrol)/(ns[i]-1)
      Q = (outtreat + null*outcontrol/(ns[i]-1) - sum.cont)*ns[i]
      if(sum(treatstrat)>1)
      {
        sum.cont = sum(outtreat)/(ns[i]-1)
        Q = -(outtreat/(ns[i]-1) + null*outcontrol - sum.cont)*ns[i]  
      }
      
      qi = Q*max.e - Q*(!max.e)
      ord = order(qi)
      qi.sort = sort(qi)
      
      
      Gamma.sens = 1
      
      
      
      mu = mean(qi.sort)
      sigma2 = mean(qi.sort^2) - mu^2
      
      mu[abs(mu) < 1e-8] = 0
      sigma2[sigma2 < 1e-8] = 0
      
      PV[symmgroup[jj]] = (sigma2)
      Diff[symmgroup[jj]] = sum(outtreat - null*outcontrol)  
    }
  }
  values[(N.vars+1):(2*N.vars)] = Diff
  b[nosymm+1] = 0
  const.dir = c(rep("=", nosymm+1))
  A = sparseMatrix(row.ind, col.ind, x=values)
  
  model = list()
  
  model$A = A  	
  model$obj = c(PV)
  model$sense = const.dir
  model$rhs = b
  model$vtype = c(rep("I", N.vars))
  if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}
  model$modelsense = "max"
  
  solm = gurobi(model, params = list(OutputFlag = 0))
  SE = sqrt(solm$objval)
  tstat = (CRR.est)/SE
  
  if(alternative == "two.sided")
  {
    pval = 2*pnorm(-abs(tstat))
  }
  if(alternative == "greater")
  {
    pval = 1 - pnorm((tstat))
  }
  if(alternative == "less")
  {
    pval = pnorm((tstat))
  }
  return(list(tstat = tstat, SE = SE, pval = pval, phihat = phihat))
}
sort.new = function(x)
{
  temp = sort(unique(x))
  new = 1:length(temp)
  ret = rep(0,length(x))
  for(i in new)
  {
    ret[x == temp[i]] = i
  }
  ret
}	

########
#sensCRR
########
sensCRR = function(index, treatment, outcome, null = 1, DE = "both", alternative = "two.sided", alpha = 0.05, Gamma.vec = 1, calculate.pval = T, continuous.relax = F)
{
  PVAL = calculate.pval
  ns = table(index)
  ms = table(index[treatment*1==1]) 
  N.total = sum(ns)
  nostratum = length(unique(index))
  treatment = 1*treatment
  
  treatment = (1*treatment == 1)
  gur = suppressWarnings(require(gurobi))
  phihat = sum((ns)*(tapply((treatment)*outcome, index, sum)/ms))/sum(ns*tapply((1-treatment)*outcome, index, sum)/(ns-ms))
  null = round(null, 2)
  wCRR.per.strat = (ns)*(tapply((treatment)*outcome, index, sum)/ms - tapply((1-treatment)*null*outcome, index, sum)/(ns-ms))
  CRR.est = sum(wCRR.per.strat)
  max.e = (CRR.est > 0)
  
  #Gamma.vec = round(Gamma.vec, 2)
  
  
  vec.11 = tapply(outcome*treatment, index, sum)
  vec.01 = tapply((!treatment)*outcome, index, sum)
  vec.10 = tapply((treatment)*(!outcome), index, sum)
  vec.00 = tapply((!treatment)*(!outcome), index, sum)
  V = cbind(vec.11, vec.10, vec.00, vec.01)
  id = apply(V, 1, paste, collapse = "-")
  
  num.id =  sort.new(xtfrm(id))
  nosymm = length(unique(num.id))
  cc = table(num.id)
  bb = tapply(ns, num.id, mean)
  N.sr = sum(cc-1)
  m.01 = vec.01+1
  m.00 = vec.00+1
  m.10 = vec.10 + 1
  m.11 = vec.11 + 1
  
  mult.00 = tapply(m.00, num.id, mean)
  mult.01 = tapply(m.01, num.id, mean)
  mult.11 = tapply(m.11, num.id, mean)
  mult.10 = tapply(m.10, num.id, mean)
  mult.all = mult.11*mult.10*mult.00*mult.01
  if(DE == "nonpositive")
  {
    mult.all = mult.10*mult.01
  }
  if(DE=="nonnegative")
  {
    mult.all = mult.11*mult.00
  }
  ns.type = tapply(ns, num.id, mean)
  N.vars = sum(mult.all*(ns.type-1))
  index.symm = rep(1:nosymm, mult.all*(ns.type-1))
  n.types = (mult.all)
  Diff = rep(0, N.vars)
  #V.list = vector("list", N.vars)
  row.ind = rep(0, 3*N.vars+1)
  col.ind = row.ind
  values = row.ind
  n.per = cc
  b = rep(0, nosymm + 2)
  for(kk in 1:nosymm)
  {
    row.ind[which(index.symm==kk)] = rep(kk, (ns.type[kk]-1)*n.types[kk])
    col.ind[which(index.symm==kk)] = which(index.symm==kk)  
    values[which(index.symm==kk)] = rep(1, (ns.type[kk]-1)*(n.types[kk]))
    b[kk] = n.per[kk]
  }
  row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
  col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
  row.ind[(2*N.vars+1):(3*N.vars+1)] = rep(nosymm + 2, N.vars+1)
  col.ind[(2*N.vars+1):(3*N.vars+1)] = 1:(N.vars+1)
  PM = rep(0, N.vars)
  PV = rep(0, N.vars)
  
  
  
  #Gamma.vec = sort(c(seq(1, 1.16, by = .04), 1.163, seq(1.06, 1.079, by = .001))) 
  #Gamma.vec = 1
  
  zscore = rep(0, length(Gamma.vec))
  Rejectvec = zscore
  pvalvec = zscore
  kappavec = zscore
  for(ee in 1:length(Gamma.vec))
  {
    Gamma.sens = Gamma.vec[ee]
    
    for(kk in 1:nosymm)
    {
      i = which(num.id==kk)[1]
      symmgroup = which(index.symm == kk)
      ind = which(index==i)
      treatstrat = treatment[ind]
      outstrat = outcome[ind]
      outsymm = c(sort(outstrat[treatstrat==F]), sort(outstrat[treatstrat==T]))
      PO = matrix(0, n.types[kk], ns[i])
      count = 1
      mt.01 = mult.01[kk]
      mt.00 = mult.00[kk]
      mt.10 = mult.10[kk]
      mt.11 = mult.11[kk]
      T.00 = matrix(1, mt.00-1, mt.00-1)
      T.00[lower.tri(T.00)] = 0
      T.00 = rbind(T.00, c(rep(0, mt.00-1)) )
      T.01 = matrix(1, mt.01-1, mt.01-1)
      T.01[lower.tri(T.01)] = 0
      T.01 = rbind(T.01, c(rep(0, mt.01-1)) )
      
      T.10 = matrix(1, mt.10-1, mt.10-1)
      T.10[lower.tri(T.10)] = 0
      T.10 = rbind(T.10, c(rep(0, mt.10-1)) )
      T.11 = matrix(1, mt.11-1, mt.11-1)
      T.11[lower.tri(T.11)] = 0
      T.11 = rbind(T.11, c(rep(0, mt.11-1)) )
      c.00 = mt.00
      c.01 = mt.01
      c.10 = mt.10
      c.11 = mt.11
      if(DE == "nonnegative")
      {
        T.10 = matrix(0, 1, mt.10-1)
        c.10 = 1
        T.01 = matrix(1, 1, mt.01-1)
        c.01 = 1
      }
      if(DE == "nonpositive")
      {
        T.11 = matrix(1, 1, mt.11-1)
        c.11 = 1
        T.00 = matrix(0, 1, mt.00-1)
        c.00 = 1
      }
      
      for(ll in 1:c.00)
      {
        for(mm in 1:c.01)
        {
          for(jj in 1:c.10)
          {
            for(oo in 1:c.11)
            {
              
              tempvec = c(T.00[ll,], T.01[mm,], T.10[jj,], T.11[oo,])
              PO[count,] = tempvec[!is.na(tempvec)]
              count = count+1
            }
          }
          
        }
      }
      for(jj in 1:n.types[kk])
      {
        ind.jj = (((jj-1)*(ns[i]-1))+1):(jj*(ns[i]-1))
        po.symm = PO[jj,]
        treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
        outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
        outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm)
        sum.cont = sum(null*outcontrol)/(ns[i]-1)
        Q = (outtreat + null*outcontrol/(ns[i]-1) - sum.cont)*ns[i]
        if(sum(treatstrat)>1)
        {
          sum.cont = sum(outtreat)/(ns[i]-1)
          Q = -(outtreat/(ns[i]-1) + null*outcontrol - sum.cont)*ns[i]  
        }
        
        
        
        
        qi = Q*max.e - Q*(!max.e)
        ord = order(qi)
        qi.sort = sort(qi)
        
        
        mu = rep(0, length(ind)-1)
        sigma2 = rep(0, length(ind)-1)
        
        
        for(j in 1:(length(ind)-1))
        {
          mu[j] = (sum(qi.sort[1:(j)]) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]))/((j) + Gamma.sens*(ns[i]-(j)))
          sigma2[j] = (sum(qi.sort[1:(j)]^2) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]^2))/((j) + Gamma.sens*(ns[i]-(j))) - mu[j]^2
        }
        mu[abs(mu) < 1e-8] = 0
        sigma2[sigma2 < 1e-8] = 0
        
        
        PM[symmgroup[ind.jj]] = mu*(max.e) - mu*(!max.e)
        PV[symmgroup[ind.jj]] = (sigma2)
        
        Diff[symmgroup[ind.jj]] = sum(outtreat - null*outcontrol)
        
      }
    }
    
    
    
    values[(N.vars+1):(2*N.vars)] = Diff
    values[(2*N.vars+1):(3*N.vars+1)] = c(-PM, 1)
    b[nosymm+1] = 0
    b[nosymm+2] = 0
    alpha.opt = alpha
    if(alternative != "two.sided")
    {
      alpha.opt = 2*alpha
    }
    
    const.dir = c(rep("=", nosymm+2))
    model = list()
    if(Gamma.sens==1)
    {
      model$A = sparseMatrix(row.ind[1:(2*N.vars)], col.ind[1:(2*N.vars)], x=values[1:(2*N.vars)])
      model$obj = c(PV)
      model$sense = const.dir[1:(nosymm+1)]
      model$rhs = b[1:(nosymm+1)]
      model$vtype = c(rep("I", N.vars))
      if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}
      
      
      model$modelsense = "max"
      
      
      solm = gurobi(model, params = list(OutputFlag = 0))
      zed = (CRR.est/sqrt(solm$objval))
      kappavec[ee] = CRR.est^2 - qchisq(1-alpha.opt, 1)*solm$objval
      x = solm$x[1:N.vars]
      tstat = zed
      pval = 0
      if(alternative == "two.sided")
      {
        pval = 2*pnorm(-abs(tstat))
      }
      if(alternative == "greater")
      {
        pval = 1 - pnorm((tstat))
      }
      if(alternative == "less")
      {
        pval = pnorm((tstat))
      }
      Reject = (pval < alpha)
    }
    if(Gamma.sens != 1)
    {
      diff = 200
      kappa = qchisq(1-alpha.opt, 1)
      while(diff > 1e-8)
      {
        Plin = -2*CRR.est*PM - kappa*PV 
        rowind.q =  1:(N.vars+1)
        colind.q = 1:(N.vars+1)
        values.q = c(rep(0, N.vars),1)
        Q = sparseMatrix(rowind.q, colind.q, x=values.q)
        model$A = sparseMatrix(row.ind, col.ind, x=values)
        model$obj = c(Plin,0)
        model$Q = Q
        model$sense = const.dir
        model$rhs = b
        
        
        model$vtype = c(rep("I", N.vars), "C")
        if(continuous.relax == T){model$vtype = c(rep("C", N.vars+1))}
        model$lb = c(rep(0, N.vars), -Inf)
        
        
        model$modelsense = "min"
        
        
        solm = gurobi(model, params = list(OutputFlag = 0))
        x = solm$x[1:N.vars]
        kappa.new = (CRR.est - sum(PM*x))^2/sum(PV*x)
        diff = abs(kappa.new - kappa)
        pval = 0
        if(PVAL == F)
        {
          diff = 0
          Reject = (kappa.new > kappa)
        }
        kappa = kappa.new
        kappavec[ee] = (CRR.est - sum(PM*x))^2 - qchisq(1-alpha.opt, 1)*sum(PV*x)
      }
      zed = sqrt((CRR.est - sum(PM*x))^2/sum(PV*x))
      
      if(alternative == "less")
      {
        zed = -zed
      }
      zscore[ee] = zed
      tstat = zed
      
      if(PVAL == T)
      {
        if(alternative == "two.sided")
        {
          pval = 2*pnorm(-abs(tstat))
        }
        if(alternative == "greater")
        {
          pval = 1 - pnorm((tstat))
        }
        if(alternative == "less")
        {
          pval = pnorm((tstat))
        }
        Reject = (pval < alpha)
      }
      
      
      if(sign(CRR.est - sum(PM*x))!=sign(CRR.est))
      {
        Reject = F
        pval = 0.5
        if(alternative == "two.sided")
        {
          pval = 1
        }
      }
      
      if(alternative == "greater" & sum(PM*x) < 0)
      {
        pval = .5
      }
      if(alternative == "less" & sum(PM*x) > 0)
      {
        pval = .5
      }
      
      
      
    }
    pvalvec[ee] = pval
    Rejectvec[ee] = Reject 
  }
  if(PVAL == F)
  {
    return(list(Gamma.vec = Gamma.vec, Reject = Rejectvec, phihat = phihat, kappa = kappavec))
  }
  if(PVAL == T)
  {
    return(list(Gamma.vec=Gamma.vec, pval = pvalvec, phihat = phihat))
  }
  
  
  
}     


##########
#ERbinary
##########

ERbinary = function(index, treatment, outcome, dose, null=0, DE = "both", MO = T, ER = T, alternative = "two.sided", continuous.relax = F)
{
  require(gurobi)
  require(Matrix)
  
  
  
  ns = table(index)
  ms = table(index[treatment*1==1]) 
  N.total = sum(ns)
  
  
  nostratum = length(unique(index))
  treatment = 1*treatment
  treatment = (1*treatment == 1)
  gur = suppressWarnings(require(gurobi))
  
  
  vec.111 = tapply((treatment)*outcome*(dose), index, sum)
  vec.110 = tapply((treatment)*outcome*(!dose), index, sum)
  vec.010 = tapply((!treatment)*outcome*(!dose), index, sum)
  vec.101 = tapply((treatment)*(!outcome)*dose, index, sum)
  vec.100 = tapply((treatment)*(!outcome)*(!dose), index, sum)
  vec.001 = tapply((!treatment)*(!outcome)*dose, index, sum)
  vec.011 = tapply((!treatment)*(outcome)*dose, index, sum)
  vec.000 = tapply((!treatment)*(!outcome)*(!dose), index, sum)
  V = cbind(vec.111, vec.110, vec.101, vec.100, vec.011, vec.010,vec.001, vec.000)
  id = apply(V, 1, paste, collapse = "-")
  
  num.id =  sort.new(xtfrm(id))
  nosymm = length(unique(num.id))
  cc = table(num.id)
  bb = tapply(ns, num.id, mean)
  N.sr = sum(cc-1)
  m.011 = vec.011+1
  m.010 = vec.010+1
  m.001 = vec.001+1
  m.000 = vec.000+1
  m.101 = vec.101 + 1
  m.100 = vec.100 + 1
  m.111 = vec.111 + 1
  m.110 = vec.110 + 1
  
  
  mult.000 = tapply(m.000, num.id, mean)
  mult.001 = tapply(m.001, num.id, mean)
  mult.010 = tapply(m.010, num.id, mean)
  mult.011 = tapply(m.011, num.id, mean)
  mult.110 = tapply(m.110, num.id, mean)
  mult.111 = tapply(m.111, num.id, mean)
  mult.101 = tapply(m.101, num.id, mean)
  mult.100 = tapply(m.100, num.id, mean)
  
  wATE.per.strat = (ns)*(tapply((treatment)*outcome, index, sum)/ms - tapply((1-treatment)*outcome, index, sum)/(ns-ms))
  ATE.est = sum(wATE.per.strat)
  wATE2.per.strat = (ns)*(tapply((treatment)*dose, index, sum)/ms - tapply((1-treatment)*dose, index, sum)/(ns-ms))
  ATE2.est = sum(wATE2.per.strat)
  
  lambdahat = ATE.est/ATE2.est
  Gamma.vec = 1
  null = round(null, 2)
  
  wTE.per.strat = (ns)*(tapply((treatment)*(outcome - null*dose), index, sum)/ms - tapply((1-treatment)*(outcome - null*dose), index, sum)/(ns-ms))
  
  
  
  TE.est = sum(wTE.per.strat)
  max.e = (TE.est > 0)
  
  pvalvec = rep(0, length(Gamma.vec))
  Rejectvec = rep(0, length(Gamma.vec))
  SE = 0
  if(ER == F)
  {
    mult.po = (mult.111)*mult.110*(mult.101)*mult.100*mult.001*(mult.000)*(mult.010)*mult.011
    mult.do = (mult.111)*mult.110*(mult.101)*mult.100*mult.001*(mult.000)*(mult.010)*mult.011
    if(DE == "nonpositive")
    {
      mult.po = mult.100*mult.101*mult.010*mult.011
    }
    if(DE=="nonnegative")
    {
      mult.po = mult.111*mult.110*mult.000*mult.001
    }
    if(MO == T)
    {
      mult.do = mult.111*mult.101*mult.000*mult.010
    }
    mult.all = mult.po*mult.do
    ns.type = tapply(ns, num.id, mean)
    N.vars = sum(mult.all)
    index.symm = rep(1:nosymm, mult.all*(ns.type-1))
    n.types = (mult.all)
    n.po = mult.po
    n.do= mult.do
    n.per = cc
    Diff = rep(0, N.vars)
    Diff2 = rep(0, N.vars)
    #V.list = vector("list", N.vars)
    row.ind = rep(0, 4*N.vars+1)
    col.ind = row.ind
    values = row.ind
    b = rep(0, nosymm + 3)
    for(kk in 1:nosymm)
    {
      row.ind[which(index.symm==kk)] = rep(kk, (ns.type[kk]-1))
      col.ind[which(index.symm==kk)] = which(index.symm==kk)  
      values[which(index.symm==kk)] = rep(1, (ns.type[kk]-1))
      b[kk] = n.per[kk]
    }
    row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
    col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
    row.ind[(2*N.vars+1):(3*N.vars)] = rep(nosymm + 2, N.vars)
    col.ind[(2*N.vars+1):(3*N.vars)] = 1:N.vars
    row.ind[(3*N.vars+1):(4*N.vars+1)] = rep(nosymm + 3, N.vars+1)
    col.ind[(3*N.vars+1):(4*N.vars+1)] = 1:(N.vars+1)
    PM = rep(0, N.vars)
    PV = rep(0, N.vars)
    
    
    
    Gamma.vec = 1
    zscore = rep(0, length(Gamma.vec))
    
    
    for(kk in 1:nosymm)
    {
      i = which(num.id==kk)[1]
      symmgroup = which(index.symm == kk)
      ind = which(index==i)
      treatstrat = treatment[ind]
      outstrat = outcome[ind]
      dosestrat = dose[ind]
      dosesymm = c(rep(0, sum(treatstrat==F&dosestrat==F)), rep(1, sum(treatstrat==F& dosestrat==T)), rep(0, sum(treatstrat==T&dosestrat==F)), rep(1, sum(treatstrat==T& dosestrat==T)))
      
      outsymm = c(sort(outstrat[treatstrat==F& dosestrat==F]), sort(outstrat[treatstrat==F& dosestrat==T]), sort(outstrat[treatstrat==T&dosestrat == F]), sort(outstrat[treatstrat==T&dosestrat == T]))
      PO = matrix(0, n.po[kk], ns[i])
      DO = matrix(0, n.do[kk],ns[i])
      
      mt.011 = mult.011[kk]
      mt.001 = mult.001[kk]
      mt.000 = mult.000[kk]
      mt.010 = mult.010[kk]
      mt.110 = mult.110[kk]
      mt.100 = mult.100[kk]
      mt.111 = mult.111[kk]
      mt.101 = mult.101[kk]
      p.011 = mt.011
      p.001 = mt.001
      p.000 = mt.000
      p.010 = mt.010
      p.110 = mt.110
      p.100 = mt.100
      p.111 = mt.111
      p.101 = mt.101
      
      T.000 = matrix(1, mt.000-1, mt.000-1)
      T.000[lower.tri(T.000)] = 0
      T.000 = rbind(T.000, c(rep(0, mt.000-1)) )
      T.001 = matrix(1, mt.001-1, mt.001-1)
      T.001[lower.tri(T.001)] = 0
      T.001 = rbind(T.001, c(rep(0, mt.001-1)) )
      T.010 = matrix(1, mt.010-1, mt.010-1)
      T.010[lower.tri(T.010)] = 0
      T.010 = rbind(T.010, c(rep(0, mt.010-1)) )
      T.011 = matrix(1, mt.011-1, mt.011-1)
      T.011[lower.tri(T.011)] = 0
      T.011 = rbind(T.011, c(rep(0, mt.011-1)))
      T.100 = matrix(1, mt.100-1, mt.100-1)
      T.100[lower.tri(T.100)] = 0
      T.100 = rbind(T.100, c(rep(0, mt.100-1)) )
      T.101 = matrix(1, mt.101-1, mt.101-1)
      T.101[lower.tri(T.101)] = 0
      T.101 = rbind(T.101, c(rep(0, mt.101-1)) )
      T.110 = matrix(1, mt.110-1, mt.110-1)
      T.110[lower.tri(T.110)] = 0
      T.110 = rbind(T.110, c(rep(0, mt.110-1)) )
      T.111 = matrix(1, mt.111-1, mt.111-1)
      T.111[lower.tri(T.111)] = 0
      T.111 = rbind(T.111, c(rep(0, mt.111-1)))
      
      if(DE == "nonnegative")
      {
        T.100 = matrix(0, 1, mt.100-1)
        p.100 = 1
        T.101 = matrix(0, 1, mt.101-1)
        p.101 = 1
        T.011 = matrix(1, 1, mt.011-1)
        p.011 = 1
        T.010 = matrix(1, 1, mt.010-1)
        p.010 = 1
      }
      if(DE == "nonpositive")
      {
        T.111 = matrix(1, 1, mt.111-1)
        p.111 = 1
        T.110 = matrix(1, 1, mt.110-1)
        p.110 = 1
        T.001 = matrix(0, 1, mt.001-1)
        p.001 = 1
        T.000 = matrix(0, 1, mt.000-1)
        p.000 = 1
      }
      count = 1
      for(ll in 1:p.000)
      {
        for(mm in 1:p.010)
        {
          for(uu in 1:p.001)
          {
            for(vv in 1:p.011)
            {
              for(ww in 1:p.100)
              {
                for(xx in 1:p.110)
                {
                  for(yy in 1:p.101)
                  {
                    for(zz in 1:p.111)	
                    {
                      tempvec = c(T.000[ll,], T.010[mm,], T.001[uu,], T.011[vv,], T.100[ww,], T.110[xx,], T.101[yy,], T.111[zz,])
                      tempvec = tempvec[!is.na(tempvec)]
                      PO[count,] = tempvec
                      count = count+1
                    }
                  }
                }
              }
              
            }
            
          }
        }
        
      }
      count =1 
      p.011 = mt.011
      p.001 = mt.001
      p.000 = mt.000
      p.010 = mt.010
      p.110 = mt.110
      p.100 = mt.100
      p.111 = mt.111
      p.101 = mt.101
      
      T.000 = matrix(1, mt.000-1, mt.000-1)
      T.000[lower.tri(T.000)] = 0
      T.000 = rbind(T.000, c(rep(0, mt.000-1)) )
      T.001 = matrix(1, mt.001-1, mt.001-1)
      T.001[lower.tri(T.001)] = 0
      T.001 = rbind(T.001, c(rep(0, mt.001-1)) )
      T.010 = matrix(1, mt.010-1, mt.010-1)
      T.010[lower.tri(T.010)] = 0
      T.010 = rbind(T.010, c(rep(0, mt.010-1)) )
      T.011 = matrix(1, mt.011-1, mt.011-1)
      T.011[lower.tri(T.011)] = 0
      T.011 = rbind(T.011, c(rep(0, mt.011-1)))
      T.100 = matrix(1, mt.100-1, mt.100-1)
      T.100[lower.tri(T.100)] = 0
      T.100 = rbind(T.100, c(rep(0, mt.100-1)) )
      T.101 = matrix(1, mt.101-1, mt.101-1)
      T.101[lower.tri(T.101)] = 0
      T.101 = rbind(T.101, c(rep(0, mt.101-1)) )
      T.110 = matrix(1, mt.110-1, mt.110-1)
      T.110[lower.tri(T.110)] = 0
      T.110 = rbind(T.110, c(rep(0, mt.110-1)) )
      T.111 = matrix(1, mt.111-1, mt.111-1)
      T.111[lower.tri(T.111)] = 0
      T.111 = rbind(T.111, c(rep(0, mt.111-1)))
      
      if(MO == T)
      {
        T.100 = matrix(0, 1, mt.100-1)
        p.100 = 1
        T.110 = matrix(0, 1, mt.110-1)
        p.110 = 1
        T.011 = matrix(1, 1, mt.011-1)
        p.011 = 1
        T.001 = matrix(1, 1, mt.001-1)
        p.001 = 1
      }
      for(ll in 1:p.000)
      {
        for(mm in 1:p.010)
        {
          for(uu in 1:p.001)
          {
            for(vv in 1:p.011)
            {
              for(ww in 1:p.100)
              {
                for(xx in 1:p.110)
                {
                  for(yy in 1:p.101)
                  {
                    for(zz in 1:p.111)	
                    {
                      tempvec = c(T.000[ll,], T.010[mm,], T.001[uu,], T.011[vv,], T.100[ww,], T.110[xx,], T.101[yy,], T.111[zz,])
                      tempvec = tempvec[!is.na(tempvec)]
                      DO[count,] = tempvec
                      count = count+1
                    }
                  }
                }
              }
              
            }
            
          }
        }
        
      }
      
      
      
      
      
      count = 1
      
      for(jj in 1:n.po[kk])
      {
        for(ll in 1:n.do[kk])
        {
          
          po.symm = PO[jj,]
          do.symm = DO[ll,]
          treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
          outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
          outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm) 
          dosecontrol = dosesymm*(1-treatsymm) + do.symm*(treatsymm)
          dosetreat = dosesymm*(treatsymm) + do.symm*(1-treatsymm) 
          
          sum.cont = sum(outcontrol - null*dosecontrol)/(ns[i]-1)
          Q = (outtreat + outcontrol/(ns[i]-1) - null*dosetreat - null*dosecontrol/(ns[i]-1)- sum.cont)*ns[i]
          if(sum(treatstrat)>1)
          {
            sum.cont = sum(outtreat - null*dosetreat)/(ns[i]-1)
            Q = -(outtreat/(ns[i]-1) + outcontrol - null*dosetreat/(ns[i]-1) - null*dosecontrol- sum.cont)*ns[i]
            
          }
          
          
          
          qi = Q*max.e - Q*(!max.e)
          ord = order(qi)
          qi.sort = sort(qi)
          
          
          Gamma.sens = 1
          
          
          mu = mean(qi.sort)
          sigma2 = mean(qi.sort^2) - mu^2
          
          mu[abs(mu) < 1e-8] = 0
          sigma2[sigma2 < 1e-8] = 0
          
          
          PV[symmgroup[count]] = (sigma2)
          Diff[symmgroup[count]] = sum(outtreat -outcontrol - (null*(dosetreat - dosecontrol)))
          Diff2[symmgroup[count]] = sum(((dosetreat - dosecontrol)))
          count = count+1
        }
      }
    }
    
    
    values[(N.vars+1):(2*N.vars)] = Diff
    values[(2*N.vars+1):(3*N.vars)] = Diff2
    values[(3*N.vars+1):(4*N.vars+1)] = c(-PM, 1)
    b[nosymm+1] = 0
    b[nosymm+2] = 1
    b[nosymm+3] = 0
    
    
    const.dir = c(rep("=", nosymm+1), ">=", "=")
    model = list()
    
    model$A = sparseMatrix(row.ind[1:(3*N.vars)], col.ind[1:(3*N.vars)], x=values[1:(3*N.vars)])
    model$obj = c(PV)
    model$sense = const.dir[1:(nosymm+2)]
    model$rhs = b[1:(nosymm+2)]
    model$vtype = c(rep("I", N.vars))
    if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}
    
    
    model$modelsense = "max"
    
    
    solm = gurobi(model, params = list(OutputFlag = 0))
    zed = (TE.est/sqrt(solm$objval))
    SE = sqrt(solm$objval)
    x = solm$x[1:N.vars]
    tstat = zed
    pval = 0
    if(alternative == "two.sided")
    {
      pval = 2*pnorm(-abs(tstat))
    }
    if(alternative == "greater")
    {
      pval = 1 - pnorm((tstat))
    }
    if(alternative == "less")
    {
      pval = pnorm((tstat))
    }
    
    
    
  }
  if(ER == T)
  {
    mult.all = (mult.101*(mult.101+1)/2)*(mult.111*(mult.111+1)/2)*(mult.010*(mult.010+1)/2)*(mult.000*(mult.000+1)/2)*(mult.100*(mult.100+1)/2)*(mult.110*(mult.110+1)/2)*(mult.011*(mult.011+1)/2)*(mult.001*(mult.001+1)/2)
    if(MO == T)
    {
      mult.all = (mult.101*(mult.101+1)/2)*(mult.111*(mult.111+1)/2)*(mult.010*(mult.010+1)/2)*(mult.000*(mult.000+1)/2)
      if(DE == "nonpositive")
      {
        mult.all = (mult.111)*(mult.101*(mult.101+1)/2)*(mult.000)*(mult.010*(mult.010+1)/2)
      }
      if(DE=="nonnegative")
      {
        mult.all = (mult.101)*(mult.111*(mult.111+1)/2)*(mult.010)*(mult.000*(mult.000+1)/2)
      }
    }
    if(MO == F)
    {
      mult.all = (mult.101*(mult.101+1)/2)*(mult.111*(mult.111+1)/2)*(mult.010*(mult.010+1)/2)*(mult.000*(mult.000+1)/2)*(mult.100*(mult.100+1)/2)*(mult.110*(mult.110+1)/2)*(mult.011*(mult.011+1)/2)*(mult.001*(mult.001+1)/2)
      if(DE == "nonpositive")
      {
        mult.all = (mult.111)*(mult.101*(mult.101+1)/2)*(mult.000)*(mult.010*(mult.010+1)/2)*mult.001*(mult.011*(mult.011+1)/2)*mult.110*(mult.100*(mult.100+1)/2)
      }
      if(DE=="nonnegative")
      {
        mult.all = (mult.101)*(mult.111*(mult.111+1)/2)*(mult.010)*(mult.000*(mult.000+1)/2)*mult.011*(mult.001*(mult.001+1)/2)*mult.100*(mult.110*(mult.110+1)/2)
      }
    }
    
    ns.type = tapply(ns, num.id, mean)
    N.vars = sum(mult.all)
    index.symm = rep(1:nosymm, mult.all)
    n.types = (mult.all)
    mult.po= mult.all
    mult.do = mult.all
    n.po = mult.po
    n.do= mult.do
    n.per = cc
    Diff = rep(0, N.vars)
    Diff2 = rep(0, N.vars)
    #V.list = vector("list", N.vars)
    row.ind = rep(0, 4*N.vars+1)
    col.ind = row.ind
    values = row.ind
    b = rep(0, nosymm + 3)
    for(kk in 1:nosymm)
    {
      row.ind[which(index.symm==kk)] = rep(kk, n.types[kk])
      col.ind[which(index.symm==kk)] = which(index.symm==kk)  
      values[which(index.symm==kk)] = rep(1, (n.types[kk]))
      b[kk] = n.per[kk]
    }
    row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
    col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
    row.ind[(2*N.vars+1):(3*N.vars)] = rep(nosymm + 2, N.vars)
    col.ind[(2*N.vars+1):(3*N.vars)] = 1:N.vars
    row.ind[(3*N.vars+1):(4*N.vars+1)] = rep(nosymm + 3, N.vars+1)
    col.ind[(3*N.vars+1):(4*N.vars+1)] = 1:(N.vars+1)
    PM = rep(0, N.vars)
    PV = rep(0, N.vars)
    
    
    Gamma.vec = 1
    Gamma.sens = 1
    
    for(kk in 1:nosymm)
    {
      i = which(num.id==kk)[1]
      symmgroup = which(index.symm == kk)
      ind = which(index==i)
      treatstrat = treatment[ind]
      outstrat = outcome[ind]
      dosestrat = dose[ind]
      dosesymm = c(rep(0, sum(treatstrat==F&dosestrat==F)), rep(1, sum(treatstrat==F& dosestrat==T)), rep(0, sum(treatstrat==T&dosestrat==F)), rep(1, sum(treatstrat==T& dosestrat==T)))
      
      outsymm = c(sort(outstrat[treatstrat==F& dosestrat==F]), sort(outstrat[treatstrat==F& dosestrat==T]), sort(outstrat[treatstrat==T&dosestrat == F]), sort(outstrat[treatstrat==T&dosestrat == T]))
      PO = matrix(0, n.po[kk], ns[i])
      DO = matrix(0, n.po[kk],ns[i])
      
      mt.011 = mult.011[kk]
      mt.001 = mult.001[kk]
      mt.000 = mult.000[kk]
      mt.010 = mult.010[kk]
      mt.110 = mult.110[kk]
      mt.100 = mult.100[kk]
      mt.111 = mult.111[kk]
      mt.101 = mult.101[kk]
      p.011 = mt.011
      p.001 = mt.001
      p.000 = mt.000
      p.010 = mt.010
      p.110 = mt.110
      p.100 = mt.100
      p.111 = mt.111
      p.101 = mt.101
      
      D.000 = matrix(1, mt.000-1, mt.000-1)
      D.000[lower.tri(D.000)] = 0
      D.000 = rbind(D.000, c(rep(0, mt.000-1)) )
      D.010 = matrix(1, mt.010-1, mt.010-1)
      D.010[lower.tri(D.010)] = 0
      D.010 = rbind(D.010, c(rep(0, mt.010-1)) )
      D.001 = matrix(1, mt.001-1, mt.001-1)
      D.001[lower.tri(D.001)] = 0
      D.001 = rbind(D.001, c(rep(0, mt.001-1)) )
      D.011 = matrix(1, mt.011-1, mt.011-1)
      D.011[lower.tri(D.011)] = 0
      D.011 = rbind(D.011, c(rep(0, mt.011-1)))
      D.100 = matrix(1, mt.100-1, mt.100-1)
      D.100[lower.tri(D.100)] = 0
      D.100 = rbind(D.100, c(rep(0, mt.100-1)) )
      D.110 = matrix(1, mt.110-1, mt.110-1)
      D.110[lower.tri(D.110)] = 0
      D.110 = rbind(D.110, c(rep(0, mt.110-1)) )
      D.101 = matrix(1, mt.101-1, mt.101-1)
      D.101[lower.tri(D.101)] = 0
      D.101 = rbind(D.101, c(rep(0, mt.101-1)) )
      
      D.111 = matrix(1, mt.111-1, mt.111-1)
      D.111[lower.tri(D.111)] = 0
      D.111 = rbind(D.111, c(rep(0, mt.111-1)))
      
      if(MO == T)
      {
        D.100 = matrix(0, 1, mt.100-1)
        p.100 = 1
        D.110 = matrix(0, 1, mt.110-1)
        p.110 = 1
        D.011 = matrix(1, 1, mt.011-1)
        p.011 = 1
        D.001 = matrix(1, 1, mt.001-1)
        p.001 = 1
      }
      count = 1
      for(ll in 1:p.000)
      {
        for(mm in 1:p.010)
        {
          for(uu in 1:p.001)
          {
            for(vv in 1:p.011)
            {
              for(ww in 1:p.100)
              {
                for(xx in 1:p.110)
                {
                  for(yy in 1:p.101)
                  {
                    for(zz in 1:p.111)	
                    {
                      T.000 = matrix(1, p.000-ll, mt.000-ll)
                      T.000[lower.tri(T.000)] = 0
                      T.000 = rbind(T.000, c(rep(0, mt.000-ll)))
                      F.000 = matrix(0, (p.000-ll+1), (ll-1))
                      T.000 = cbind(F.000, T.000)
                      
                      T.010 = matrix(1, p.010-mm, mt.010-mm)
                      T.010[lower.tri(T.010)] = 0
                      T.010 = rbind(T.010, c(rep(0, mt.010-mm)))
                      F.010 = matrix(1, (p.010-mm+1), (mm-1))
                      T.010 = cbind(F.010, T.010)
                      
                      T.001 = matrix(1, uu-1, uu-1)
                      T.001[lower.tri(T.001)] = 0
                      T.001 = rbind(T.001, c(rep(0, uu-1)))
                      F.001 = matrix(0, (uu), (mt.001-uu))
                      T.001 = cbind(F.001, T.001)
                      if(MO == T)
                      {
                        T.001 = matrix(0, 1, mt.001-1)
                      }
                      
                      T.011 = matrix(1, vv-1, vv-1)
                      T.011[lower.tri(T.011)] = 0
                      T.011 = rbind(T.011, c(rep(0, vv-1)))
                      F.011 = matrix(1, (vv), (mt.011-vv))
                      T.011 = cbind(F.011, T.011)
                      if(MO == T)
                      {
                        T.011 = matrix(1, 1, mt.011-1)
                      }
                      
                      
                      T.100 = matrix(1, p.100-ww, mt.100-ww)
                      T.100[lower.tri(T.100)] = 0
                      T.100 = rbind(T.100, c(rep(0, mt.100-ww)))
                      F.100 = matrix(0, (p.100-ww+1), (ww-1))
                      T.100 = cbind(F.100, T.100)
                      if(MO == T)
                      {
                        T.100 = matrix(0, 1,mt.100-1)
                      }
                      
                      T.110 = matrix(1, p.110-xx, mt.110-xx)
                      T.110[lower.tri(T.110)] = 0
                      T.110 = rbind(T.110, c(rep(0, mt.110-xx)))
                      F.110 = matrix(1, (p.110-xx+1), (xx-1))
                      T.110 = cbind(F.110, T.110)
                      if(MO == T)
                      {
                        T.110 = matrix(1, 1, mt.110-1)
                      }
                      
                      T.101 = matrix(1, yy-1, yy-1)
                      T.101[lower.tri(T.101)] = 0
                      T.101 = rbind(T.101, c(rep(0, yy-1)))
                      F.101 = matrix(0, (yy), (mt.101-yy))
                      T.101 = cbind(F.101, T.101)
                      
                      T.111 = matrix(1, zz-1, zz-1)
                      T.111[lower.tri(T.111)] = 0
                      T.111 = rbind(T.111, c(rep(0, zz-1)))
                      F.111 = matrix(1, (zz), (mt.111-zz))
                      T.111 = cbind(F.111, T.111) 
                      
                      
                      
                      
                      
                      if(DE == "nonnegative")
                      {
                        T.100 = matrix(0, 1, mt.100-1)
                        T.101 = matrix(0, 1, mt.101-1)
                        T.011 = matrix(1, 1, mt.011-1)
                        T.010 = matrix(1, 1, mt.010-1)
                      }
                      if(DE == "nonpositive")
                      {
                        T.111 = matrix(1, 1, mt.111-1)
                        T.110 = matrix(1, 1, mt.110-1)
                        T.001 = matrix(0, 1, mt.001-1)
                        T.000 = matrix(0, 1, mt.000-1)
                        
                      }
                      
                      
                      
                      dosevec = c(D.000[ll,], D.010[mm,], D.001[uu,], D.011[vv,], D.100[ww,], D.110[xx,], D.101[yy,], D.111[zz,])
                      dosevec = dosevec[!is.na(dosevec)]
                      for(lll in 1:nrow(T.000))
                      {
                        for(mmm in 1:nrow(T.010))
                        {
                          for(uuu in 1:nrow(T.001))
                          {
                            for(vvv in 1:nrow(T.011))
                            {
                              for(www in 1:nrow(T.100))
                              {
                                for(xxx in 1:nrow(T.110))
                                {
                                  for(yyy in 1:nrow(T.101))
                                  {
                                    for(zzz in 1:nrow(T.111))	
                                    {
                                      DO[count,] = dosevec
                                      tempvec = c(T.000[lll,], T.010[mmm,], T.001[uuu,], T.011[vvv,], T.100[www,], T.110[xxx,], T.101[yyy,], T.111[zzz,])
                                      tempvec = tempvec[!is.na(tempvec)]
                                      PO[count,] = tempvec
                                      count = count+1
                                    }}}}}}}}}}}}}}}}             
      
      
      
      count = 1
      
      for(jj in 1:n.po[kk])
      {
        
        po.symm = PO[jj,]
        do.symm = DO[jj,]
        treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
        outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
        outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm) 
        dosecontrol = dosesymm*(1-treatsymm) + do.symm*(treatsymm)
        dosetreat = dosesymm*(treatsymm) + do.symm*(1-treatsymm) 
        sum.cont = sum(outcontrol - null*dosecontrol)/(ns[i]-1)
        Q = (outtreat + outcontrol/(ns[i]-1) - null*dosetreat - null*dosecontrol/(ns[i]-1)- sum.cont)*ns[i]
        if(sum(treatstrat)>1)
        {
          sum.cont = sum(outtreat - null*dosetreat)/(ns[i]-1)
          Q = -(outtreat/(ns[i]-1) + outcontrol - null*dosetreat/(ns[i]-1) - null*dosecontrol- sum.cont)*ns[i]
          
        }
        
        
        
        qi = Q*max.e - Q*(!max.e)
        ord = order(qi)
        qi.sort = sort(qi)
        
        
        
        
        mu = mean(qi.sort)        
        sigma2 = mean(qi.sort^2) - mu^2
        
        mu[abs(mu) < 1e-8] = 0
        sigma2[sigma2 < 1e-8] = 0
        
        
        
        PV[symmgroup[jj]] = (sigma2)
        
        Diff[symmgroup[jj]] = sum(outtreat -outcontrol - (null*(dosetreat - dosecontrol)))
        Diff2[symmgroup[jj]] = sum(dosetreat-dosecontrol)     
        
      }
    }
    
    
    values[(N.vars+1):(2*N.vars)] = Diff
    values[(2*N.vars+1):(3*N.vars)] = Diff2
    values[(3*N.vars+1):(4*N.vars+1)] = c(-PM, 1)
    b[nosymm+1] = 0
    b[nosymm+2] = 1
    b[nosymm+3] = 0
    
    
    const.dir = c(rep("=", nosymm+1), ">=", "=")
    model = list()
    
    model$A = sparseMatrix(row.ind[1:(3*N.vars)], col.ind[1:(3*N.vars)], x=values[1:(3*N.vars)])
    model$obj = c(PV)
    model$sense = const.dir[1:(nosymm+2)]
    model$rhs = b[1:(nosymm+2)]
    model$vtype = c(rep("I", N.vars))
    if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}
    
    
    model$modelsense = "max"
    
    
    solm = gurobi(model, params = list(OutputFlag = 0))
    zed = (TE.est/sqrt(solm$objval))
    x = solm$x[1:N.vars]
    SE = sqrt(solm$objval)
    tstat = zed
    pval = 0
    if(alternative == "two.sided")
    {
      pval = 2*pnorm(-abs(tstat))
    }
    if(alternative == "greater")
    {
      pval = 1 - pnorm((tstat))
    }
    if(alternative == "less")
    {
      pval = pnorm((tstat))
    }
    
  }
  
  
  
  return(list(SE = SE, pval = pval, tstat = tstat, lambdahat = round(lambdahat,2)))
  
}


###########
#sensEffectRatio
###########
sensEffectRatio = function(index, treatment, outcome, dose, null=0, DE = "both", MO = T, ER = T, alternative = "two.sided", alpha = 0.05, Gamma.vec = 1, calculate.pval = T, continuous.relax = F)
{
  PVAL = calculate.pval
  require(gurobi)
  require(Matrix)
  
  
  
  ns = table(index)
  ms = table(index[treatment*1==1]) 
  N.total = sum(ns)
  
  
  nostratum = length(unique(index))
  treatment = 1*treatment
  treatment = (1*treatment == 1)
  gur = suppressWarnings(require(gurobi))
  
  
  vec.111 = tapply((treatment)*outcome*(dose), index, sum)
  vec.110 = tapply((treatment)*outcome*(!dose), index, sum)
  vec.010 = tapply((!treatment)*outcome*(!dose), index, sum)
  vec.101 = tapply((treatment)*(!outcome)*dose, index, sum)
  vec.100 = tapply((treatment)*(!outcome)*(!dose), index, sum)
  vec.001 = tapply((!treatment)*(!outcome)*dose, index, sum)
  vec.011 = tapply((!treatment)*(outcome)*dose, index, sum)
  vec.000 = tapply((!treatment)*(!outcome)*(!dose), index, sum)
  V = cbind(vec.111, vec.110, vec.101, vec.100, vec.011, vec.010,vec.001, vec.000)
  id = apply(V, 1, paste, collapse = "-")
  
  num.id =  sort.new(xtfrm(id))
  nosymm = length(unique(num.id))
  cc = table(num.id)
  bb = tapply(ns, num.id, mean)
  N.sr = sum(cc-1)
  m.011 = vec.011+1
  m.010 = vec.010+1
  m.001 = vec.001+1
  m.000 = vec.000+1
  m.101 = vec.101 + 1
  m.100 = vec.100 + 1
  m.111 = vec.111 + 1
  m.110 = vec.110 + 1
  
  
  mult.000 = tapply(m.000, num.id, mean)
  mult.001 = tapply(m.001, num.id, mean)
  mult.010 = tapply(m.010, num.id, mean)
  mult.011 = tapply(m.011, num.id, mean)
  mult.110 = tapply(m.110, num.id, mean)
  mult.111 = tapply(m.111, num.id, mean)
  mult.101 = tapply(m.101, num.id, mean)
  mult.100 = tapply(m.100, num.id, mean)
  
  wATE.per.strat = (ns)*(tapply((treatment)*outcome, index, sum)/ms - tapply((1-treatment)*outcome, index, sum)/(ns-ms))
  ATE.est = sum(wATE.per.strat)
  wATE2.per.strat = (ns)*(tapply((treatment)*dose, index, sum)/ms - tapply((1-treatment)*dose, index, sum)/(ns-ms))
  ATE2.est = sum(wATE2.per.strat)
  
  lambdahat = ATE.est/ATE2.est
  #Gamma.vec = round(Gamma.vec, 2)
  null = round(null, 2)
  
  
  wTE.per.strat = (ns)*(tapply((treatment)*(outcome - null*dose), index, sum)/ms - tapply((1-treatment)*(outcome - null*dose), index, sum)/(ns-ms))
  
  
  
  TE.est = sum(wTE.per.strat)
  max.e = (TE.est > 0)
  
  pvalvec = rep(0, length(Gamma.vec))
  Rejectvec = rep(0, length(Gamma.vec))
  kappavec = rep(0, length(Gamma.vec))
  SE = 0
  if(ER == F)
  {
    mult.po = (mult.111)*mult.110*(mult.101)*mult.100*mult.001*(mult.000)*(mult.010)*mult.011
    mult.do = (mult.111)*mult.110*(mult.101)*mult.100*mult.001*(mult.000)*(mult.010)*mult.011
    if(DE == "nonpositive")
    {
      mult.po = mult.100*mult.101*mult.010*mult.011
    }
    if(DE=="nonnegative")
    {
      mult.po = mult.111*mult.110*mult.000*mult.001
    }
    if(MO == T)
    {
      mult.do = mult.111*mult.101*mult.000*mult.010
    }
    mult.all = mult.po*mult.do
    ns.type = tapply(ns, num.id, mean)
    N.vars = sum(mult.all*(ns.type-1))
    index.symm = rep(1:nosymm, mult.all*(ns.type-1))
    n.types = (mult.all)
    n.po = mult.po
    n.do= mult.do
    n.per = cc
    Diff = rep(0, N.vars)
    Diff2 = Diff
    #V.list = vector("list", N.vars)
    row.ind = rep(0, 4*N.vars+1)
    col.ind = row.ind
    values = row.ind
    b = rep(0, nosymm + 2)
    for(kk in 1:nosymm)
    {
      row.ind[which(index.symm==kk)] = rep(kk, (ns.type[kk]-1)*n.types[kk])
      col.ind[which(index.symm==kk)] = which(index.symm==kk)  
      values[which(index.symm==kk)] = rep(1, (ns.type[kk]-1)*(n.types[kk]))
      b[kk] = n.per[kk]
    }
    row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
    col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
    row.ind[(2*N.vars+1):(3*N.vars)] = rep(nosymm + 2, N.vars)
    col.ind[(2*N.vars+1):(3*N.vars)] = 1:N.vars
    row.ind[(3*N.vars+1):(4*N.vars+1)] = rep(nosymm + 3, N.vars+1)
    col.ind[(3*N.vars+1):(4*N.vars+1)] = 1:(N.vars+1)
    PM = rep(0, N.vars)
    PV = rep(0, N.vars)
    
    
    
    #Gamma.vec = seq(1.2, 1.22, by = .001)
    
    zscore = rep(0, length(Gamma.vec))
   
    for(ee in 1:length(Gamma.vec))
    {
      Gamma.sens = Gamma.vec[ee]
      
      for(kk in 1:nosymm)
      {
        i = which(num.id==kk)[1]
        symmgroup = which(index.symm == kk)
        ind = which(index==i)
        treatstrat = treatment[ind]
        outstrat = outcome[ind]
        dosestrat = dose[ind]
        dosesymm = c(rep(0, sum(treatstrat==F&dosestrat==F)), rep(1, sum(treatstrat==F& dosestrat==T)), rep(0, sum(treatstrat==T&dosestrat==F)), rep(1, sum(treatstrat==T& dosestrat==T)))
        
        outsymm = c(sort(outstrat[treatstrat==F& dosestrat==F]), sort(outstrat[treatstrat==F& dosestrat==T]), sort(outstrat[treatstrat==T&dosestrat == F]), sort(outstrat[treatstrat==T&dosestrat == T]))
        PO = matrix(0, n.po[kk], ns[i])
        DO = matrix(0, n.do[kk],ns[i])
        
        mt.011 = mult.011[kk]
        mt.001 = mult.001[kk]
        mt.000 = mult.000[kk]
        mt.010 = mult.010[kk]
        mt.110 = mult.110[kk]
        mt.100 = mult.100[kk]
        mt.111 = mult.111[kk]
        mt.101 = mult.101[kk]
        p.011 = mt.011
        p.001 = mt.001
        p.000 = mt.000
        p.010 = mt.010
        p.110 = mt.110
        p.100 = mt.100
        p.111 = mt.111
        p.101 = mt.101
        
        T.000 = matrix(1, mt.000-1, mt.000-1)
        T.000[lower.tri(T.000)] = 0
        T.000 = rbind(T.000, c(rep(0, mt.000-1)) )
        T.001 = matrix(1, mt.001-1, mt.001-1)
        T.001[lower.tri(T.001)] = 0
        T.001 = rbind(T.001, c(rep(0, mt.001-1)) )
        T.010 = matrix(1, mt.010-1, mt.010-1)
        T.010[lower.tri(T.010)] = 0
        T.010 = rbind(T.010, c(rep(0, mt.010-1)) )
        T.011 = matrix(1, mt.011-1, mt.011-1)
        T.011[lower.tri(T.011)] = 0
        T.011 = rbind(T.011, c(rep(0, mt.011-1)))
        T.100 = matrix(1, mt.100-1, mt.100-1)
        T.100[lower.tri(T.100)] = 0
        T.100 = rbind(T.100, c(rep(0, mt.100-1)) )
        T.101 = matrix(1, mt.101-1, mt.101-1)
        T.101[lower.tri(T.101)] = 0
        T.101 = rbind(T.101, c(rep(0, mt.101-1)) )
        T.110 = matrix(1, mt.110-1, mt.110-1)
        T.110[lower.tri(T.110)] = 0
        T.110 = rbind(T.110, c(rep(0, mt.110-1)) )
        T.111 = matrix(1, mt.111-1, mt.111-1)
        T.111[lower.tri(T.111)] = 0
        T.111 = rbind(T.111, c(rep(0, mt.111-1)))
        
        if(DE == "nonnegative")
        {
          T.100 = matrix(0, 1, mt.100-1)
          p.100 = 1
          T.101 = matrix(0, 1, mt.101-1)
          p.101 = 1
          T.011 = matrix(1, 1, mt.011-1)
          p.011 = 1
          T.010 = matrix(1, 1, mt.010-1)
          p.010 = 1
        }
        if(DE == "nonpositive")
        {
          T.111 = matrix(1, 1, mt.111-1)
          p.111 = 1
          T.110 = matrix(1, 1, mt.110-1)
          p.110 = 1
          T.001 = matrix(0, 1, mt.001-1)
          p.001 = 1
          T.000 = matrix(0, 1, mt.000-1)
          p.000 = 1
        }
        count = 1
        for(ll in 1:p.000)
        {
          for(mm in 1:p.010)
          {
            for(uu in 1:p.001)
            {
              for(vv in 1:p.011)
              {
                for(ww in 1:p.100)
                {
                  for(xx in 1:p.110)
                  {
                    for(yy in 1:p.101)
                    {
                      for(zz in 1:p.111)	
                      {
                        tempvec = c(T.000[ll,], T.010[mm,], T.001[uu,], T.011[vv,], T.100[ww,], T.110[xx,], T.101[yy,], T.111[zz,])
                        tempvec = tempvec[!is.na(tempvec)]
                        PO[count,] = tempvec
                        count = count+1
                      }
                    }
                  }
                }
                
              }
              
            }
          }
          
        }
        count =1 
        p.011 = mt.011
        p.001 = mt.001
        p.000 = mt.000
        p.010 = mt.010
        p.110 = mt.110
        p.100 = mt.100
        p.111 = mt.111
        p.101 = mt.101
        
        T.000 = matrix(1, mt.000-1, mt.000-1)
        T.000[lower.tri(T.000)] = 0
        T.000 = rbind(T.000, c(rep(0, mt.000-1)) )
        T.001 = matrix(1, mt.001-1, mt.001-1)
        T.001[lower.tri(T.001)] = 0
        T.001 = rbind(T.001, c(rep(0, mt.001-1)) )
        T.010 = matrix(1, mt.010-1, mt.010-1)
        T.010[lower.tri(T.010)] = 0
        T.010 = rbind(T.010, c(rep(0, mt.010-1)) )
        T.011 = matrix(1, mt.011-1, mt.011-1)
        T.011[lower.tri(T.011)] = 0
        T.011 = rbind(T.011, c(rep(0, mt.011-1)))
        T.100 = matrix(1, mt.100-1, mt.100-1)
        T.100[lower.tri(T.100)] = 0
        T.100 = rbind(T.100, c(rep(0, mt.100-1)) )
        T.101 = matrix(1, mt.101-1, mt.101-1)
        T.101[lower.tri(T.101)] = 0
        T.101 = rbind(T.101, c(rep(0, mt.101-1)) )
        T.110 = matrix(1, mt.110-1, mt.110-1)
        T.110[lower.tri(T.110)] = 0
        T.110 = rbind(T.110, c(rep(0, mt.110-1)) )
        T.111 = matrix(1, mt.111-1, mt.111-1)
        T.111[lower.tri(T.111)] = 0
        T.111 = rbind(T.111, c(rep(0, mt.111-1)))
        
        if(MO == T)
        {
          T.100 = matrix(0, 1, mt.100-1)
          p.100 = 1
          T.110 = matrix(0, 1, mt.110-1)
          p.110 = 1
          T.011 = matrix(1, 1, mt.011-1)
          p.011 = 1
          T.001 = matrix(1, 1, mt.001-1)
          p.001 = 1
        }
        for(ll in 1:p.000)
        {
          for(mm in 1:p.010)
          {
            for(uu in 1:p.001)
            {
              for(vv in 1:p.011)
              {
                for(ww in 1:p.100)
                {
                  for(xx in 1:p.110)
                  {
                    for(yy in 1:p.101)
                    {
                      for(zz in 1:p.111)	
                      {
                        tempvec = c(T.000[ll,], T.010[mm,], T.001[uu,], T.011[vv,], T.100[ww,], T.110[xx,], T.101[yy,], T.111[zz,])
                        tempvec = tempvec[!is.na(tempvec)]
                        DO[count,] = tempvec
                        count = count+1
                      }
                    }
                  }
                }
                
              }
              
            }
          }
          
        }
        
        
        
        
        
        count = 1
        
        for(jj in 1:n.po[kk])
        {
          for(ll in 1:n.do[kk])
          {
            ind.jj = (((count-1)*(ns[i]-1))+1):(count*(ns[i]-1))
            po.symm = PO[jj,]
            do.symm = DO[ll,]
            treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
            outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
            outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm) 
            dosecontrol = dosesymm*(1-treatsymm) + do.symm*(treatsymm)
            dosetreat = dosesymm*(treatsymm) + do.symm*(1-treatsymm) 
            
            sum.cont = sum(outcontrol - null*dosecontrol)/(ns[i]-1)
            Q = (outtreat + outcontrol/(ns[i]-1) - null*dosetreat - null*dosecontrol/(ns[i]-1)- sum.cont)*ns[i]
            if(sum(treatstrat)>1)
            {
              sum.cont = sum(outtreat - null*dosetreat)/(ns[i]-1)
              Q = -(outtreat/(ns[i]-1) + outcontrol - null*dosetreat/(ns[i]-1) - null*dosecontrol- sum.cont)*ns[i]
              
            }
            
            
            
            qi = Q*max.e - Q*(!max.e)
            ord = order(qi)
            qi.sort = sort(qi)
            
            
            mu = rep(0, length(ind)-1)
            sigma2 = rep(0, length(ind)-1)
            
            
            for(j in 1:(length(ind)-1))
            {
              mu[j] = (sum(qi.sort[1:(j)]) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]))/((j) + Gamma.sens*(ns[i]-(j)))
              sigma2[j] = (sum(qi.sort[1:(j)]^2) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]^2))/((j) + Gamma.sens*(ns[i]-(j))) - mu[j]^2
            }
            mu[abs(mu) < 1e-8] = 0
            sigma2[sigma2 < 1e-8] = 0
            
            
            PM[symmgroup[ind.jj]] = mu*(max.e) - mu*(!max.e)
            PV[symmgroup[ind.jj]] = (sigma2)
            
            Diff[symmgroup[ind.jj]] = sum(outtreat -outcontrol - (null*(dosetreat - dosecontrol)))
            Diff2[symmgroup[count]] = sum(((dosetreat - dosecontrol)))
            count = count+1
          }
        }
      }
      
      
      values[(N.vars+1):(2*N.vars)] = Diff
      values[(2*N.vars+1):(3*N.vars)] = Diff2
      values[(3*N.vars+1):(4*N.vars+1)] = c(-PM, 1)
      b[nosymm+1] = 0
      b[nosymm+2] = 1
      b[nosymm+3] = 0
      alpha.opt = alpha
      if(alternative != "two.sided")
      {
        alpha.opt = 2*alpha
      }
      
      const.dir = c(rep("=", nosymm+1), ">=", "=")
      model = list()
      if(Gamma.sens==1)
      {
        model$A = sparseMatrix(row.ind[1:(3*N.vars)], col.ind[1:(3*N.vars)], x=values[1:(3*N.vars)])
        model$obj = c(PV)
        model$sense = const.dir[1:(nosymm+2)]
        model$rhs = b[1:(nosymm+2)]
        model$vtype = c(rep("I", N.vars))
        if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}
        
        
        model$modelsense = "max"
        
        
        solm = gurobi(model, params = list(OutputFlag = 0))
        zed = (TE.est/sqrt(solm$objval))
        kappavec[ee] = (TE.est)^2 - qchisq(1-alpha.opt, 1)*solm$objval
        SE = sqrt(solm$objval)
        x = solm$x[1:N.vars]
        tstat = zed
        pval = 0
        if(alternative == "two.sided")
        {
          pval = 2*pnorm(-abs(tstat))
        }
        if(alternative == "greater")
        {
          pval = 1 - pnorm((tstat))
        }
        if(alternative == "less")
        {
          pval = pnorm((tstat))
        }
        Reject = (pval < alpha)
      }
      if(Gamma.sens != 1)
      {
        diff = 200
        kappa = qchisq(1-alpha.opt, 1)
        while(diff > 1e-8)
        {
          Plin = -2*TE.est*PM - kappa*PV 
          rowind.q =  1:(N.vars+1)
          colind.q = 1:(N.vars+1)
          values.q = c(rep(0, N.vars),1)
          Q = sparseMatrix(rowind.q, colind.q, x=values.q)
          model$A = sparseMatrix(row.ind, col.ind, x=values)
          model$obj = c(Plin,0)
          model$Q = Q
          model$sense = const.dir
          model$rhs = b
          
          
          model$vtype = c(rep("I", N.vars), "C")
          if(continuous.relax == T){model$vtype = c(rep("C", N.vars), "C")}
          model$lb = c(rep(0, N.vars), -Inf)
          
          
          model$modelsense = "min"
          
          
          solm = gurobi(model, params = list(OutputFlag = 0))
          x = solm$x[1:N.vars]
          kappa.new = (TE.est - sum(PM*x))^2/sum(PV*x)
          diff = abs(kappa.new - kappa)
          pval = 0
          if(PVAL == F)
          {
            diff = 0
            Reject = (kappa.new > kappa)
          }
          kappavec[ee] = (TE.est - sum(PM*x))^2 - qchisq(1-alpha.opt, 1)*sum(PV*x)
          kappa = kappa.new
        }
        zed = sqrt((TE.est - sum(PM*x))^2/sum(PV*x))
        
        if(alternative == "less")
        {
          zed = -zed
        }
        zscore[ee] = zed
        tstat = zed
        
        
        if(alternative == "two.sided")
        {
          pval = 2*pnorm(-abs(tstat))
        }
        if(alternative == "greater")
        {
          pval = 1 - pnorm((tstat))
        }
        if(alternative == "less")
        {
          pval = pnorm((tstat))
        }
        Reject = (pval < alpha)
        
        
        
        if(sign(TE.est - sum(PM*x))!=sign(TE.est))
        {
          Reject = F
          pval = 0.5
          if(alternative == "two.sided")
          {
            pval = 1
          }
        }
        
      }
      pvalvec[ee] = pval
      Rejectvec[ee] = Reject  
    }
  }
  if(ER == T)
  {
    mult.all = (mult.101*(mult.101+1)/2)*(mult.111*(mult.111+1)/2)*(mult.010*(mult.010+1)/2)*(mult.000*(mult.000+1)/2)*(mult.100*(mult.100+1)/2)*(mult.110*(mult.110+1)/2)*(mult.011*(mult.011+1)/2)*(mult.001*(mult.001+1)/2)
    if(MO == T)
    {
      mult.all = (mult.101*(mult.101+1)/2)*(mult.111*(mult.111+1)/2)*(mult.010*(mult.010+1)/2)*(mult.000*(mult.000+1)/2)
      if(DE == "nonpositive")
      {
        mult.all = (mult.111)*(mult.101*(mult.101+1)/2)*(mult.000)*(mult.010*(mult.010+1)/2)
      }
      if(DE=="nonnegative")
      {
        mult.all = (mult.101)*(mult.111*(mult.111+1)/2)*(mult.010)*(mult.000*(mult.000+1)/2)
      }
    }
    if(MO == F)
    {
      mult.all = (mult.101*(mult.101+1)/2)*(mult.111*(mult.111+1)/2)*(mult.010*(mult.010+1)/2)*(mult.000*(mult.000+1)/2)*(mult.100*(mult.100+1)/2)*(mult.110*(mult.110+1)/2)*(mult.011*(mult.011+1)/2)*(mult.001*(mult.001+1)/2)
      if(DE == "nonpositive")
      {
        mult.all = (mult.111)*(mult.101*(mult.101+1)/2)*(mult.000)*(mult.010*(mult.010+1)/2)*mult.001*(mult.011*(mult.011+1)/2)*mult.110*(mult.100*(mult.100+1)/2)
      }
      if(DE=="nonnegative")
      {
        mult.all = (mult.101)*(mult.111*(mult.111+1)/2)*(mult.010)*(mult.000*(mult.000+1)/2)*mult.011*(mult.001*(mult.001+1)/2)*mult.100*(mult.110*(mult.110+1)/2)
      }
    }
    
    ns.type = tapply(ns, num.id, mean)
    N.vars = sum(mult.all*(ns.type-1))
    index.symm = rep(1:nosymm, mult.all*(ns.type-1))
    n.types = (mult.all)
    mult.po= mult.all
    mult.do = mult.all
    n.po = mult.po
    n.do= mult.do
    n.per = cc
    Diff = rep(0, N.vars)
    Diff2 = Diff
    #V.list = vector("list", N.vars)
    row.ind = rep(0, 3*N.vars+1)
    col.ind = row.ind
    values = row.ind
    b = rep(0, nosymm + 2)
    for(kk in 1:nosymm)
    {
      row.ind[which(index.symm==kk)] = rep(kk, (ns.type[kk]-1)*n.types[kk])
      col.ind[which(index.symm==kk)] = which(index.symm==kk)  
      values[which(index.symm==kk)] = rep(1, (ns.type[kk]-1)*(n.types[kk]))
      b[kk] = n.per[kk]
    }
    row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
    col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
    row.ind[(2*N.vars+1):(3*N.vars)] = rep(nosymm + 2, N.vars)
    col.ind[(2*N.vars+1):(3*N.vars)] = 1:N.vars
    row.ind[(3*N.vars+1):(4*N.vars+1)] = rep(nosymm + 3, N.vars+1)
    col.ind[(3*N.vars+1):(4*N.vars+1)] = 1:(N.vars+1)
    PM = rep(0, N.vars)
    PV = rep(0, N.vars)
    
    
    zscore = rep(0, length(Gamma.vec))
    
    for(ee in 1:length(Gamma.vec))
    {
      Gamma.sens = Gamma.vec[ee]
      
      for(kk in 1:nosymm)
      {
        i = which(num.id==kk)[1]
        symmgroup = which(index.symm == kk)
        ind = which(index==i)
        treatstrat = treatment[ind]
        outstrat = outcome[ind]
        dosestrat = dose[ind]
        dosesymm = c(rep(0, sum(treatstrat==F&dosestrat==F)), rep(1, sum(treatstrat==F& dosestrat==T)), rep(0, sum(treatstrat==T&dosestrat==F)), rep(1, sum(treatstrat==T& dosestrat==T)))
        
        outsymm = c(sort(outstrat[treatstrat==F& dosestrat==F]), sort(outstrat[treatstrat==F& dosestrat==T]), sort(outstrat[treatstrat==T&dosestrat == F]), sort(outstrat[treatstrat==T&dosestrat == T]))
        PO = matrix(0, n.po[kk], ns[i])
        DO = matrix(0, n.po[kk],ns[i])
        
        mt.011 = mult.011[kk]
        mt.001 = mult.001[kk]
        mt.000 = mult.000[kk]
        mt.010 = mult.010[kk]
        mt.110 = mult.110[kk]
        mt.100 = mult.100[kk]
        mt.111 = mult.111[kk]
        mt.101 = mult.101[kk]
        p.011 = mt.011
        p.001 = mt.001
        p.000 = mt.000
        p.010 = mt.010
        p.110 = mt.110
        p.100 = mt.100
        p.111 = mt.111
        p.101 = mt.101
        
        D.000 = matrix(1, mt.000-1, mt.000-1)
        D.000[lower.tri(D.000)] = 0
        D.000 = rbind(D.000, c(rep(0, mt.000-1)) )
        D.010 = matrix(1, mt.010-1, mt.010-1)
        D.010[lower.tri(D.010)] = 0
        D.010 = rbind(D.010, c(rep(0, mt.010-1)) )
        D.001 = matrix(1, mt.001-1, mt.001-1)
        D.001[lower.tri(D.001)] = 0
        D.001 = rbind(D.001, c(rep(0, mt.001-1)) )
        D.011 = matrix(1, mt.011-1, mt.011-1)
        D.011[lower.tri(D.011)] = 0
        D.011 = rbind(D.011, c(rep(0, mt.011-1)))
        D.100 = matrix(1, mt.100-1, mt.100-1)
        D.100[lower.tri(D.100)] = 0
        D.100 = rbind(D.100, c(rep(0, mt.100-1)) )
        D.110 = matrix(1, mt.110-1, mt.110-1)
        D.110[lower.tri(D.110)] = 0
        D.110 = rbind(D.110, c(rep(0, mt.110-1)) )
        D.101 = matrix(1, mt.101-1, mt.101-1)
        D.101[lower.tri(D.101)] = 0
        D.101 = rbind(D.101, c(rep(0, mt.101-1)) )
        
        D.111 = matrix(1, mt.111-1, mt.111-1)
        D.111[lower.tri(D.111)] = 0
        D.111 = rbind(D.111, c(rep(0, mt.111-1)))
        
        if(MO == T)
        {
          D.100 = matrix(0, 1, mt.100-1)
          p.100 = 1
          D.110 = matrix(0, 1, mt.110-1)
          p.110 = 1
          D.011 = matrix(1, 1, mt.011-1)
          p.011 = 1
          D.001 = matrix(1, 1, mt.001-1)
          p.001 = 1
        }
        count = 1
        for(ll in 1:p.000)
        {
          for(mm in 1:p.010)
          {
            for(uu in 1:p.001)
            {
              for(vv in 1:p.011)
              {
                for(ww in 1:p.100)
                {
                  for(xx in 1:p.110)
                  {
                    for(yy in 1:p.101)
                    {
                      for(zz in 1:p.111)	
                      {
                        T.000 = matrix(1, p.000-ll, mt.000-ll)
                        T.000[lower.tri(T.000)] = 0
                        T.000 = rbind(T.000, c(rep(0, mt.000-ll)))
                        F.000 = matrix(0, (p.000-ll+1), (ll-1))
                        T.000 = cbind(F.000, T.000)
                        
                        T.010 = matrix(1, p.010-mm, mt.010-mm)
                        T.010[lower.tri(T.010)] = 0
                        T.010 = rbind(T.010, c(rep(0, mt.010-mm)))
                        F.010 = matrix(1, (p.010-mm+1), (mm-1))
                        T.010 = cbind(F.010, T.010)
                        
                        T.001 = matrix(1, uu-1, uu-1)
                        T.001[lower.tri(T.001)] = 0
                        T.001 = rbind(T.001, c(rep(0, uu-1)))
                        F.001 = matrix(0, (uu), (mt.001-uu))
                        T.001 = cbind(F.001, T.001)
                        if(MO == T)
                        {
                          T.001 = matrix(0, 1, mt.001-1)
                        }
                        
                        T.011 = matrix(1, vv-1, vv-1)
                        T.011[lower.tri(T.011)] = 0
                        T.011 = rbind(T.011, c(rep(0, vv-1)))
                        F.011 = matrix(1, (vv), (mt.011-vv))
                        T.011 = cbind(F.011, T.011)
                        if(MO == T)
                        {
                          T.011 = matrix(1, 1, mt.011-1)
                        }
                        
                        
                        T.100 = matrix(1, p.100-ww, mt.100-ww)
                        T.100[lower.tri(T.100)] = 0
                        T.100 = rbind(T.100, c(rep(0, mt.100-ww)))
                        F.100 = matrix(0, (p.100-ww+1), (ww-1))
                        T.100 = cbind(F.100, T.100)
                        if(MO == T)
                        {
                          T.100 = matrix(0, 1,mt.100-1)
                        }
                        
                        T.110 = matrix(1, p.110-xx, mt.110-xx)
                        T.110[lower.tri(T.110)] = 0
                        T.110 = rbind(T.110, c(rep(0, mt.110-xx)))
                        F.110 = matrix(1, (p.110-xx+1), (xx-1))
                        T.110 = cbind(F.110, T.110)
                        if(MO == T)
                        {
                          T.110 = matrix(1, 1, mt.110-1)
                        }
                        
                        T.101 = matrix(1, yy-1, yy-1)
                        T.101[lower.tri(T.101)] = 0
                        T.101 = rbind(T.101, c(rep(0, yy-1)))
                        F.101 = matrix(0, (yy), (mt.101-yy))
                        T.101 = cbind(F.101, T.101)
                        
                        T.111 = matrix(1, zz-1, zz-1)
                        T.111[lower.tri(T.111)] = 0
                        T.111 = rbind(T.111, c(rep(0, zz-1)))
                        F.111 = matrix(1, (zz), (mt.111-zz))
                        T.111 = cbind(F.111, T.111) 
                        
                        
                        
                        
                        
                        if(DE == "nonnegative")
                        {
                          T.100 = matrix(0, 1, mt.100-1)
                          T.101 = matrix(0, 1, mt.101-1)
                          T.011 = matrix(1, 1, mt.011-1)
                          T.010 = matrix(1, 1, mt.010-1)
                        }
                        if(DE == "nonpositive")
                        {
                          T.111 = matrix(1, 1, mt.111-1)
                          T.110 = matrix(1, 1, mt.110-1)
                          T.001 = matrix(0, 1, mt.001-1)
                          T.000 = matrix(0, 1, mt.000-1)
                          
                        }
                        
                        
                        
                        dosevec = c(D.000[ll,], D.010[mm,], D.001[uu,], D.011[vv,], D.100[ww,], D.110[xx,], D.101[yy,], D.111[zz,])
                        dosevec = dosevec[!is.na(dosevec)]
                        for(lll in 1:nrow(T.000))
                        {
                          for(mmm in 1:nrow(T.010))
                          {
                            for(uuu in 1:nrow(T.001))
                            {
                              for(vvv in 1:nrow(T.011))
                              {
                                for(www in 1:nrow(T.100))
                                {
                                  for(xxx in 1:nrow(T.110))
                                  {
                                    for(yyy in 1:nrow(T.101))
                                    {
                                      for(zzz in 1:nrow(T.111))	
                                      {
                                        DO[count,] = dosevec
                                        tempvec = c(T.000[lll,], T.010[mmm,], T.001[uuu,], T.011[vvv,], T.100[www,], T.110[xxx,], T.101[yyy,], T.111[zzz,])
                                        tempvec = tempvec[!is.na(tempvec)]
                                        PO[count,] = tempvec
                                        count = count+1
                                      }}}}}}}}}}}}}}}}             
        
        
        
        count = 1
        
        for(jj in 1:n.po[kk])
        {
          
          ind.jj = (((jj-1)*(ns[i]-1))+1):(jj*(ns[i]-1))
          po.symm = PO[jj,]
          do.symm = DO[jj,]
          treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
          outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
          outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm) 
          dosecontrol = dosesymm*(1-treatsymm) + do.symm*(treatsymm)
          dosetreat = dosesymm*(treatsymm) + do.symm*(1-treatsymm) 
          sum.cont = sum(outcontrol - null*dosecontrol)/(ns[i]-1)
          Q = (outtreat + outcontrol/(ns[i]-1) - null*dosetreat - null*dosecontrol/(ns[i]-1)- sum.cont)*ns[i]
          if(sum(treatstrat)>1)
          {
            sum.cont = sum(outtreat - null*dosetreat)/(ns[i]-1)
            Q = -(outtreat/(ns[i]-1) + outcontrol - null*dosetreat/(ns[i]-1) - null*dosecontrol- sum.cont)*ns[i]
            
          }
          
          
          
          qi = Q*max.e - Q*(!max.e)
          ord = order(qi)
          qi.sort = sort(qi)
          
          
          mu = rep(0, length(ind)-1)
          sigma2 = rep(0, length(ind)-1)
          
          
          for(j in 1:(length(ind)-1))
          {
            mu[j] = (sum(qi.sort[1:(j)]) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]))/((j) + Gamma.sens*(ns[i]-(j)))
            sigma2[j] = (sum(qi.sort[1:(j)]^2) + Gamma.sens*sum(qi.sort[(j+1):(length(ind))]^2))/((j) + Gamma.sens*(ns[i]-(j))) - mu[j]^2
          }
          mu[abs(mu) < 1e-8] = 0
          sigma2[sigma2 < 1e-8] = 0
          
          
          PM[symmgroup[ind.jj]] = mu*(max.e) - mu*(!max.e)
          PV[symmgroup[ind.jj]] = (sigma2)
          Diff[symmgroup[ind.jj]] = sum(outtreat -outcontrol - (null*(dosetreat - dosecontrol)))
          
          Diff2[symmgroup[count]] = sum(((dosetreat - dosecontrol)))
          count = count+1
          
        }
      }
      
      
      values[(N.vars+1):(2*N.vars)] = Diff
      values[(2*N.vars+1):(3*N.vars)] = Diff2
      values[(3*N.vars+1):(4*N.vars+1)] = c(-PM, 1)
      b[nosymm+1] = 0
      b[nosymm+2] = 1
      b[nosymm+3] = 0  
      alpha.opt = alpha
      if(alternative != "two.sided")
      {
        alpha.opt = 2*alpha
      }
      
      const.dir = c(rep("=", nosymm+1), ">=", "=")
      model = list()
      if(Gamma.sens==1)
      {
        model$A = sparseMatrix(row.ind[1:(3*N.vars)], col.ind[1:(3*N.vars)], x=values[1:(3*N.vars)])
        model$obj = c(PV)
        model$sense = const.dir[1:(nosymm+2)]
        model$rhs = b[1:(nosymm+2)]
        model$vtype = c(rep("I", N.vars))
        if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}
        
        
        model$modelsense = "max"
        
        
        solm = gurobi(model, params = list(OutputFlag = 0))
        zed = (TE.est/sqrt(solm$objval))
        x = solm$x[1:N.vars]
        SE = sqrt(solm$objval)
        kappavec[ee] = (TE.est)^2 - qchisq(1-alpha.opt, 1)*solm$objval
        tstat = zed
        pval = 0
        if(alternative == "two.sided")
        {
          pval = 2*pnorm(-abs(tstat))
        }
        if(alternative == "greater")
        {
          pval = 1 - pnorm((tstat))
        }
        if(alternative == "less")
        {
          pval = pnorm((tstat))
        }
        Reject = (pval < alpha)
      }
      if(Gamma.sens != 1)
      {
        diff = 200
        kappa = qchisq(1-alpha.opt, 1)
        while(diff > 1e-8)
        {
          Plin = -2*TE.est*PM - kappa*PV 
          rowind.q =  1:(N.vars+1)
          colind.q = 1:(N.vars+1)
          values.q = c(rep(0, N.vars),1)
          Q = sparseMatrix(rowind.q, colind.q, x=values.q)
          model$A = sparseMatrix(row.ind, col.ind, x=values)
          model$obj = c(Plin,0)
          model$Q = Q
          model$sense = const.dir
          model$rhs = b
          
          
          model$vtype = c(rep("I", N.vars), "C")
          if(continuous.relax == T){model$vtype = c(rep("C", N.vars+1))}
          model$lb = c(rep(0, N.vars), -Inf)
          
          
          model$modelsense = "min"
          
          
          solm = gurobi(model, params = list(OutputFlag = 0))
          x = solm$x[1:N.vars]
          kappa.new = (TE.est - sum(PM*x))^2/sum(PV*x)
          kappavec[ee] = (TE.est - sum(PM*x))^2 - qchisq(1-alpha.opt, 1)*sum(PV*x)
          diff = abs(kappa.new - kappa)
          pval = 0
          if(PVAL == F)
          {
            diff = 0
            Reject = kappa.new > kappa
          }
          kappa = kappa.new
        }
        zed = sqrt((TE.est - sum(PM*x))^2/sum(PV*x))
        
        if(alternative == "less")
        {
          zed = -zed
        }
        zscore[ee] = zed
        tstat = zed
        
        
        if(alternative == "two.sided")
        {
          pval = 2*pnorm(-abs(tstat))
        }
        if(alternative == "greater")
        {
          pval = 1 - pnorm((tstat))
        }
        if(alternative == "less")
        {
          pval = pnorm((tstat))
        }
        Reject = (pval < alpha)
        
        
        
        if(sign(TE.est - sum(PM*x))!=sign(TE.est))
        {
          Reject = F
          pval = 0.5
          
          if(alternative == "two.sided")
          {
            pval = 1
          }
        }
        
        if(alternative == "greater" & sum(PM*x) < 0)
        {
          pval = .5
        }
        if(alternative == "less" & sum(PM*x) > 0)
        {
          pval = .5
        }
        
      }
      pvalvec[ee] = pval
      Rejectvec[ee] = Reject 
    } 
  }
  if(PVAL == F)
  {
    return(list(Gamma.vec = Gamma.vec, Reject = Rejectvec, lambdahat = lambdahat, kappa = kappavec))
  }
  if(PVAL == T)
  {
    return(list(Gamma.vec = Gamma.vec, pval = pvalvec, lambdahat = lambdahat))
  }
  
  
}  

sensCRD2 = function(Gamma.vec, index, treatment, outcome, null = 0, DE = "both", alternative = "two.sided", alpha = 0.05, continuous.relax = F)
{
  sensCRD(index, treatment, outcome, null, DE, alternative, alpha, Gamma.vec = Gamma.vec, calculate.pval=F, continuous.relax)$kappa
}
sensCRR2  = function(Gamma.vec, index, treatment, outcome, null = 1, DE = "both", alternative = "two.sided", alpha = 0.05, continuous.relax = F)
{
  sensCRR(index, treatment, outcome, null, DE, alternative, alpha, Gamma.vec = Gamma.vec, calculate.pval = F, continuous.relax)$kappa
}
sensEffectRatio2  = function(Gamma.vec, index, treatment, outcome, dose, null = 0, DE = "both", MO = T, ER = T, alternative = "two.sided", alpha = 0.05, continuous.relax = F)
{
  sensEffectRatio(index, treatment, outcome, dose, null, DE, MO, ER, alternative, alpha, Gamma.vec = Gamma.vec, calculate.pval = F, continuous.relax)$kappa
}