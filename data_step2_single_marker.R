#=== run data analysis====
#have to put it with "data_step0_func.R", "data_step1_extra_data_Breast_include cheme et al.R" and "data_step1_extra_data_Heart.R" in the same working path

rm(list=ls())

set.seed(123456)

dis.e= dis.c = "weibull" # if use Weibull distribution in parametric models for T or C.
RK1=30 # RK1-1 is the number of interpolation between adjacent time points in single integrals.
RK2=1 # RK2-1 is the number of interpolation between adjacent time points in double integrals.


{ #---Analysis Breast Data---
  source("data_step1_extra_data_Breast_include cheme et al.R")
  xt="esr"
  X1 = -W.org$esr1 #biomarker
  tau = 10 #c(5,10)
}

# { #--- Analysis Hearst Data---
#   source("data_step1_extra_data_Heart.R")
#   xt = "Creatinine" #c("Sodium","Creatinine") #biomarker
#   if (xt == "Sodium"){X1 = -W.org$Sodium}else{X1 = W.org$Creatinine}
#   tau =180 #c(90,180)
# }

source("data_step0_func.R")


#==== estimate =====

theta1 = f.theta.par(U, D, X1, W.e, W.c, tau, dis.e=dis.e, dis.c=dis.c, e.IPW=T, e.DR=T, e.OR=T) 

theta2 = f.theta.cox.ph.pli0(U, D, X1, W.e, W.c, tau, e.IPW=T, e.DR=T, e.OR=T, su=su, RK1=RK1, dis.c=dis.c)

theta3 = f.Ctau.par(U, D, X1, W.e, W.c, tau, dis.e=dis.e, dis.c=dis.c, e.IPW=T, e.DR=T, e.OR=T)

theta4 = f.Ctau.cox.ph.pli0(U, D, X1, W.e, W.c, tau, e.IPW=T, e.DR=T, e.OR=T, RK1, RK2)

est = c(theta1,theta2,theta3,theta4)

#==== bootstrap=====

nboots = 200

f.boots.one.single = function(seed, U, D, X1, W.e, W.c, tau, f1, f2, f3, f4, RK1, RK2){
  set.seed(seed)
  dis.e=dis.c="weibull"
  # print(cat(seed, date()))
  
  library(survival)
  n=length(U)
  idx = sort(sample(1:n,n,replace=T))
  U.b=U[idx]; D.b=D[idx]; X1.b=X1[idx]; W.e.b=W.e[idx,]; W.c.b=W.c[idx,]
  res1 = f1(U.b, D.b, X1.b, W.e.b, W.c.b, tau, dis.e=dis.e, dis.c=dis.c, e.IPW=T, e.DR=T, e.OR=T)
  res2 = f2(U.b, D.b, X1.b, W.e.b, W.c.b, tau, e.IPW=T, e.DR=T, e.OR=T,su=su, RK1=RK1, dis.c=dis.c)
  res3 = f3(U.b, D.b, X1.b, W.e.b, W.c.b, tau, dis.e=dis.e, dis.c=dis.c, e.IPW=T, e.DR=T, e.OR=T, RK1=RK1)
  res4 = f4(U.b, D.b, X1.b, W.e.b, W.c.b, tau, e.IPW=T, e.DR=T, e.OR=T, RK1=RK1, RK2=RK2)
  return(c(res1,res2,res3,res4))
}

res = c()
for (iboots in 1: nboots){
  ires = f.boots.one.single(iboots, U=U, D=D, X1=X1, W.e=W.e, W.c=W.c,
                  tau=tau, f1 = f.theta.par, f2 = f.theta.cox.ph.pli1,
                  f3 = f.Ctau.par, f4 = f.Ctau.cox.ph.pli1, RK1=RK1, RK2=RK2)
  res = rbind(res, ires)
}



