#===== All estimator methods===


#===== AUC function:  fit using parametric models 
f.theta.par = function(U, D, X1, W.e, W.c, tau, dis.e, dis.c, e.IPW, e.DR, e.OR,...){
  f.theta.par.ph = function(U, D, X1, W.e, W.c, tau, dis.e, dis.c, e.IPW, e.DR, e.OR){
    #---- function
    {
      f.rand = function(n,w,eta,dis){
        p = length(eta)
        mn = w%*%matrix(eta[1:(p-1)],ncol=1)
        rsurvreg(n, mn, eta[p], dist = dis)}
      f.surv = function(t,w,eta,dis){
        p = length(eta)
        mn = w%*%matrix(eta[1:(p-1)],ncol=1)
        1-psurvreg(t, mn, eta[p], dist = dis)}
      f.dens = function(t,w,eta,dis){
        p = length(eta)
        mn = w%*%matrix(eta[1:(p-1)],ncol=1)
        dsurvreg(t, mn, eta[p], dist = dis)}
      f.chazd = function(t,w,eta,dis){
        p = length(eta)
        mn = w%*%matrix(eta[1:(p-1)],ncol=1)
        -log(1-psurvreg(t, mn, eta[p], dist = dis))}
      f.hazd = function(t,w,eta,dis){
        p = length(eta)
        mn = w%*%matrix(eta[1:(p-1)],ncol=1)
        dsurvreg(t, mn, eta[p], dist = dis)/(1-psurvreg(t, mn, eta[p], dist = dis))}
      f.beta = function(eta,dis){
        # beta is coefficients vector of 
        p=length(eta)
        if(dis=="exponential"){bet = -eta[-p]}
        if(dis=="weibull"){alph = 1/eta[p]; bet = -alph*eta[-p]}
        if(dis=="lognormal"){bet = eta[-p]}
        return(bet)}
    }
    
    
    n = length(U)
    W.e = cbind(1,W.e)
    W.c = cbind(1,W.c)
    
    #---IPW: AUC ----
    
    modl.c = survreg(Surv(U,!D)~0+W.c)
    eta.c = c(coef(modl.c), modl.c$scale)
    
    surv.c = f.surv(U, W.c, eta.c, dis.c)
    f.h = function(x1,x2){(x1>x2)+(x1==x2)/2}
    h.ij = outer(X1,X1,f.h)
    if (e.IPW){
      
      theta.IPW.1 = outer(D*(U<=tau)/surv.c,(U>tau)/f.surv(tau, W.c, eta.c,dis.c)) *h.ij
      theta.IPW.2 =  outer(D*(U<=tau)/surv.c,(U>tau)/f.surv(tau, W.c, eta.c,dis.c))
      theta.a.IPW = sum(theta.IPW.1, na.rm = T)/sum(theta.IPW.2, na.rm = T)
    }else{theta.a.IPW=NA}
    
    #---DR: AUC ----
    modl.e = survreg(Surv(U,D)~0+W.e)
    eta.e = c(coef(modl.e), modl.e$scale)
    surv.e.tau = f.surv(tau, W.e, eta.e,dis.e)
    surv.e = f.surv(U, W.e, eta.e, dis.e)
    
    if (e.DR){
      t1 = D/surv.c
      t2 = outer(t1*(U<=tau), t1*(U>tau))
      DR1.1 =  sum(t2*h.ij)
      DR2.1 =  sum(t2)
      
      f.int1 = function(j){ 
        f.t = function(t){
          surv.c.t = sapply(t, f.surv, W.c[j,], eta.c, dis.c)
          surv.e.t = sapply(t, f.surv, W.e[j,], eta.e, dis.e)
          pmin(surv.e.t, surv.e.tau[j])/surv.e.t /surv.c.t*sapply(t,f.hazd,w=W.c[j,], eta=eta.c, dis=dis.c)
        }
        integrate(f.t,0,U[j])$value #, abs.tol = 1e-9
      }
      f.int2 = function(j){ 
        f.t = function(t){
          surv.c.t = sapply(t, f.surv, W.c[j,], eta.c, dis.c)
          surv.e.t = sapply(t, f.surv, W.e[j,], eta.e, dis.e)
          (surv.e.t-surv.e.tau[j])/surv.e.t /surv.c.t*sapply(t,f.hazd,w=W.c[j,], eta=eta.c, dis=dis.c)
        }
        integrate(f.t,0,min(tau,U[j]))$value #, abs.tol = 1e-9
      }
      
      t1.2 = sapply(1:n, f.int1)
      t2.2 = sapply(1:n, f.int2)
      
      t1.1 = (1-D)*pmin(surv.e,surv.e.tau)/surv.e/surv.c
      t2.1 = (1-D)*(U<=tau)*(surv.e -surv.e.tau)/surv.e/surv.c
      
      int1 = t1.1-t1.2
      int2 = t2.1-t2.2
      
      tt1 = outer(t1*(U<=tau), int1)
      tt2 = outer(t1*(U>tau), int2)
      DR1.2 = sum(tt1*h.ij + tt2*t(h.ij))
      DR2.2 = sum(tt1 + tt2)   
      
      tt3 = outer(int2,int1)
      DR1.3 = sum(tt3*h.ij)
      DR2.3 = sum(tt3)
      
      theta.a.DR = sum(DR1.1+DR1.2+DR1.3)/sum(DR2.1+DR2.2+DR2.3)
      theta.a.DR.asy = sum(DR1.1+DR1.2)/sum(DR2.1+DR2.2)
      
    }else{theta.a.DR=theta.a.DR.asy =NA}
    
    #---OR: AUC ----
    if(e.OR){
      B.OR = (do.call(rbind,replicate(n,surv.e.tau,simplify = F))*
                (1-do.call(cbind,replicate(n,surv.e.tau,simplify = F))))
      diag(B.OR)=0
      A.OR = B.OR*h.ij
      theta.a.OR = sum(A.OR)/sum(B.OR)

    }else{theta.a.OR = NA}
    
    #---All ----
    res = c(a.OR = theta.a.OR,
            a.IPW = theta.a.IPW,
            a.DR = theta.a.DR)
    return(res)
  }
  theta.par = try(f.theta.par.ph(U, D, X1, W.e, W.c, tau, dis.e, dis.c, e.IPW, e.DR, e.OR),  silent = TRUE)
  if (class(theta.par)=="try-error"){ theta.par=rep(NA,3)} 
  return(theta.par)
}


#===== AUC function: fit using semi-parametric Cox models with the properties of PH model 
#---  f.theta.cox.ph.pli0
f.theta.cox.ph.pli0 = function(U, D, X1, W.e, W.c, tau, e.IPW, e.DR, e.OR, RK1,...){
  #---- function
  
  n = length(U)
  idx1 = order(U)
  U=U[idx1];D = D[idx1]; W.e = W.e[idx1,]; W.c = W.c[idx1,];X1 = X1[idx1]
  dat.c = data.frame(U,D,W.c)
  dat.e = data.frame(U,D,W.e)
  t1 = table(U)
  ties = as.numeric(rep(1:length(t1),t1))
  
  f.h = function(x1,x2){(x1>x2)+(x1==x2)/2}
  h.ij = outer(X1,X1,f.h)
  
  #---IPW: AUC ----
  modl.c = coxph(Surv(U,!D)~., data = dat.c)
  eta.c = matrix(coef(modl.c),ncol=1)
  if (e.IPW){
    surv.c.all = sapply(1:n, function(j){survfit(modl.c, newdata = W.c[j,])$surv})[ties,]
    surv.c = diag(surv.c.all)
    surv.c.tau = surv.c.all[sum(tau>=U),] #sum(tau>U)+1?
    
    theta.IPW.2 = outer(D*(U<=tau)/surv.c,(U>tau)/surv.c.tau)
    theta.IPW.1 = theta.IPW.2*h.ij
    theta.a.IPW = sum(theta.IPW.1, na.rm=T)/sum(theta.IPW.2, na.rm=T)
  }else{theta.a.IPW=NA}
  
  #---DR: AUC ----
  
  modl.e = coxph(Surv(U,D)~., data=dat.e)
  eta.e =  matrix(coef(modl.e),ncol=1)
  
  if (e.DR){
    #--- by pli (Piecewise linear interpolation)---
    linterp = function (x1, y1, x2, y2){
      m = (y2 - y1)/(x2 - x1)
      b = y2 - m * x2
      return(c(b, m))
    }
    f.pli = function(s, x, y){
      y = y[order(x)]
      x = x[order(x)]
      n = length(x)
      ylieq = c(y[1],y[n])
      fi = function(s){
        i = sum(s>=x)
        if(i==0){ p=matrix(c(ylieq[1],0),2)
        }else if(i==n){p=matrix(c(ylieq[2],0),2)
        }else{p = as.vector(linterp( x[i], y[i],  x[i+1], y[i+1]))}
        p
      }
      p = sapply(s,fi)
      tt = p[1,]+p[2,]*s
      return(as.vector(tt))
    }
    
    f.surv.pli = function(s, w0, eta, x, y){
      sr = f.pli(s, x, y)
      as.vector(sr^exp(as.vector(w0%*%eta)))
    }
    rk = RK1
    su = c(rep(U[-n],each=rk) + c(rep(diff(U)/rk,each=rk) *rep(0:(rk-1), n-1)),U[n])
 
    
    W.c0 = sweep(W.c, 2, colMeans(W.c))
    surv.c0= survfit(modl.c, newdata=colMeans(W.c))$surv[ties]
    surv.c.tau = f.surv.pli(tau, W.c0, eta.c, U, surv.c0)
    chazd.c0 = survfit(modl.c, newdata=colMeans(W.c))$cumhaz[ties] #,  lieq =list(matrix(c(0, 0),1)), monotone = 1)#  xki =c(0, max(U)),
    
    chazd.c1 = f.pli(su, U, chazd.c0)
    chazd.c = sapply(1:n, function(j){chazd.c1*exp(as.vector(W.c0[j,]%*%eta.c))})
    hazd.c.t = diff(rbind(0,chazd.c))
    
    W.e0 = sweep(W.e, 2, colMeans(W.e))
    surv.e0 = survfit(modl.e, newdata=colMeans(W.e))$surv[ties]
    surv.e.tau = f.surv.pli(tau, W.e0, eta.e, U, surv.e0)
    
    surv.e = (surv.e0)^exp(as.vector(W.e0%*%eta.e))
    surv.c = (surv.c0)^exp(as.vector(W.c0%*%eta.c))
    
    t1 = D/surv.c
    t2 = outer(t1*(U<=tau), t1*(U>tau))
    DR1.1 =  sum(t2*h.ij)
    DR2.1 =  sum(t2) 
    
    f.int = function(j){ 
      surv.e.t = f.surv.pli(su, w0=W.e0[j,], eta=eta.e, x=U, y=surv.e0)
      surv.c.t = f.surv.pli(su, w0=W.c0[j,], eta=eta.c, x=U, y=surv.c0)
      
      tj1 = pmin(surv.e.t, surv.e.tau[j])/surv.e.t /surv.c.t*hazd.c.t[,j]
      int1 = sum(tj1[su<=U[j]], na.rm = T)
      
      tj2 = (surv.e.t-surv.e.tau[j])/surv.e.t /surv.c.t* hazd.c.t[,j]
      int2 = sum(tj2[(su<=U[j])&(su<tau)], na.rm = T) #, abs.tol = 1e-9
      c(int1, int2)
    }
    
    t1.2 = sapply(1:n, f.int)
    
    t1.1 = (1-D)*pmin(surv.e,surv.e.tau)/surv.e/surv.c
    t2.1 = (1-D)*(U<=tau)*(surv.e -surv.e.tau)/surv.e/surv.c
    
    int1 = t1.1-t1.2[1,]; int2 = t2.1-t1.2[2,]
    
    tt1 = outer(t1*(U<=tau), int1); tt2 = outer(t1*(U>tau), int2)
    DR1.2 = sum(tt1*h.ij + tt2*t(h.ij))
    DR2.2 = sum(tt1 + tt2)   
    
    tt3 = outer(int2,int1)
    DR1.3 = sum(tt3*h.ij)
    DR2.3 = sum(tt3)
    
    theta.a.DR = (DR1.1+DR1.2+DR1.3)/(DR2.1+DR2.2+DR2.3)
    theta.a.DR.asy = (DR1.1+DR1.2)/(DR2.1+DR2.2)
  }else{theta.a.DR = NA}
  
  #---OR: AUC ----
  if(e.OR){
    surv.e.all = sapply(1:n, function(j){survfit(modl.e, newdata = W.e[j,])$surv})[ties,]
    surv.e.tau = surv.e.all[sum(tau>=U),] 
    
    B.OR = (do.call(rbind,replicate(n,surv.e.tau,simplify = F))*
              (1-do.call(cbind,replicate(n,surv.e.tau,simplify = F))))
    diag(B.OR)=0
    A.OR = B.OR*h.ij
    theta.a.OR = sum(A.OR)/sum(B.OR) 
  }else{theta.a.OR = NA}
  
  #---All ----
  res = c(a.OR = theta.a.OR,
          a.IPW = theta.a.IPW,
          a.DR = theta.a.DR)
  
  return(res)
}
#---  f.theta.cox.ph.pli1: using g^OR_ave() instead of g^DR_ave
f.theta.cox.ph.pli1 = function(U, D, X1, W.e, W.c, tau, e.IPW, e.DR, e.OR, RK1,...){
  #---- function
  n = length(U)
  idx1 = order(U)
  U=U[idx1];D = D[idx1]; W.e = W.e[idx1,]; W.c = W.c[idx1,];X1 = X1[idx1]
  dat.c = data.frame(U,D,W.c)
  dat.e = data.frame(U,D,W.e)
  t1 = table(U)
  ties = as.numeric(rep(1:length(t1),t1))
  
  f.h = function(x1,x2){(x1>x2)+(x1==x2)/2}
  h.ij = outer(X1,X1,f.h)
  
  #---IPW: AUC ----
  modl.c = coxph(Surv(U,!D)~., data = dat.c)
  eta.c = matrix(coef(modl.c),ncol=1)
  if (e.IPW){
    surv.c.all = sapply(1:n, function(j){survfit(modl.c, newdata = W.c[j,])$surv})[ties,]
    surv.c = diag(surv.c.all)
    surv.c.tau = surv.c.all[sum(tau>=U),] #sum(tau>U)+1?
    
    theta.IPW.2 = outer(D*(U<=tau)/surv.c,(U>tau)/surv.c.tau)
    theta.IPW.1 = theta.IPW.2*h.ij
    theta.a.IPW = sum(theta.IPW.1, na.rm=T)/sum(theta.IPW.2, na.rm=T)
  }else{theta.a.IPW=NA}
  
  #---DR: AUC ----
  
  modl.e = coxph(Surv(U,D)~., data=dat.e)
  eta.e =  matrix(coef(modl.e),ncol=1)
  
  if (e.DR){
    #--- by pli (Piecewise linear interpolation)---
    linterp = function (x1, y1, x2, y2){
      m = (y2 - y1)/(x2 - x1)
      b = y2 - m * x2
      return(c(b, m))
    }
    f.pli = function(s, x, y, ylieq){
      y = y[order(x)]
      x = x[order(x)]
      n = length(x)
      if(missing(ylieq))ylieq = c(y[1],y[n])
      fi = function(s){
        i = sum(s>=x)
        if(i==0){ p=matrix(c(ylieq[1],0),2)
        }else if(i==n){p=matrix(c(ylieq[2],0),2)
        }else{p = as.vector(linterp( x[i], y[i],  x[i+1], y[i+1]))}
        p
      }
      p = sapply(s,fi)
      tt = p[1,]+p[2,]*s
      return(as.vector(tt))
    }
    
    f.surv.pli = function(s, w0, eta, x, y, ylieq){
      sr = f.pli(s, x, y, ylieq)
      as.vector(sr^exp(as.vector(w0%*%eta)))
    }
    rk = RK1
    su = c(rep(U[-n],each=rk) + c(rep(diff(U)/rk,each=rk) *rep(0:(rk-1), n-1)),U[n])

    y.lieq=c(1,0)

    W.c0 = sweep(W.c, 2, colMeans(W.c))
    surv.c0= survfit(modl.c, newdata=colMeans(W.c))$surv[ties]
    surv.c.tau = f.surv.pli(tau, W.c0, eta.c, U, surv.c0, y.lieq)
    chazd.c0 = survfit(modl.c, newdata=colMeans(W.c))$cumhaz[ties] #,  lieq =list(matrix(c(0, 0),1)), monotone = 1)#  xki =c(0, max(U)),
    
    chazd.c1 = f.pli(su, U, chazd.c0)
    chazd.c = sapply(1:n, function(j){chazd.c1*exp(as.vector(W.c0[j,]%*%eta.c))})
    hazd.c.t = diff(rbind(0,chazd.c))
    
    W.e0 = sweep(W.e, 2, colMeans(W.e))
    surv.e0 = survfit(modl.e, newdata=colMeans(W.e))$surv[ties]
    surv.e.tau = f.surv.pli(tau, W.e0, eta.e, U, surv.e0, y.lieq)

    U.ind = seq(from = 1, by=rk, length.out=n) # index of data in su -1
    suL = su[c(1,U.ind[-1]-1)]
    surv.eL = f.surv.pli(suL, W.e0, eta.e, U, surv.e0, y.lieq)
    surv.e = (surv.e0)^exp(as.vector(W.e0%*%eta.e))
    surv.c = (surv.c0)^exp(as.vector(W.c0%*%eta.c))
    
    surv.e.t = sapply(1:n,function(j){f.surv.pli(su, w0=W.e0[j,], eta=eta.e, x=U, y=surv.e0, ylieq=y.lieq)})
    surv.eL.t = rbind(surv.e.t[1,],surv.e.t[-length(su),])
    surv.c.t = sapply(1:n,function(j){f.surv.pli(su, w0=W.c0[j,], eta=eta.c, x=U, y=surv.c0, ylieq=y.lieq)})
    
    t1 = D/surv.c
    t2 = outer((U<=tau), surv.c.tau)*h.ij + outer((U>tau), (1-surv.c.tau))*t(h.ij)
    DR1.1 =  sum(rowSums(t2)*t1) #/2/n/n
    
    t3 = outer((U<=tau), surv.c.tau) + outer((U>tau), (1-surv.c.tau))
    DR2.1 =  sum(rowSums(t3)*t1) #/2/n/n
    
    
    f.int = function(j){ 

      tj1 = pmin(surv.e.t[,j], surv.e.tau[j])/surv.eL.t[,j] /surv.c.t[,j] #*t1
      int1 = sum(tj1*hazd.c.t[,j]*(su<=U[j]), na.rm = T)
      
      tj2 = (surv.e.t[,j]-surv.e.tau[j])/surv.eL.t[,j] /surv.c.t[,j]*(su<=tau) #*t2
      int2 = sum(tj2* hazd.c.t[,j] *((su<=U[j])), na.rm = T)
      
      c(int1, int2, tj1[U.ind[j]], tj2[U.ind[j]])
    }
    
    t.int = sapply(1:n, f.int)
    
    t4 = colSums(h.ij*(1-surv.e.tau),na.rm=T)
    t5 = colSums(t(h.ij)*surv.e.tau,na.rm=T)

    DR1.2 = sum( ( (1-D)*t.int[3,]-t.int[1,])*t4 + ( (1-D)*t.int[4,]-t.int[2,])*t5,na.rm=T)
    
    t4 = sum((1-surv.e.tau),na.rm=T)
    t5 = sum(surv.e.tau,na.rm=T)
    DR2.2 = sum( ( (1-D)*t.int[3,]-t.int[1,])*t4 + ( (1-D)*t.int[4,]-t.int[2,])*t5,na.rm=T)

    theta.a.DR = (DR1.1+DR1.2)/(DR2.1+DR2.2)
  }else{theta.a.DR = NA}
  
  #---OR: AUC ----
  if(e.OR){
    surv.e.all = sapply(1:n, function(j){survfit(modl.e, newdata = W.e[j,])$surv})[ties,]
    surv.e.tau = surv.e.all[sum(tau>=U),] 
    
    B.OR = (do.call(rbind,replicate(n,surv.e.tau,simplify = F))*
              (1-do.call(cbind,replicate(n,surv.e.tau,simplify = F))))
    diag(B.OR)=0
    A.OR = B.OR*h.ij
    theta.a.OR = sum(A.OR)/sum(B.OR) 
  }else{theta.a.OR = NA}
  
  #---All ----
  res = c(a.OR = theta.a.OR,
          a.IPW = theta.a.IPW,
          a.DR = theta.a.DR)
  
  return(res)
}


#===== C-index function: fit using parametri models and original formula
f.Ctau.par = function(U, D, X1, W.e, W.c, tau, dis.e, dis.c, e.IPW, e.DR, e.OR,...){
  {
    f.rand = function(n,w,eta,dis){
      p = length(eta)
      mn = w%*%matrix(eta[1:(p-1)],ncol=1)
      rsurvreg(n, mn, eta[p], dist = dis)}
    f.surv = function(t,w,eta,dis){
      p = length(eta)
      mn = w%*%matrix(eta[1:(p-1)],ncol=1)
      1-psurvreg(t, mn, eta[p], dist = dis)}
    f.dens = function(t,w,eta,dis){
      p = length(eta)
      mn = w%*%matrix(eta[1:(p-1)],ncol=1)
      dsurvreg(t, mn, eta[p], dist = dis)}
    f.hazd = function(t,w,eta,dis){
      p = length(eta)
      mn = w%*%matrix(eta[1:(p-1)],ncol=1)
      dsurvreg(t, mn, eta[p], dist = dis)/(1-psurvreg(t, mn, eta[p], dist = dis))}
    f.beta = function(eta,dis){
      # beta is coefficients vector of 
      p=length(eta)
      if(dis=="exponential"){bet = -eta[-p]}
      if(dis=="weibull"){alph = 1/eta[p]; bet = -alph*eta[-p]}
      if(dis=="lognormal"){bet = eta[-p]}
      
      return(bet)}
  }
  
  f.Ctau.par.org = function(U, D, X1, W.e, W.c, tau, dis.e, dis.c, e.IPW, e.DR, e.OR,...){
    
    n = length(U)
    W.e = cbind(1,W.e)
    W.c = cbind(1,W.c)
    Ut = U
    Ut[Ut>tau] = tau
    
    
    f.h = function(x1,x2){(x1>x2)+(x1==x2)/2}
    h.ij = sapply(1:n, function(j){f.h(X1[j],X1)}) #column for fixed j
    
    #---IPW: C-index ----
    modl.c =  survreg(Surv(U,!D)~0+W.c)
    eta.c =c(coef(modl.c), modl.c$scale)
    surv.c = f.surv(U, W.c, eta.c, dis.c) # is diag(surv.c) in code based on  Cox model.
    
    if (e.IPW){
      Ds = D*(U<tau) # --D.star 
      
      theta.IPW.1 = sapply(1:n, function(i){Ds[i]*(U>U[i])*h.ij[,i]/surv.c[i]/f.surv(U[i], W.c, eta.c,dis.c)})
      theta.IPW.2 = sapply(1:n, function(i){Ds[i]*(U>U[i])/surv.c[i]/f.surv(U[i], W.c, eta.c, dis.c)})  
      theta.c.IPW = sum(theta.IPW.1, na.rm=T)/sum(theta.IPW.2, na.rm=T)
    }else{theta.c.IPW=NA}
    
    #---DR: C-index ----
    modl.e = survreg(Surv(U,D)~0+W.e)
    eta.e = c(coef(modl.e), modl.e$scale)
    surv.e = f.surv(U, W.e, eta.e, dis.e) # is diag(surv.e) in code based on  Cox model.
    ETA.e = f.beta(eta.e, dis.e)
    
    if (e.DR){
      A1.1j = sapply(1:n, function(j){mean( h.ij[,j]*f.surv(U[j], W.e, eta.e, dis.e) )*(U[j]<tau)})
      A1.2j = sapply(1:n, function(j){mean( h.ij[j,]*(1-f.surv(Ut[j], W.e, eta.e,dis.e)) )})
      a2 = D*(A1.1j+ A1.2j)/surv.c
      tm1= sum(a2, na.rm=T)
      
      c.tij = function(s,j){
        t1 = f.surv(s, W.e, eta.e,dis.e)*f.surv(s, W.e[j,], eta.e,dis.e)
        t2 = f.surv(tau, W.e, eta.e,dis.e)*f.surv(tau, W.e[j,], eta.e,dis.e)
        (t1-t2)/(1+exp(as.matrix(W.e)%*%ETA.e-c(matrix(W.e[j,],nrow=1)%*%ETA.e)))
      }
      B.1j.c = function(s,j,th.ij){
        if(length(s)==length(j)){
          surv.c.t = f.surv(s, W.c[j,], eta.c, dis.c)*f.surv(s, W.e[j,], eta.e, dis.e)
        }else{
          surv.c.t1 = sapply(s, f.surv, w= W.c[j,], eta = eta.c, dis = dis.c)
          surv.c.t2 = sapply(s, f.surv, w= W.e[j,], eta = eta.e, dis = dis.e)
          surv.c.t = surv.c.t1*surv.c.t2
        }
        
        x =rep(NA,length(s))
        for (k in 1:length(s)){
          ctj1 = c.tij(s[k],j)
          ctj2 = f.surv(s[k], W.e[j,], eta.e, dis.e)-f.surv(tau, W.e[j,], eta.e, dis.e)
          ctj3 = (1-f.surv(tau, W.e, eta.e, dis.e))* f.surv(tau, W.e[j,], eta.e, dis.e)
          ctj4 = (1-f.surv(tau, W.e, eta.e, dis.e))* f.surv(s[k], W.e[j,], eta.e, dis.e)
          
          if(s[k]<tau){
            t1 =  ctj2-ctj1+ctj3; t2 = ctj1
          }else{t1 =ctj4;  t2 = 0}
          x[k]=mean(th.ij[,j]*t2+ th.ij[j,]*t1)/surv.c.t[k]
        }
        x
      }
      
      
      f.B1 = function(j, th.ij){ 
        f.t = function(s){B.1j.c(s,j,th.ij)*sapply(s, f.hazd, w=W.c[j,], eta=eta.c, dis=dis.c)}
        integrate(f.t, 0, U[j])$value
      }
      
      int.B2.a = sapply(1:n, f.B1, th.ij=h.ij)
      tm2 = sum((1-D)*sapply(1:n, function(j){B.1j.c(U[j],j,th.ij=h.ij)})-int.B2.a)
      theta.c1.DR = (tm1+tm2)/2/n 
      
      #--
      A1.1j = sapply(1:n, function(j){mean( f.surv(U[j], W.e, eta.e, dis.e) )*(U[j]<tau)})
      A1.2j = sapply(1:n, function(j){mean( (1-f.surv(Ut[j], W.e, eta.e,dis.e)) )})
      a2 = D*(A1.1j+ A1.2j)/surv.c
      tm1= sum(a2, na.rm=T)
      tt = matrix(1,n,n)
      int.B2.a = sapply(1:n, f.B1, th.ij= tt )
      tm2 = sum((1-D)*sapply(1:n, function(j){B.1j.c(U[j],j,th.ij=tt)})-int.B2.a)
      theta.c2.DR = (tm1+tm2)/2/n #sum(tm1+(1-D)*sapply(1:n, function(j){B.1j.c(U[j],j)}) - int.B2.a)/2/n
      
      theta.c.DR = theta.c1.DR/theta.c2.DR
      
    }else{theta.c.DR = NA}
    
    if(e.OR){
      c.tj = function(j){
        t2 = f.surv(tau, W.e, eta.e,dis.e)*f.surv(tau, W.e[j,], eta.e,dis.e)
        (1-t2)/(1+exp(as.matrix(W.e)%*%ETA.e-c(matrix(W.e[j,],nrow=1)%*%ETA.e)))
      }
      A = sapply(1:n, c.tj)
      diag(A)=0
      theta.c.OR = sum(A*h.ij)/sum(A) #( sum(A*h.ij)/n/(n-1)  )/( sum(A)/n/(n-1)  )
    }else{theta.c.OR = NA}
    res = c(ctau.OR = theta.c.OR,
            ctau.IPW = theta.c.IPW,
            ctau.DR = theta.c.DR)
    return(res)
  }
  
  
  theta.par = try(f.Ctau.par.org(U, D, X1, W.e, W.c, tau, dis.e, dis.c, e.IPW, e.DR, e.OR,...),  silent = TRUE)
  
  if (class(theta.par)=="try-error"){ theta.par=rep(NA,3)} 
  return(theta.par)
}


#===== C-index function: fit using semi-parametric Cox-PH models and original formula
#---  f.Ctau.cox.ph.pli0
f.Ctau.cox.ph.pli0 = function(U, D, X1, W.e, W.c, tau, e.IPW, e.DR, e.OR, RK1, RK2, ...){
  {
    f.rand = function(n,w,eta,dis){
      p = length(eta)
      mn = w%*%matrix(eta[1:(p-1)],ncol=1)
      rsurvreg(n, mn, eta[p], dist = dis)}
    f.surv = function(t,w,eta,dis){
      p = length(eta)
      mn = w%*%matrix(eta[1:(p-1)],ncol=1)
      1-psurvreg(t, mn, eta[p], dist = dis)}
    f.dens = function(t,w,eta,dis){
      p = length(eta)
      mn = w%*%matrix(eta[1:(p-1)],ncol=1)
      dsurvreg(t, mn, eta[p], dist = dis)}
    f.chazd = function(t,w,eta,dis){
      p = length(eta)
      mn = w%*%matrix(eta[1:(p-1)],ncol=1)
      -log(1-psurvreg(t, mn, eta[p], dist = dis))}
    f.hazd = function(t,w,eta,dis){
      p = length(eta)
      mn = w%*%matrix(eta[1:(p-1)],ncol=1)
      dsurvreg(t, mn, eta[p], dist = dis)/(1-psurvreg(t, mn, eta[p], dist = dis))}
    f.beta = function(eta,dis){
      # beta is coefficients vector of 
      p=length(eta)
      if(dis=="exponential"){bet = -eta[-p]}
      if(dis=="weibull"){alph = 1/eta[p]; bet = -alph*eta[-p]}
      if(dis=="lognormal"){bet = eta[-p]}
      return(bet)}
  }
  
  n = length(U)
  idx1 = order(U)
  U=U[idx1];D = D[idx1]; W.e = W.e[idx1,]; W.c = W.c[idx1,];X1 = X1[idx1]
  
  nt = sum(U<tau)
  nt1= pmin(1:n,nt+1)
  
  dat.e = data.frame(U,D,W.e)
  dat.c = data.frame(U,D,W.c)
  t1 = table(U)
  ties = as.numeric(rep(1:length(t1),t1))
  
  f.h = function(x1,x2){(x1>x2)+(x1==x2)/2}
  h.ij = outer(X1,X1,f.h)
  
  #---IPW: C-index ----
  modl.c = coxph(Surv(U,!D)~., data = dat.c)
  eta.c = matrix(coef(modl.c),ncol=1)
  surv.c.all = sapply(1:n, function(j){survfit(modl.c, newdata = W.c[j,])$surv})[ties,]
  
  if (e.IPW){
    Ds = D*(U<tau) # --D.star
    theta.IPW.1 = sapply(1:n, function(i){Ds[i]*(U>U[i])*h.ij[i,]/surv.c.all[i,i]/surv.c.all[i,]}) #surv.c[i,j]=Sc(Ui|Wj)
    theta.IPW.2 = sapply(1:n, function(i){Ds[i]*(U>U[i])/surv.c.all[i,i]/surv.c.all[i,]})  
    theta.c.IPW = mean(theta.IPW.1, na.rm=T)/mean(theta.IPW.2, na.rm=T)
  }else{theta.c.IPW=NA}
  
  
  #---DR: C-index ----
  modl.e = coxph(Surv(U,D)~., data=dat.e)
  eta.e =  matrix(coef(modl.e),ncol=1)
  surv.e.all = sapply(1:n, function(j){survfit(modl.e, newdata = W.e[j,])$surv})[ties,]
  
  if (e.DR){
    linterp = function (x1, y1, x2, y2){
      m = (y2 - y1)/(x2 - x1)
      b = y2 - m * x2
      return(c(b, m))
    }
    
    f.pli = function(s, x, y){
      y = y[order(x)]
      x = x[order(x)]
      n = length(x)
      ylieq = c(y[1],y[n])
      fi = function(s){
        i = sum(s>=x)
        if(i==0){ p=matrix(c(ylieq[1],0),2)
        }else if(i==n){p=matrix(c(ylieq[2],0),2)
        }else{p = as.vector(linterp( x[i], y[i],  x[i+1], y[i+1]))}
        p
      }
      p = sapply(s,fi)
      tt = p[1,]+p[2,]*s
      return(as.vector(tt))
    }
    
    f.surv.pli = function(s, w0, eta, x, y){
      sr = f.pli(s, x, y)
      as.vector(sr^exp(as.vector(w0%*%eta)))
    }
    
    # y.lieq=c(1,0)
    {
      rk = RK1 
      
      su = c(rep(U[-n],each=rk) + c(rep(diff(U)/rk,each=rk) *rep(0:(rk-1), n-1)),U[n])
      ns = length(su)

      W.c0 = sweep(W.c, 2, colMeans(W.c))
      surv.c0= survfit(modl.c, newdata=colMeans(W.c))$surv[ties]
      surv.c.tau = f.surv.pli(tau, W.c0, eta.c, U, surv.c0)
      chazd.c0 = survfit(modl.c, newdata=colMeans(W.c))$cumhaz[ties] #,  lieq =list(matrix(c(0, 0),1)), monotone = 1)#  xki =c(0, max(U)),
      
      chazd.c1 = f.pli(su, U, chazd.c0)
      chazd.c = sapply(1:n, function(j){chazd.c1*exp(as.vector(W.c0[j,]%*%eta.c))})
      hazd.c.t = diff(rbind(0,chazd.c))
      
      W.e0 = sweep(W.e, 2, colMeans(W.e))
      surv.e0= survfit(modl.e, newdata=colMeans(W.e))$surv[ties]
      surv.e.tau = f.surv.pli(tau, W.e0, eta.e, U, surv.e0)
      
      surv.e = diag(surv.e.all) 
      surv.c = diag(surv.c.all)
      
      nt = sum(U<tau)

      t1 = D/surv.c
      t2 = outer(t1, t1)* outer(U, U, function(t1,t2){t1<pmin(t2,tau)})
      
      DR1.1 =  sum(t2*h.ij)
      DR2.1 =  sum(t2)
      
      
      surv.e.t =  sapply(1:n, function(j){f.surv.pli(su, w0=W.e0[j,], eta=eta.e, x=U, y=surv.e0)}) # St(u|Wj) for the j-column
      surv.c.t =  sapply(1:n, function(j){f.surv.pli(su, w0=W.c0[j,], eta=eta.c, x=U, y=surv.c0)}) # Sc(u|Wj) for the j-column
      
      zero = rep(0,ns)
      onen = rep(1,n) 
      onens = rep(1,ns)
      A.indx = outer(su,U,"<=")
      
      f.int1 = function(i){
        se.ij = outer(onens,surv.e.all[i,])
        colSums( (U[i]<tau)*pmin(se.ij,surv.e.t)/surv.e.t/surv.c.t*hazd.c.t*A.indx, na.rm=T) 
      }
      
      f.int2 = function(i){
        t1 = outer(su<min(U[i],tau),onen)
        se.ij = outer(onens,pmax(surv.e.all[i,],surv.e.tau))
        colSums(t1*(surv.e.t- se.ij)/surv.e.t/surv.c.t*hazd.c.t*A.indx, na.rm=T)
      }
      
      t1.2 = sapply(1:n,f.int1)
      t2.2 = sapply(1:n,f.int2)
      
      t1.1 = sapply(1:n, function(i){(1-D)*(U[i]<tau)*pmin(surv.e, surv.e.all[i,])/surv.e/surv.c})
      t2.1 = sapply(1:n, function(i){(1-D)*(U<=tau)*(U<U[i])*(surv.e -pmax(surv.e.tau,surv.e.all[i,]))/surv.e/surv.c})
      
      int1 = t(t1.1-t1.2); int2 = t(t2.1-t2.2)
      
      DR1.2 = sum( (int1*h.ij+int2*t(h.ij) ) *outer(D/surv.c,rep(1,n)) )
      DR2.2 = sum( (int1+int2) *outer(D/surv.c,rep(1,n)))
      
    }
    
    {
      rk = RK2 
      
      su = c(rep(U[-n],each=rk) + c(rep(diff(U)/rk,each=rk) *rep(0:(rk-1), n-1)),U[n])
      ns = length(su)

      W.c0 = sweep(W.c, 2, colMeans(W.c))
      surv.c0= survfit(modl.c, newdata=colMeans(W.c))$surv[ties]
      surv.c.tau = f.surv.pli(tau, W.c0, eta.c, U, surv.c0)
      chazd.c0 = survfit(modl.c, newdata=colMeans(W.c))$cumhaz[ties] #,  lieq =list(matrix(c(0, 0),1)), monotone = 1)#  xki =c(0, max(U)),
      
      chazd.c1 = f.pli(su, U, chazd.c0)
      chazd.c = sapply(1:n, function(j){chazd.c1*exp(as.vector(W.c0[j,]%*%eta.c))})
      hazd.c.t = diff(rbind(0,chazd.c))
      
      W.e0 = sweep(W.e, 2, colMeans(W.e))
      surv.e0= survfit(modl.e, newdata=colMeans(W.e))$surv[ties]
      surv.e.tau = f.surv.pli(tau, W.e0, eta.e, U, surv.e0)
      
      surv.e = diag(surv.e.all)#(surv.e0)^exp(as.vector(W.e0%*%eta.e))
      surv.c = diag(surv.c.all)#(surv.c0)^exp(as.vector(W.c0%*%eta.c))
      
      # ts = surv.e.all
      nt = sum(U<tau)
      # ts[U>tau,] = do.call(rbind, replicate(n-nt, surv.e.tau, simplify=F))
      
      
      surv.e.t =  sapply(1:n, function(j){f.surv.pli(su, w0=W.e0[j,], eta=eta.e, x=U, y=surv.e0)}) # St(u|Wj) for the j-column
      surv.c.t =  sapply(1:n, function(j){f.surv.pli(su, w0=W.c0[j,], eta=eta.c, x=U, y=surv.c0)}) # Sc(u|Wj) for the j-column
      
      zero = rep(0,ns)
      onen = rep(1,n) 
      onens = rep(1,ns)
      A.indx = outer(su,U,"<=")
      
      t1 = outer(onen,surv.e)-pmax(outer(onen,surv.e.tau),surv.e.all)
      t2 = outer(surv.e, (U<tau))*outer(U,U,">")*t1
      
      t3 = ( pmin(outer(onen, surv.e), surv.e.all) * pmin(outer(surv.e, onen), t(surv.e.all))
             -outer(surv.e.tau,surv.e.tau) )*outer(U<tau,U<tau)
      mn = c(matrix(W.e,nrow=n)%*%eta.e)
      t4 = t3/(1+exp(outer(mn,mn,"-")) )
      
      tt = (1-D)/surv.e/surv.c
      item1 = (t2+t4)*outer(tt,tt)
      
      f.item2 = function(i){
        t1 = surv.e.t*(surv.e[i]-pmax(surv.e.tau[i], surv.e.t[,i]))*(U[i]<tau)*(U[i]<su)
        
        t2 = ((pmin(surv.e[i],surv.e.t[,i])*pmin(outer(onens,surv.e.all[i,]), surv.e.t)
               -surv.e.tau[i]*outer(onens, surv.e.tau))*(su<tau)*(U[i]<tau))
        t3 = outer(onens, 1+exp( c(matrix(W.e,nrow=n)%*%eta.e) -c(matrix(W.e[i,],nrow=1)%*%eta.e) ) )
        t4 = t1+t2/t3
        
        (1-D[i])*colSums(t4/surv.e.t/surv.c.t/surv.e[i]/surv.c[i]*hazd.c.t*A.indx, na.rm=T)
      }
      item2 = sapply(1:n,f.item2)
      
      f.item3 = function(j){
        
        t1 = surv.e[j]*(surv.e.t- outer(onens,pmax(surv.e.all[j,] ,surv.e.tau)))*(su<min(tau,U[j]))
        
        t2 = ((pmin(surv.e[j],surv.e.t[,j])*pmin(outer(onens,surv.e.all[j,]), surv.e.t)
               -surv.e.tau[j]*outer(onens, surv.e.tau))*(su<tau)*(U[j]<tau))
        t3 = outer(onens, 1+exp( c(matrix(W.e[j,],nrow=1)%*%eta.e) -c(matrix(W.e,nrow=n)%*%eta.e)) )
        t4 = t1+t2/t3
        
        (1-D[j])*colSums(t4/surv.e.t/surv.c.t/surv.e[j]/surv.c[j]*hazd.c.t*A.indx, na.rm=T)
      }
      item3 = t(sapply(1:n,f.item3))
      
      f.item4 = function(i,j){
        # outer--(u,t)
        t0 = outer((su<=U[j])/surv.e.t[,j]/surv.c.t[,j]*hazd.c.t[,j],
                   (su<=U[i])/surv.e.t[,i]/surv.c.t[,i]*hazd.c.t[,i])
        t1 = outer(surv.e.t[,j], onens)*outer( -pmax(surv.e.tau[i], surv.e.t[,i]), surv.e.t[,i], "+")*outer(pmin(tau,su),su,">")
        t2 = (outer(surv.e.t[,i], surv.e.t[,i], pmin)*outer(surv.e.t[,j], surv.e.t[,j], pmin)-
                surv.e.tau[i]*surv.e.tau[j])*(outer(su, su, pmax)<tau)
        t3 = 1+exp(c(matrix(W.e[j,]-W.e[i,],nrow=1)%*%eta.e)) 
        sum(t0*(t1+t2/t3), na.rm =T)
      }
      item4 = sapply(1:n, function(i){sapply(1:n, f.item4,i=i)})
      
      t1 = t(item1-item2-item3+item4)
      DR2.3 = sum(t1)
      DR1.3 = sum(t1*h.ij)
    }
    
    # print(c(DR1.1,DR2.1,DR1.2,DR2.2,DR1.3,DR2.3))
    
    theta.c.DR =  (DR1.1+DR1.2+DR1.3)/ (DR2.1+DR2.2+DR2.3)

  }else{theta.c.DR =  NA}
  
  
  #---OR: C-index ----
  if(e.OR){
    surv.e = sapply(1:n, function(j){survfit(modl.e, newdata = W.e[j,])$surv})[ties,]
    d.dis.e = -diff(rbind(1,surv.e))
    
    A = B = matrix(NA,n,n)
    
    for(j in 1:n){
      B0 = A0 = matrix(NA,n,n)
      for (i in 1:n){
        B0[,i] = surv.e[,i]*d.dis.e[,j]
        A0[,i] = h.ij[j,i]*B0[,i]
      }
      A0[(nt+1):n,]=0
      B0[(nt+1):n,]=0
      A[,j] = colSums(A0)
      B[,j] = colSums(B0)
    }
    diag(A)=0
    diag(B)=0
    
    theta.c.OR = sum(A)/sum(B)
    
  }else{theta.c.OR = NA}
  
  #--- all------
  
  res = c(ctau.OR = theta.c.OR,
          ctau.IPW = theta.c.IPW,
          ctau.DR = theta.c.DR)
  
  return(res)
}
#---  f.Ctau.cox.ph.pli1: using g^OR_ave() instead of g^DR_ave
f.Ctau.cox.ph.pli1 = function(U, D, X1, W.e, W.c, tau, e.IPW, e.DR, e.OR, RK1, ...){
  
  n = length(U)
  idx1 = order(U)
  U=U[idx1];D = D[idx1]; W.e = W.e[idx1,]; W.c = W.c[idx1,];X1 = X1[idx1]
  
  nt = sum(U<tau)
  nt1= pmin(1:n,nt+1)
  
  dat.e = data.frame(U,D,W.e)
  dat.c = data.frame(U,D,W.c)
  t1 = table(U)
  ties = as.numeric(rep(1:length(t1),t1))
  
  f.h = function(x1,x2){(x1>x2)+(x1==x2)/2}
  h.ij = outer(X1,X1,f.h)
  
  #---IPW: C-index ----
  modl.c = coxph(Surv(U,!D)~., data = dat.c)
  eta.c = matrix(coef(modl.c),ncol=1)
  surv.c.all = survfit(modl.c, newdata = data.frame(W.c))$surv[ties,]

  if (e.IPW){
    Ds = D*(U<tau) # --D.star
    theta.IPW.1 = sapply(1:n, function(i){Ds[i]*(U>U[i])*h.ij[i,]/surv.c.all[i,i]/surv.c.all[i,]}) #surv.c[i,j]=Sc(Ui|Wj)
    theta.IPW.2 = sapply(1:n, function(i){Ds[i]*(U>U[i])/surv.c.all[i,i]/surv.c.all[i,]})  
    theta.c.IPW = mean(theta.IPW.1, na.rm=T)/mean(theta.IPW.2, na.rm=T)
  }else{theta.c.IPW=NA}
  
  
  #---DR: C-index ----
  modl.e = coxph(Surv(U,D)~., data=dat.e)
  eta.e =  matrix(coef(modl.e),ncol=1)
  surv.e.all = survfit(modl.e, newdata = data.frame(W.e))$surv[ties,]

  if (e.DR){
    linterp = function (x1, y1, x2, y2){
      m = (y2 - y1)/(x2 - x1)
      b = y2 - m * x2
      return(c(b, m))
    }
    
    f.pli = function(s, x, y, ylieq){
      y = y[order(x)]
      x = x[order(x)]
      n = length(x)
      if(missing(ylieq))ylieq = c(y[1],y[n])
      fi = function(s){
        i = sum(s>=x)
        if(i==0){ p=matrix(c(ylieq[1],0),2)
        }else if(i==n){p=matrix(c(ylieq[2],0),2)
        }else{p = as.vector(linterp( x[i], y[i],  x[i+1], y[i+1]))}
        p
      }
      p = sapply(s,fi)
      tt = p[1,]+p[2,]*s
      return(as.vector(tt))
    }
    
    f.surv.pli = function(s, w0, eta, x, y, ylieq){
      sr = f.pli(s, x, y, ylieq)
      as.vector(sr^exp(as.vector(w0%*%eta)))
    }
    
    y.lieq=c(1,0)
    
    rk = RK1 
    
    su = c(rep(U[-n],each=rk) + c(rep(diff(U)/rk,each=rk) *rep(0:(rk-1), n-1)),U[n])
    ns = length(su)
    # indx.u = (1+rk*(0:(n-1))) # the position of U in su 
    
    W.c0 = sweep(W.c, 2, colMeans(W.c))
    surv.c0= survfit(modl.c, newdata=colMeans(W.c))$surv[ties]
    surv.c.tau = f.surv.pli(tau, W.c0, eta.c, U, surv.c0)
    chazd.c0 = survfit(modl.c, newdata=colMeans(W.c))$cumhaz[ties] 

    chazd.c1 = f.pli(su, U, chazd.c0)
    chazd.c = sapply(1:n, function(j){chazd.c1*exp(as.vector(W.c0[j,]%*%eta.c))})
    hazd.c.t = diff(rbind(0,chazd.c))
    
    W.e0 = sweep(W.e, 2, colMeans(W.e))
    surv.e0= survfit(modl.e, newdata=colMeans(W.e))$surv[ties]
    surv.e.tau = f.surv.pli(tau, W.e0, eta.e, U, surv.e0)
    
    surv.e = diag(surv.e.all)
    surv.c = diag(surv.c.all)

    
    t1 = D/surv.c
    
    t3 = (U<tau)*surv.e.all+outer((U>=tau), surv.e.tau)
    t2 = h.ij*(U<tau)*surv.e.all + t(h.ij)*(1-t3)
    DR1.1 =  sum(rowSums(t2)*t1) #/2/n/n
    
    t2 = (U<tau)*surv.e.all + (1- t3)
    DR2.1 =  sum(rowSums(t2)*t1)
    
    
    
    mu = as.vector((W.e%*%eta.e))
    prp = exp( outer(-mu, mu,"+") )+1
    t6.1 = (U<tau)*(surv.e*surv.e.all-outer(surv.e.tau, surv.e.tau))/prp
    t6.2 =  surv.e - ( t6.1 + outer((U>=tau)*surv.e+(U<tau)*surv.e.tau, surv.e.tau) )
    
    t7 = (1-D)/surv.e/surv.c
    DR1.2 = sum( t7*rowSums(h.ij*t6.1 + t(h.ij)*t6.2) )
    DR2.2 = sum( t7*rowSums(t6.1 + t6.2) )
    
    
    surv.e.t =  sapply(1:n, function(j){f.surv.pli(su, w0=W.e0[j,], eta=eta.e, x=U, y=surv.e0)}) # St(u|Wj) for the j-column
    surv.c.t =  sapply(1:n, function(j){f.surv.pli(su, w0=W.c0[j,], eta=eta.c, x=U, y=surv.c0)}) # Sc(u|Wj) for the j-column
    
    f.int = function(i){ 

      t8 = surv.e.t[,i]*surv.c.t[,i]
      
      t6.1 = (su<tau)*(surv.e.t[,i]*surv.e.t- outer(rep(surv.e.tau[i],ns), surv.e.tau))/outer(rep(1,ns),prp[i,])
      
      t6.2 = surv.e.t[,i]-( t6.1 + outer((su>=tau)*surv.e.t[,i]+(su<tau)*surv.e.tau[i], surv.e.tau) )
      
      int1 = sum( rowSums((outer(rep(1,ns),h.ij[i,])*t6.1+outer(rep(1,ns),h.ij[,i])*t6.2)/t8)*(su<=U[i]) )
      int2 = sum( rowSums((t6.1 + t6.2)/t8)*(su<=U[i]) )
      return(c(int1,int2))
    }
    t9 = sapply(1:n, f.int)
    DR1.3 = sum(t9[1,])
    DR2.3 = sum(t9[2,])
    
    theta.c.DR = (DR1.1+DR1.2-DR1.3)/ (DR2.1+DR2.2-DR2.3)
    
  }else{theta.c.DR =  NA}
  
  
  #---OR: C-index ----
  if(e.OR){
    surv.e = sapply(1:n, function(j){survfit(modl.e, newdata = W.e[j,])$surv})[ties,]
    d.dis.e = -diff(rbind(1,surv.e))
    
    A = B = matrix(NA,n,n)
    
    for(j in 1:n){
      B0 = A0 = matrix(NA,n,n)
      for (i in 1:n){
        B0[,i] = surv.e[,i]*d.dis.e[,j]
        A0[,i] = h.ij[j,i]*B0[,i]
      }
      A0[(nt+1):n,]=0
      B0[(nt+1):n,]=0
      A[,j] = colSums(A0)
      B[,j] = colSums(B0)
    }
    diag(A)=0
    diag(B)=0
    
    theta.c.OR = sum(A)/sum(B)
    
  }else{theta.c.OR = NA}
  
  #--- all------
  
  res = c(ctau.OR = theta.c.OR,
          ctau.IPW = theta.c.IPW,
          ctau.DR = theta.c.DR)
  
  return(res)
}


