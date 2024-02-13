#=======Extra variables of Breast Data====


dat0 = read.csv("NKI_cleaned_short.csv")[,-c(1,2,10,15,16)]
t(t(names(dat0)))
str(dat0)

any(is.na(dat0))

sapply(dat0,is.integer)
apply(dat0[sapply(dat0,is.integer)],2,table)

dat=dat0
dat$posnodes = factor(2-((dat0$posnodes==0)+(dat0$posnodes<=3)))
dat$grade = factor(dat0$grade)
dat$angioinv = factor(dat0$angioinv)


dat$TIME = pmin(dat0$timerecurrence,dat0$survival)
dat$EVENT = (ifelse(dat0$survival>dat0$timerecurrence,T,F)|dat0$eventdeath)

#== here to determine whether chemo,hormonal and amputation included as baseline covariates in the analysis.
dat = subset(dat,select = -c(eventdeath,survival,timerecurrence))


# data "dat" include the following variables
# time: denoted by "U"
# even: denoted by "D"
# covariates : denoted by "W.e" for event time model, "W.c" for censoring time model
# biomr:  denoted by "X1"

dat = dat[order(dat$TIME),] #important!!!! for survfit$surv

U = dat$TIME #+ rnorm(nrow(dat),0,0.05)

D = dat$EVENT

W.org = subset(dat,select = -c(TIME,EVENT))

Y.c = Surv(U,!D)

W.c = model.matrix(~.,W.org)[,-1]

dat.c = data.frame(U,D,W.c)

Y.e = Surv(U,D)

W.e = model.matrix(~.,W.org)[,-1]

dat.e = data.frame(U,D,W.e)

n = length(U)

p.c = ncol(W.c)
p.e = ncol(W.e)






