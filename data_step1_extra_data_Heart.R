
dat0 = read.csv("S1Data.csv")

names(dat0)[9]="EF"

table(dat0[,"Event"])

any(is.na(dat0))

apply(dat0[,2:7],2,FUN=table)
summary(dat0[,c(1,8:13)])

any(is.na(dat0))

dat = dat0
dat$CPK = log(dat0$CPK)
dat$Pletelets = log(dat0$Pletelets)
dat$EF = factor((dat0$EF>45)+(dat0$EF>30))

table(dat$BP[dat$Event==1])/sum(dat$Event==1)
table(dat$BP[dat$Event==0])/sum(dat$Event==0)


# data "dat" include the following variables
# time: denoted by "U"
# even: denoted by "D"
# covariates : denoted by "W.e" for event time model, "W.c" for censoring time model
# biomr:  denoted by "X1"

dat = dat[order(dat$TIME),] #important!!!! for survfit$surv

U = dat$TIME #+ rnorm(nrow(dat),0,0.05)

D = dat$Event

W.org = subset(dat,select = -c(TIME,Event))

Y.c = Surv(U,!D)

W.c = model.matrix(~.,W.org)[,-1]

dat.c = data.frame(U,D,W.c)

Y.e = Surv(U,D)

W.e = model.matrix(~.,W.org)[,-1]

dat.e = data.frame(U,D,W.e)

n = length(U)

p.c = ncol(W.c)
p.e = ncol(W.e)



