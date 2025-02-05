
library(KMsurv)
library(survival)
library(ggplot2)
library(survminer)


data("diabetic")
Diabetic <- diabetic
attach(Diabetic)

# Subdataset for no treatment
data_no_trt <- Diabetic[Diabetic$trt==0,]

# Subdataset for laser
data_laser <- Diabetic[Diabetic$trt==1,]



##########################################################

# bar plot treatment

# Creating dataframe
table1 <- table(trt)
table1
tabl<- c("No Treatment", "Laser")
Freq <- c(197,197)
DF1 <- data.frame(tabl, Freq)

# Plotting

barplot1<-ggplot(data=DF1, aes(x=tabl, y=Freq)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()+ggtitle("Treatments") +
  xlab("Tratment Type") + ylab("Frequency")
barplot1


# bar plot status

# Creating dataframe
table2 <- table(status)
table2
tabl2<- c("Censored", "Event")
Freq <- c(239,155)
Df_cen <- data.frame(tabl2, Freq)

# Plotting

barplot2<-ggplot(data=Df_cen, aes(x=tabl2, y=Freq)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()+ggtitle("Treatments") +
  xlab("Tratment Type") + ylab("Frequency")
barplot2


##########################################################
# KM Estimates and survival curve wrt treatment


Diabetic_KM_Fit <- survfit(Surv(time, status) ~ trt,data=Diabetic)

summary(Diabetic_KM_Fit)


# Survical Curve


surv_curve <- ggsurvplot(
  Diabetic_KM_Fit,                     # survfit object with calculated statistics.
  data = Diabetic,             # data used to fit survival curves.
  risk.table = FALSE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#0072B2", "#CC79A7"),
  xlim = c(0,80),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in Months",   # customize X axis label.
  break.time.by = 10,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =
    c("No Treatment", "Laser")    # change legend labels.
)


surv_curve



# Survival curve with Risk table and censor plot


risktbl_cens <- ggsurvplot(
  Diabetic_KM_Fit,                     # survfit object with calculated statistics.
  data = Diabetic,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("#0072B2", "#CC79A7"),
  xlim = c(0,80),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 10,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =
    c("No Treatment", "Laser")    # change legend labels.
)


risktbl_cens






##########################################################




# Plot two hazrd functions for two treatment groups

Cum_Haz_Curve <- ggsurvplot(Diabetic_KM_Fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#0072B2", "#CC79A7"),
           fun = "cumhaz",legend.labs =
             c("No Treatment", "Laser") )

Cum_Haz_Curve


##########################################################

# CI for median


km_fit_sur_linear <- survfit(Surv(time, status) ~ trt,data=Diabetic,conf.type="plain")
print(km_fit_sur_linear,  rmedian = "individual" )

km_fit_sur_log <- survfit(Surv(time, status) ~ trt,data=Diabetic,conf.type="log-log")
print(km_fit_sur_log,  rmedian = "individual" )

km_fit_sur_arcsin <- survfit(Surv(time, status) ~ trt,data=Diabetic,conf.type="arcsin")
print(km_fit_sur_arcsin,  rmedian = "individual" )


##########################################################

# CI for survival function

# USing written functions in tutorial 4

alpha=0.05
critical.Z=qnorm(1-alpha/2)
KMPCIestimator <- function (Z,delta)
{ UZ<-unique(Z) #Unique observed times;
N<-length(UZ)
UZ.order<-order(UZ)
UZ<-UZ[UZ.order] #sort data;
KM <- rep(0,N)
Y<-rep(0,N)
D<-rep(0,N)
D[1]<-sum(Z[delta==1]==UZ[1])
Y[1]<-sum(Z >= UZ[1])
KM[1] <- 1-D[1]/Y[1] #this is for right continuous value
for (i in 2: N){
  D[i]<-sum(Z[delta==1]==UZ[i])
  Y[i]<-sum(Z >= UZ[i])
  KM[i] <- KM[i-1]*(1-D[i]/Y[i])
}
#

#Calculate variance and standard error;
sigma2.s<-rep(0,N)
for (i in 1: N){
  ## sigma2.s[i]<-sum( (UZ<=UZ[i])*(D/(Y*(Y-D))) ) #old version
  sigma2.s[i]<-sum( (UZ[1:i]<=UZ[i])*(D[1:i]/(Y[1:i]*(Y[1:i]-D[1:i]))))
  ## Note the data is sorted by UZ;
  ## Using this to avoid NaN for times smaller than the largest observation;
}
KM.var<-KM^2*sigma2.s
KM.se<-sqrt(KM.var)
sigma.s<-sqrt(sigma2.s)

# 100(1-alpha)% log-transformed confidence interval;

theta<-exp((critical.Z*sigma.s)/log(KM))
###according to lecture note in chapter 4 
logtransL<-KM^(1/theta)
logtransU<-KM^theta
# 100(1-alpha)% arcsine-square-root-transformed confidence interval;
###according to lecture note in chapter 4 
ass<-asin(KM^(1/2))
zss<-critical.Z*sigma.s*(KM/(1-KM))^(1/2)
azssL<-ass-0.5*zss
azssU<-ass+0.5*zss
asrtransL<-(sin(pmax(0,azssL)))^2
asrtransU<-(sin(pmin(pi/2,azssU)))^2
Full<-data.frame(UZ,D,Y,KM,sigma2.s,KM.var,KM.se)
Reduced<-subset(Full, (Full$D>0))
PCI.full<-data.frame(UZ,D, logtransL, logtransU, asrtransL, asrtransU)
PCI.reduced<-subset(PCI.full, (PCI.full$D>0))
list(Full, Reduced, PCI.full, PCI.reduced)
}


# log and arcsin cI for no treatment
KMPCI.est_no_trt<-KMPCIestimator(data_no_trt$time, data_no_trt$status)
KMPCI.est_no_trt

# log and arcsin cI for laser
KMPCI.est_laser<-KMPCIestimator(data_laser$time, data_laser$status)
KMPCI.est_laser

##############################################
# Two sample test 
# Hypothesis
#h1(t): Hazard function for no treatment
#h2(t): Haard function for laser
#H0: h1(t) = h2(t) = h(t) for all t  vs.
#Ha: h1(t) = h2(t) for some t 


# Breslow method for handling ties and the Wald test

cox_fit1<-coxph(Surv(time, status)~trt, method=c("breslow"), data = Diabetic) 
summary(cox_fit1)

##log rank
logranktest<-survdiff(Surv(time,status)~trt,data=Diabetic,rho=0)
logranktest

#### Peto

test==4
RiskSet <- function (Z,delta,Z1,delta1, test)   #(Z1,delta1) represent obs. in a group;
{ UZ<-unique(Z[delta==1])      #Unique observed times;
N<-length(UZ)
UZ.order<-order(UZ)
UZ<-UZ[UZ.order]            #sort data;
KM <- rep(0,N)
KM.tilde<-rep(0,N)           #Used in Peto-Peto weights;
Y<-rep(0,N)
D<-rep(0,N)
Y1<-rep(0,N)
D1<-rep(0,N)
YY<-rep(0,N)
for (i in 1 : N){
  D[i]<-sum(Z[delta==1]==UZ[i])
  Y[i]<-sum(Z >= UZ[i])
  D1[i]<-sum(Z1[delta1==1]==UZ[i])
  Y1[i]<-sum(Z1 >= UZ[i])
  YY[i]<-Y[i]+1              #Used in Peto-Peto weights;
  if (i==1){
    KM[i]=1
    KM.tilde[i]=1-D[1]/YY[1]
  }
  else{
    KM[i]<-KM[i-1]*(1-D[i-1]/Y[i-1])              #left continuous version of K-M estimator for Fleiming-Harrington weights;
    KM.tilde[i]<-KM.tilde[i-1]*(1-D[i]/YY[i])   #Modified K-M estimator;
  }
}

D2<-D-D1
Y2<-Y-Y1
Y1DY<-Y1*(D/Y)
Y2DY<-Y2*(D/Y)
DY1DY<-D1-Y1DY

YDY=((Y-D)/((Y-1)*(Y!=1)+1*(Y==1)))*(Y!=1)+ 0*(Y==1)  
#This means if Y=1, (Y-D)/(Y-1)=0, which will not affect results;
Y1YYYD<-(Y1/Y)*(1-Y1/Y)*YDY*D
# print(UZ)
dyTable<-data.frame(UZ, Y1, D1, Y2, D2, Y, D, Y1DY, DY1DY, Y1YYYD, KM, KM.tilde)

# Change test method here, e.g. 1->log-rank test, 2->Gehan test;
if (test==1){
  W=1
  cat("The test is Log-rank test:", "\n")
}
if (test==2){
  W<- Y    
  cat("The test is Gehan test:", "\n")
}
if (test==3){
  W<- KM*(1-KM)       # KM*(1-KM)-> Fleming-Harrinnton test with p=1 and q=1;
  cat("The test is Fleming-Harrinnton test with p=1 and q=1:", "\n")
}
if (test==4){
  W<- KM.tilde          
  cat("The test is Peto-Peto test:", "\n")
}
if (test==5){
  W<- Y^(1/2)       
  cat("The test is Tarone and Ware  test:", "\n")
}


#Calculate Z(tau) values;

Z1Tau<-sum(W*(D1-Y1DY))
Z2Tau<-sum(W*(D2-Y2DY))
ZTau<-matrix(c(Z1Tau,Z2Tau),2,1)


# Calculate variance and covariance;

sigma11<-sum(W^2*(Y1/Y)*(1-Y1/Y)*YDY*D)
sigma22<-sum(W^2*(Y2/Y)*(1-Y2/Y)*YDY*D)
sigma12<--sum(W^2*(Y1/Y)*(Y2/Y)*YDY*D)
sigma21<-sigma12

varmatrix<-matrix(c(sigma11,sigma12, sigma21, sigma22), 2,2)

# Calculate Chi-Squred statistic and p-value using Z_1(tau);
chiSq<-Z1Tau^2/sigma11
pval<-1-pchisq(chiSq, 1)
testval<-c(chiSq,pval) 

list(dyTable, ZTau, varmatrix, testval)
}

Z<-Diabetic$time
delta<-Diabetic$status
Z1<-Diabetic$time[Diabetic$trt==0]          # Use no trt as trt=0;
delta1<-Diabetic$status[Diabetic$trt==0]

# Set  a value for test for a specific test;
# test=1-> Log-rank, 2-> Gehan, 3-> Fleming-Harrington, p=1, q=1, 4-> Peto-Peto;
# test=5-> Tarone and Ware;

Results<-RiskSet(Z,delta,Z1,delta1, test)  

##the test is peto test
cat( "\n \n")
print(Results[[1]])
cat("\n \n")
cat( "\n \n")
print(Results[[2]])
cat("\n \n")
cat("Z(tau)=", Results[[2]][1], "\n \n")
cat("\n \n")
cat("variance and covariance matrix, use K-1 values only:", "\n \n")
print(Results[[3]])
cat("\n \n")
cat("sigma_11^2=", Results[[3]][1,1], "\n \n")
cat("\n \n")
cat("test statistic and p-value:", "\n \n")
print(Results[[4]])
cat("\n \n")
cat("Chi-squared test statistic=", Results[[4]][1], "\n \n")
cat("p-value of the test       =", Results[[4]][2], "\n \n")

######################################################################
######################################################################
######################################################################
######################################################################
# Model development


#Function for the local test for possible confounders, adjusted for  Treatment type;
#coxph.fit is the object from coxph() and df is the df of the confounders;

#This local test is based on the Wald test;

coxph.localtest<-function(coxph.fit, df){
  coef<-coxph.fit$coef
  var<-coxph.fit$var
  loglik<-coxph.fit$loglik[2]
  p<-length(coef)
  AIC<- -2*loglik+2*p    #Using the formula on p.277;
  var.mat<-solve(as.matrix(var[(p-df+1):p, (p-df+1):p]))
  coe.mat<-matrix(c(coef[(p-df+1):p]), df, 1)
  WaldChiSq<-t(coe.mat)%*%var.mat%*%coe.mat
  pvalue<-1-pchisq(WaldChiSq,df)
  results<-c(df, WaldChiSq, pvalue, AIC)
  list(results)
}


options(width=60,length=200,digits=5)   # Set number of digits to be printed

Table1<-matrix(0,4,4)

coxphTable1.laser<-coxph(Surv(time, status)~trt+laser, method=c("breslow"), data = Diabetic)
coxphTable1.age<-coxph(Surv(time, status)~trt+age, method=c("breslow"), data = Diabetic)
coxphTable1.eye<-coxph(Surv(time, status)~trt+eye, method=c("breslow"), data = Diabetic)
coxphTable1.risk<-coxph(Surv(time, status)~trt+risk, method=c("breslow"), data = Diabetic)

Table1[1,]<-c(coxph.localtest(coxphTable1.laser, df=1)[[1]])
Table1[2,]<-c(coxph.localtest(coxphTable1.age, df=1)[[1]])
Table1[3,]<-c(coxph.localtest(coxphTable1.eye, df=1)[[1]])
Table1[4,]<-c(coxph.localtest(coxphTable1.risk, df=1)[[1]])


cat("Table 1: Local test for possible confounders, adjusted for  treatment type", "\n")
cat("    DF, Wald Chi-Squre, p-value, AIC:", "\n")
print(Table1)



# risk has lowest p value hence include to the model 
# now Local test for possible confounders, adjusted for treatment type and risk




Table2<-matrix(0,3,4)

coxphTable2.laser<-coxph(Surv(time, status)~trt+risk+laser, method=c("breslow"), data = Diabetic)
coxphTable2.age<-coxph(Surv(time, status)~trt+risk+age, method=c("breslow"), data = Diabetic)
coxphTable2.eye<-coxph(Surv(time, status)~trt+risk+eye, method=c("breslow"), data = Diabetic)


Table2[1,]<-c(coxph.localtest(coxphTable2.laser, df=1)[[1]])
Table2[2,]<-c(coxph.localtest(coxphTable2.age, df=1)[[1]])
Table2[3,]<-c(coxph.localtest(coxphTable2.eye, df=1)[[1]])



cat("Table 3: Local test for possible confounders, adjusted for  treatment type and risk", "\n")
cat("    DF, Wald Chi-Squre, p-value, AIC:", "\n")
print(Table2)



# eye has lowest p value hence include to the model 
# now Local test for possible confounders, adjusted for treatment type,risk and eye


Table3<-matrix(0,2,4)

coxphTable3.laser<-coxph(Surv(time, status)~trt+risk+eye+laser, method=c("breslow"), data = Diabetic)
coxphTable3.age<-coxph(Surv(time, status)~trt+risk+eye+age, method=c("breslow"), data = Diabetic)


Table3[1,]<-c(coxph.localtest(coxphTable3.laser, df=1)[[1]])
Table3[2,]<-c(coxph.localtest(coxphTable3.age, df=1)[[1]])


cat("Table 3: Local test for possible confounders, adjusted for  treatment type and risk", "\n")
cat("    DF, Wald Chi-Squre, p-value, AIC:", "\n")
print(Table3)


##########################Final Model#########################
Interaction_Model<-coxph(Surv(time, status)~trt+risk+eye+risk*eye, method=c("breslow"), data = Diabetic)

print(Interaction_Model) 



Final_Model<-coxph(Surv(time, status)~trt+risk+eye, method=c("breslow"), data = Diabetic)

cat("FInal Model: Analysis of Variance Table for the Final Model for \  
       diabetic retinopathy", "\n")
print(Final_Model) 











################### Cox PH assumption Checking############

Diabetic <- diabetic
#To include time-dependent covariates in the Cox PH model, we use the following function
#expand.breakpoints () to expand the data set first;
rm(list=ls())
options(width=60,length=200,digits=5) # Set number of digits to be printed
expand.breakpoints <-
  function(dataset, index = "patnum", status = "status", tevent = "time",
           Zvar = F, breakpoints = NULL){
    # Expand <dataset> so that each individual has a row for each
    # unique failure time that is <= his/her/its own failure time.
    #
    # Original version written by Kathy Russell, for
    # Kenneth Hess, Department of Biomathematics,
    # University of Texas M. D. Anderson Cancer Centre.
    #
    # Substantially modified by JHM, Nov 23, 1998.
    #
    # ERROR checking
    onceonly <- paste("must be a character string matching exactly",
                      "one name in the first parameter, dataset")
    # Require dataset to be of type data.frame
    if((missing(dataset) || !is.data.frame(dataset)))
      stop("\n Parameter, dataset, must be a data frame")
    varset <- names(dataset)
    covset <- varset
    lvset <- 1:length(varset) # Require dataset to have unique names
    if(any(duplicated(varset))) {
      stop(paste("\n Parameter, dataset, must have uniquely defined",
                 "column names"))
    }
    # Require index to match exactly 1 dataset name
    if(length((indexloc <- lvset[varset == index])) != 1)
      stop(paste("\n Parameter, index,", onceonly))
    covset <- covset[covset != index]
    # Require status to match exactly 1 dataset name
    if(length((statusloc <- lvset[varset == status])) != 1)
      stop(paste("\n Parameter, status,", onceonly))
    covset <- covset[covset != status]
    # Require tevent to match exactly 1 dataset name
    if(length((teventloc <- lvset[varset == tevent])) != 1)
      stop(paste("\n Parameter, tevent,", onceonly))
    covset <- covset[covset != tevent] # -----------------------------
    # Form vector of breakpoints, if necessary
    if(is.null(breakpoints)) {
      times <- dataset[, tevent]
      breakpoints <- sort(unique(times[dataset[, status] == 1]))
    }
    # Require breakpoints to be a vector of length >= 1
    if((is.null(breakpoints)) || (!(is.numeric(
      breakpoints))) || (!(is.vector(
        breakpoints))) || (length(breakpoints) < 1)) stop(paste(
          "\n Parameter, breakpoints, must be a numeric vector",
          "with at least 1 element")) #*****************
    #Begin
    #*****************
    n <- nrow(dataset)
    temp.stop <- c(breakpoints, 0) # The 0 is a place-filler
    if(breakpoints[1] > 0)
      temp.start <- c(0, breakpoints)
    else temp.start <- c(-1, breakpoints)
    n.break <- length(breakpoints) + 1
    t.event <- dataset[, tevent] # ---------------------------
    ## Begin calculate n.epochs
    n.epochs <- sapply(1:n, function(m, breaks, times)
      sum(breaks < times[m]) + 1, breaks = breakpoints, times = t.event)
    # End n.epochs
    id <- rep(1:n, n.epochs) # Index into supplied dataset
    last.epoch <- cumsum(n.epochs)
    epoch <- unlist(sapply(n.epochs, function(x)
      1:x)) # Index into vectors of interval start & stop points
    if(Zvar) {
      Zmat <- diag(rep(1, n.break))[epoch, ]
      Z.name <- paste("Z", seq(1, n.break), sep = "")
    }
    Tstop <- temp.stop[epoch]
    Tstop[last.epoch] <- dataset[, tevent]
    status2 <- rep(0, length(epoch))
    status2[last.epoch] <- dataset[, status]
    new.dataset <- data.frame(dataset[id, index, drop = F], temp.start[epoch],
                              Tstop, status2, epoch, dataset[id, covset, drop = F])
    if(Zvar)
      new.dataset <- data.frame(new.dataset, Zmat)
    nam <- c(index, "Tstart", "Tstop", status, "epoch", covset)
    if(Zvar)
      nam <- c(nam, Z.name)
    dimnames(new.dataset) <- list(1:length(epoch), nam)
    return(new.dataset)
  }

#We are going to add the time-dependent covariate Z_p(t) in the above final model in Table 8.9, to obta#Use the R function expand.breakpoints () to expand the data set first;
#An example to use the R function expand.breakpoints ();
#patnum<-seq(1:nobs) #Patients's ID;
#OX<-X%*%gamma #gamma1*baseline+gamma2*processor

Diabetic <- diabetic
Diabetic 
nobs<-dim(diabetic)[1]
patnum<-1:nobs
##Need more covariates in new.data;
new.data<-data.frame(time=Diabetic$time, delta=Diabetic$status, patnum, trt)
data.expand<-expand.breakpoints(new.data, index="patnum", status="delta", tevent="time", Zvar=F)
#Create time-dependent covariate Z2t=Z1*ln(t); Z1=one of fixed factors here;
lnt<-log(data.expand$Tstop)
Z1.t<-data.expand$trt*lnt

coxph.model<-coxph(Surv(Tstart,Tstop, delta)~trt+Z1.t,
                   data=data.expand,method=c("breslow"))
summary(coxph.model)



# Assumption Satisfied



#################### Model Diagnostics######################

fit1<-coxph(Surv(time, status)~trt+risk+eye, method=c("breslow"), data = Diabetic)


# COX snell diagnostic plot


library(survival)

mres = resid(fit1, type="martingale")
# type can be martingale, deviance, score, schoenfeld, dfbeta, dfbetas,
# scaledsch, partial
csres = status-mres
r.surv = survfit(Surv(csres,status)~1,type="fleming-harrington")
#pdf("Figure11_1.pdf")
par(las=1, mfrow=c(1,1), mai = c(0.5,1.0,1,0.1), omi = c(1.0,0,0.5,0))
plot(0,0,lty=1,type='n',xlim = c(0,1.5),ylim = c(0,1.5),xlab="Residual",
     ylab="Estimated Cum Hazards", main="Cox Snell diagnostic plot")
box()
lines(r.surv$time, -log(r.surv$surv), type='s')
lines(c(0,3),c(0,3))




##################### Martingale Residuals for Form of Covariates

plot(0,0,lty=1,type='n',xlim = c(0,20),ylim = c(-1.3,1),xlab="Waiting time", 
     ylab="Martingale residuals")
box()
points(risk, mres)
lines(lowess(mres~risk)$x,lowess(mres~risk)$y)









