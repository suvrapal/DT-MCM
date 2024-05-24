# R code for data generation

mydata = function(n,alpha,beta1,beta2,cens,setting){ 
  z1 = rnorm(n,mean=0,sd=1) # z1 and z2 are the two continuous covariates 
  # generated independently from standard normal distributions
  z2 = rnorm(n,mean=0,sd=1)
  
  piz = rep(NA,n) # this is the uncured probability pi(z)
  if(setting==1){ # linear classification boundary
    piz = (exp(0.3-(5*z1)-(3*z2))/(1+exp(0.3-(5*z1)-(3*z2))))
  }
  if(setting==2){ # non-linear classification boundary
    piz =  (exp(0.3-(5*z1*z2)-(3*z1*z2))/(1+exp(0.3-(5*z1*z2)-(3*z1*z2))))
  }
  
  C = runif(n,0,cens) # censoring time 
  U = runif(n,0,1)
  Y = rep(NA,n) # observed lifetime
  D = rep(NA,n) # censoring indicator
  J = rep(NA,n) # cured indicator (J=0 implies cured)
  Sp = rep(NA,n) # overall (population) survival function
  S1 = rep(NA,n) # survival function of susceptible group
  for(i in 1:n){
    if(U[i]<= 1-piz[i]){
      Y[i] = C[i]
      D[i] = 0
      J[i] = 0
      S1[i] = exp(-((Y[i]/exp(-((beta1*z1[i])+(beta2*z2[i]))/alpha))^alpha))
      Sp[i] = (1-piz[i]) + (piz[i]*exp(-((Y[i]/exp(-((beta1*z1[i])
                                                     +(beta2*z2[i]))/alpha))^alpha)))
    }
    else{
      T =  rweibull(1,shape=alpha,scale= exp(-((beta1*z1[i])+(beta2*z2[i]))/alpha) )
      Y[i] = min(T,C[i])
      J[i] = 1
      S1[i] = exp(-((Y[i]/exp(-((beta1*z1[i])+(beta2*z2[i]))/alpha))^alpha))
      Sp[i] = (1-piz[i]) + (piz[i]*exp(-((Y[i]/exp(-((beta1*z1[i])
                                                     +(beta2*z2[i]))/alpha))^alpha)))
      if(Y[i]==C[i]){
        D[i] = 0
      }
      else{
        D[i] = 1
      }
    }
  }
  return(data.frame(Y,D,z1,z2,J,uncure=piz,S1,Sp=Sp))
} # function to generate data under 3 different classification boundaries


library(e1071)
library(survival)
library(rpart)

smsurv <-function(Time,Status,X,beta,w,model){    
  
  death_point <- sort(unique(subset(Time, Status==1)))
  if(model=='ph') coxexp <- exp((beta)%*%t(X[,-1]))  
  lambda <- numeric()
  event <- numeric()
  for(i in 1: length(death_point)){
    event[i] <- sum(Status*as.numeric(Time==death_point[i]))
    if(model=='ph')  temp <- sum(as.numeric(Time>=death_point[i])*w*drop(coxexp))
    if(model=='aft')  temp <- sum(as.numeric(Time>=death_point[i])*w)
    temp1 <- event[i]
    lambda[i] <- temp1/temp
  }
  HHazard <- numeric()
  for(i in 1:length(Time)){
    HHazard[i] <- sum(as.numeric(Time[i]>=death_point)*lambda)
    if(Time[i]>max(death_point))HHazard[i] <- Inf
    if(Time[i]<min(death_point))HHazard[i] <- 0
  }
  survival <- exp(-HHazard)
  list(survival=survival)
}


# EM code
em.dt.RC = function(Time,Status,X,X1,Z,Z1,offsetvar,uncureprob,beta,emmax,eps,data,testdata){     
  Time1<-testdata$Y
  Status1<-testdata$D
  
  w <- Status
  w1 <- Status1
  n <- length(Status)
  m <- length(Status1)
  
  s <- smsurv(Time,Status,X,beta,w,model="ph")$survival
  
  convergence<- 1000;i <-1
  while (convergence > eps & i < emmax){ 
    #print(i)
    survival<-drop(s^(exp((beta)%*%t(X[,-1]))))
    
    ## E step 
    w <- Status+(1-Status)*(uncureprob*survival)/((1-uncureprob)+uncureprob*survival)
    
    ## M step
    multipleuncureprob=matrix(1:5*n, nrow=n,ncol=5)
    for (j in 1:n){multipleuncureprob[j,]<-rbinom(5,size = 1,prob=w[j])}
    
    uncureprob1<-c(1,1)
    uncureprob2<-c(1,1)
    uncureprob3<-c(1,1)
    uncureprob4<-c(1,1)
    uncureprob5<-c(1,1)
    
    for (j in 1:n){uncureprob1[j]=multipleuncureprob[j,1]}
    for (j in 1:n){uncureprob2[j]=multipleuncureprob[j,2]}
    for (j in 1:n){uncureprob3[j]=multipleuncureprob[j,3]}
    for (j in 1:n){uncureprob4[j]=multipleuncureprob[j,4]}
    for (j in 1:n){uncureprob5[j]=multipleuncureprob[j,5]}
    
    for (j in 1:n){uncureprob1[j]=uncureprob1[j]}
    for (j in 1:n){uncureprob2[j]=uncureprob2[j]}
    for (j in 1:n){uncureprob3[j]=uncureprob3[j]}
    for (j in 1:n){uncureprob4[j]=uncureprob4[j]}
    for (j in 1:n){uncureprob5[j]=uncureprob5[j]}
    
    uncureprob1<-as.factor(uncureprob1)
    uncureprob2<-as.factor(uncureprob2)
    uncureprob3<-as.factor(uncureprob3)
    uncureprob4<-as.factor(uncureprob4)
    uncureprob5<-as.factor(uncureprob5)
    
    update_cureb<-c(1,1)
    update_pred<-c(1,1)
    
    daata=data.frame(uncureprob1,Z[,-1])
    obj3 <- tune.rpart(uncureprob1~z1+z2, data = daata, minsplit =c(11,20,25), cp = c(0.001,0.005,0.01))
    bc<-obj3$best.parameters[1]  
    bg<-obj3$best.parameters[2]
    
    mod1 <- rpart(formula = uncureprob1~z1+z2, data = data,method = "class", 
                  control = rpart.control(minsplit = bc[[1]],minbucket = round(bc[[1]]/3), cp = bg[[1]]),xval = 10,parms = list(split="gini"))
    
    cp.min1<-mod1$cptable[which.min(mod1$cptable[,"xerror"]),"CP"]
    tree1<-prune(mod1, cp=cp.min1)
    
    proba1 <- predict(tree1, newdata = data,type = "prob")
    cproba1 <- predict(tree1, newdata = testdata,type = "prob")
    
    update_cureb1<-c(1,1)
    update_pred1<-c(1,1)
    for (z in 1:n){update_cureb1[z]<-proba1[z,colnames(proba1)==1]}
    for (d in 1:m){update_pred1[d]<-cproba1[d,colnames(cproba1)==1]}
    uncureprob1<-as.numeric(as.character(uncureprob1))
    
    mod2 <- rpart(formula = uncureprob2~z1+z2, data = data,method = "class",
                  control = rpart.control(minsplit = bc[[1]],minbucket = round(bc[[1]]/3), cp = bg[[1]]),xval = 10,parms = list(split="gini"))
    
    cp.min2<-mod2$cptable[which.min(mod2$cptable[,"xerror"]),"CP"]
    tree2<-prune(mod2, cp=cp.min2)
    
    proba2 <- predict(tree2, newdata = data,type = "prob")
    cproba2 <- predict(tree2, newdata = testdata,type = "prob")
    
    update_cureb2<-c(1,1)
    update_pred2<-c(1,1)
    for (z in 1:n){update_cureb2[z]<-proba2[z,colnames(proba2)==1]}
    for (d in 1:m){update_pred2[d]<-cproba2[d,colnames(cproba2)==1]}
    uncureprob2<-as.numeric(as.character(uncureprob2))
    
    mod3 <- rpart(formula = uncureprob3~z1+z2, data = data,method = "class",
                  control = rpart.control(minsplit = bc[[1]],minbucket = round(bc[[1]]/3), cp = bg[[1]]),xval = 10,parms = list(split="gini"))
    
    cp.min3<-mod3$cptable[which.min(mod3$cptable[,"xerror"]),"CP"]
    tree3<-prune(mod3, cp=cp.min3)
    
    proba3 <- predict(tree3, newdata = data,type = "prob")
    cproba3 <- predict(tree3, newdata = testdata,type = "prob")
    
    update_cureb3<-c(1,1)
    update_pred3<-c(1,1)
    for (z in 1:n){update_cureb3[z]<-proba3[z,colnames(proba3)==1]}
    for (d in 1:m){update_pred3[d]<-cproba3[d,colnames(cproba3)==1]}
    uncureprob3<-as.numeric(as.character(uncureprob3))
    
    mod4 <- rpart(formula = uncureprob4~z1+z2, data = data,method = "class",
                  control = rpart.control(minsplit = bc[[1]],minbucket = round(bc[[1]]/3), cp = bg[[1]]),xval = 10,parms = list(split="gini"))
    cp.min4<-mod4$cptable[which.min(mod4$cptable[,"xerror"]),"CP"]
    tree4<-prune(mod4, cp=cp.min4)
    
    proba4 <- predict(tree4, newdata = data,type = "prob")
    cproba4 <- predict(tree4, newdata = testdata,type = "prob")
    
    update_cureb4<-c(1,1)
    update_pred4<-c(1,1)
    for (z in 1:n){update_cureb4[z]<-proba4[z,colnames(proba4)==1]}
    for (d in 1:m){update_pred4[d]<-cproba4[d,colnames(cproba4)==1]}
    uncureprob4<-as.numeric(as.character(uncureprob4))
    
    mod5 <- rpart(formula = uncureprob5~z1+z2, data = data,method = "class",
                  control = rpart.control(minsplit = bc[[1]],minbucket = round(bc[[1]]/3), cp = bg[[1]]),xval = 10,parms = list(split="gini"))
    
    cp.min5<-mod5$cptable[which.min(mod5$cptable[,"xerror"]),"CP"]
    tree5<-prune(mod5, cp=cp.min5)
    
    proba5 <- predict(tree5, newdata = data,type = "prob")
    cproba5 <- predict(tree5, newdata = testdata,type = "prob")
    
    update_cureb5<-c(1,1)
    update_pred5<-c(1,1)
    for (z in 1:n){update_cureb5[z]<-proba5[z,colnames(proba5)==1]}
    for (d in 1:m){update_pred5[d]<-cproba5[d,colnames(cproba5)==1]}
    uncureprob5<-as.numeric(as.character(uncureprob5))
    
    for (z in 1:n){update_cureb[z]<-(update_cureb1[z]+update_cureb2[z]+update_cureb3[z]+update_cureb4[z]+update_cureb5[z])/5}
    for (d in 1:m){update_pred[d]<-(update_pred1[d]+update_pred2[d]+update_pred3[d]+update_pred4[d]+update_pred5[d])/5}
    
    update_beta <- coxph(Surv(Time, Status)~X[,-1]+offset(log(w)),subset=w!=0, method="breslow")$coef
    update_s <-smsurv(Time,Status,X,beta,w,model="ph")$survival
    update_s1 <-smsurv(Time1,Status1,X1,beta,w1,model="ph")$survival
    
    convergence<-sum(c(mean(update_cureb)-mean(uncureprob),update_beta-beta,mean(update_s)-mean(s))^2)
    
    uncureprob <- update_cureb
    uncurepred <- update_pred
    beta <- update_beta 
    s<-update_s
    i <- i+1
    #print(i)
  }
  S1 = drop(s^(exp((beta)%*%t(X[,-1])))) # survival function of susceptible group
  Sp = (1-uncureprob)+(uncureprob*S1)
  em.dt <- list(latencyfit=beta,Uncureprob=uncureprob,Uncurepred=uncurepred,S0=s,S1=S1,Sp=Sp,tau=convergence,Mod1=mod1, Mod2=mod2, Mod3=mod3, Mod4=mod4,Mod5=mod5)
}