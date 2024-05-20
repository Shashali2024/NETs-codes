library(survival)
library(randomForestSRC)
library(randomForest)
library(glmnet)
library(plsRcox)
library(superpc)
library(survcomp)
library(devtools)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(mixOmics)
library(dplyr)
library(tibble)
library(BART)
library(gbm)
allcohort<- lapply(allcohort,function(x){
  x[,-c(1:3)] <- scale(x[,-c(1:3)])
  return(x)})
##################################
#rm(result)
result <- data.frame()
est_data <- allcohort$TCGA
val_data_list <- allcohort
pre_var <- colnames(est_data)[-c(1:3)]
est_dd <- est_data[,c('futime','fustat',pre_var)]
val_dd_list <- lapply(val_data_list,function(x){x[,c('futime','fustat',pre_var)]})
rf_nodesize <- 5
seed <- 1
##################################
library(readxl)
setwd("D:\\生信\\machine learning\\2机器学习")
TCGA=read.table("TCGA-geoExpTime.txt",header=T,sep="\t",check.names=F,row.names=1)
CGGA_693=read.table("CGGA-geoExpTime.txt",header=T,sep="\t",check.names=F,row.names=1)
Rembrandt=read.table("RE-geoExpTime.txt",header=T,sep="\t",check.names=F,row.names=1)

TCGA=as.data.frame(TCGA)
CGGA_693=as.data.frame(CGGA_693)
Rembrandt=as.data.frame(Rembrandt)

allcohort <- list(TCGA,CGGA_693,Rembrandt)
names(allcohort ) <- c("TCGA","CGGA_693","Rembrandt")


saveRDS(allcohort, file="allcohort.rds")  
gene=read.table("gene.txt",header=T,sep="\t",check.names=F,row.names=1)
lassoGene=gene$name2
est_data <- allcohort$TCGA
pre_var <- colnames(TCGA)[-c(1:3)]
est_dd <- est_data[,c('futime','fustat',pre_var)]

signature=read.table("signature.txt",header=T,sep="\t",check.names=F,row.names=1)
TCGA2=allcohort$TCGA
TCGA <- TCGA[,c('futime','fustat',lassoGene)]
CGGA_693 <- CGGA_693[,c('futime','fustat',lassoGene)]
Rembrandt <- Rembrandt[,c('futime','fustat',lassoGene)]
TCGA1 <- TCGA[,gene2]


#### 1.RSF ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(futime,fustat)~.,data = est_dd,
             ntree = 1000,nodesize = rf_nodesize,  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance)}))%>%
  rownames_to_column('id')
cc$Model <- 'RSF'
plot(fit)
imp= as.data.frame(fit$importance) %>% rownames_to_column("id")
imp= arrange(imp, -imp[,2])
imp= imp%>% head(20)
lassoGene=imp$id
est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
fit2 <- rfsrc(Surv(futime,fustat)~.,data = est_dd2,
              ntree = 1000,nodesize = rf_nodesize,
              splitrule = 'logrank',
              importance = T,
              proximity = T,
              forest = T,
              seed = seed)
rs2 <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS2=predict(fit2,newdata = x)$predicted)})
cc2 <- data.frame(Cindex=sapply(rs2,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS2,x))$concordance)}))%>%
  rownames_to_column('id')
cc2$Model <- 'RSF'
result <- rbind(result,cc2)
write.table(result,file="RSF-result.txt",sep="\t",quote=F,col.names=T)  
##################################-arrange(dat, -dat[,1], -dat[,2])   
rm(result)
result <- data.frame()
#### 2.Enet ####
##################################

x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$futime,est_dd$fustat))
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}
##################################


#### 3.StepCox ####
##################################

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(as.numeric(futime,fustat))~.,est_dd),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('StepCox','[',direction,']')
  result <- rbind(result,cc)
}
##################################

#### 4.CoxBoost ####
##################################

set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[,'futime'],est_dd[,'fustat'],as.matrix(est_dd[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'futime'],est_dd[,'fustat'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'futime'],est_dd[,'fustat'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('CoxBoost')
result <- rbind(result,cc)
write.table(result,file="CoxBoost-result.txt",sep="\t",quote=F,col.names=T)   
##################################

#### 5.plsRcox####
##################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd[,pre_var],time=est_dd$futime,status=est_dd$fustat),nt=10,verbose = FALSE)
fit <- plsRcox(est_dd[,pre_var],time=est_dd$futime,event=est_dd$fustat,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('plsRcox')
result <- rbind(result,cc)

##################################


#### 6.superpc####
##################################

data <- list(x=t(est_dd[,-c(1,2)]),y=est_dd$futime,censoring.status=est_dd$fustat,featurenames=colnames(est_dd)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$futime,censoring.status=w$fustat,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('SuperPC')
result <- rbind(result,cc)
##################################


#### 7.GBM ####
##################################

set.seed(seed)
fit <- gbm(formula = Surv(futime,fustat)~.,data = est_dd,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(futime,fustat)~.,data = est_dd,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('GBM')
result <- rbind(result,cc)

##################################

#### 8.survivalsvm ####
##################################
fit = survivalsvm(Surv(futime,fustat)~., data= est_dd, gamma.mu = 1)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('survival-SVM')
result <- rbind(result,cc)
result2 <- result
result2$Model <- gsub('α','a',result2$Model)
##################################

#### 9.Ridge  ##### 
##################################
x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$futime,est_dd$fustat))

set.seed(seed)
fit = cv.glmnet(x1, x2,family = "cox",alpha=0,nfolds = 10)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('Ridge')
result <- rbind(result,cc)
##################################

#### 10.Lasso  ##### 
##################################
x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$futime,est_dd$fustat))

set.seed(seed)
fit = cv.glmnet(x1, x2,family = "cox",alpha=1,nfolds = 10)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('Lasso')
result <- rbind(result,cc)
##################################
#### 1.CoxBoost+survival-SVM ####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[,'futime'],est_dd[,'fustat'],as.matrix(est_dd[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'futime'],est_dd[,'fustat'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'futime'],est_dd[,'fustat'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
coef=coef(fit)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=names(coef)[index]
est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
#est_dd2 <- est_data[,c('futime','fustat',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
fit2 = survivalsvm(Surv(futime,fustat)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('CoxBoost',' + survival-SVM')
result <- rbind(result,cc)
##################################

#### 2.CoxBoost+SuperPC ####
##################################
data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$futime,censoring.status=est_dd2$fustat,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$futime,censoring.status=w$fustat,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('CoxBoost',' + SuperPC')
result <- rbind(result,cc)
##################################

#### 3.CoxBoost+plsRcox ####
##################################
pre_var2 <- colnames(est_dd2)[-c(1:2)]
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,pre_var2],time=est_dd2$futime,status=est_dd2$fustat),nt=10,verbose = FALSE)
fit2 <- plsRcox(est_dd2[,pre_var2],time=est_dd2$futime,event=est_dd2$fustat,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,type="lp",newdata=x[,-c(1,2)])))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('CoxBoost',' + plsRcox')
result <- rbind(result,cc)
##################################

#### 4.CoxBoost+Enet ####
##################################
x1 <- as.matrix(est_dd2[,pre_var2])
x2 <- as.matrix(Surv(est_dd2$futime,est_dd2$fustat))
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit2 = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit2$lambda.min)))})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('CoxBoost',' + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}
##################################

#### 5.CoxBoost+StepCox  ####
##################################
set.seed(seed)
for (direction in c("both", "backward", "forward")) {
  fit2 <- step(coxph(Surv(as.numeric(futime,fustat))~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit2,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('CoxBoost',' +StepCox','[',direction,']')
  result <- rbind(result,cc)
}
##################################

#### 6.CoxBoost+GBM  ####
##################################
fit2 <- gbm(formula = Surv(futime,fustat)~.,data = est_dd2,distribution = 'coxph',
            n.trees = 10000,
            interaction.depth = 3,
            n.minobsinnode = 10,
            shrinkage = 0.001,
            cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit2$cv.error)
fit2 <- gbm(formula = Surv(futime,fustat)~.,data = est_dd2,distribution = 'coxph',
            n.trees = best,
            interaction.depth = 3,
            n.minobsinnode = 10,
            shrinkage = 0.001,
            cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,x,n.trees = best,type = 'link')))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('CoxBoost','+GBM')
result <- rbind(result,cc)
##################################
write.table(cc,file="CoxBoost1.txt",sep="\t",quote=F,col.names=T)  

#### 1.Lasso+survival-SVM ####
##################################
pre_var <- colnames(est_data)[-c(1:3)]
x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$futime,est_dd$fustat))
set.seed(seed)
fit = cv.glmnet(x1, x2,family = "cox",alpha=1,nfolds = 10)
coef=coef(fit)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})

fit2 = survivalsvm(Surv(futime,fustat)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('Lasso',' + survival-SVM')
result <- rbind(result,cc)
##################################

#### 2.Lasso + SuperPC##### 
##################################

data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$futime,censoring.status=est_dd2$fustat,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$futime,censoring.status=w$fustat,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('Lasso',' + SuperPC')
result <- rbind(result,cc)
##################################

#### 3.Lasso+plsRcox ####
##################################

pre_var2 <- colnames(est_dd2)[-c(1:2)]
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,pre_var2],time=est_dd2$futime,status=est_dd2$fustat),nt=10,verbose = FALSE)
fit <- plsRcox(est_dd2[,pre_var2],time=est_dd2$futime,event=est_dd2$fustat,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('Lasso',' + plsRcox')
result <- rbind(result,cc)
##################################

#### 4.Lasso+CoxBoost ####
##################################

pen <- optimCoxBoostPenalty(est_dd2[,'futime'],est_dd2[,'fustat'],as.matrix(est_dd2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd2[,'futime'],est_dd2[,'fustat'],as.matrix(est_dd2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit2 <- CoxBoost(est_dd2[,'futime'],est_dd2[,'fustat'],as.matrix(est_dd2[,-c(1,2)]),
                 stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('Lasso',' + CoxBoost')
result <- rbind(result,cc)
##################################

#### 5.Lasso+StepCox  ####
##################################

for (direction in c("both", "backward", "forward")) {
  fit2 <- step(coxph(Surv(as.numeric(futime,fustat))~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit2,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('Lasso',' +StepCox','[',direction,']')
  result <- rbind(result,cc)
}
##################################
#### 6.Lasso+GBM  ####
##################################

fit2 <- gbm(formula = Surv(futime,fustat)~.,data = est_dd2,distribution = 'coxph',
            n.trees = 10000,
            interaction.depth = 3,
            n.minobsinnode = 10,
            shrinkage = 0.001,
            cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit2$cv.error)
set.seed(seed)
fit2 <- gbm(formula = Surv(futime,fustat)~.,data = est_dd2,distribution = 'coxph',
            n.trees = best,
            interaction.depth = 3,
            n.minobsinnode = 10,
            shrinkage = 0.001,
            cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,x,n.trees = best,type = 'link')))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('Lasso','+GBM')
result <- rbind(result,cc)
write.table(result,file="Lasso1.txt",sep="\t",quote=F,col.names=T)   
##################################

#### 1.StepCox+survival-SVM ####
##################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(as.numeric(futime,fustat))~.,est_dd),direction = direction)
  coef=coef(fit)
  lassoGene=names(coef)
  est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
  
  fit2 = survivalsvm(Surv(futime,fustat)~., data= est_dd2, gamma.mu = 1)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2, x)$predicted))})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('StepCox','[',direction,']',' + survival-SVM')
  result <- rbind(result,cc)
}
##################################
#### 2.StepCox+SuperPC  ####
##################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(as.numeric(futime,fustat))~.,est_dd),direction = direction)
  coef=coef(fit)
  lassoGene=names(coef)
  est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
  
  data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$futime,censoring.status=est_dd2$fustat,featurenames=colnames(est_dd2)[-c(1,2)])
  set.seed(seed)
  fit2 <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
  cv.fit <- superpc.cv(fit2,data,n.threshold = 20,#default 
                       n.fold = 10,
                       n.components=3,
                       min.features=5,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  rs <- lapply(val_dd_list2,function(w){
    test <- list(x=t(w[,-c(1,2)]),y=w$futime,censoring.status=w$fustat,featurenames=colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit2,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2],RS=rr)
    return(rr2)
  })
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('StepCox','[',direction,']',' + SuperPC')
  result <- rbind(result,cc)
}
##################################

#### 3.StepCox+plsRcox ####
##################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(as.numeric(futime,fustat))~.,est_dd),direction = direction)
  coef=coef(fit)
  lassoGene=names(coef)
  est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
  pre_var2 <- colnames(est_dd2)[-c(1:2)]
  cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,pre_var2],time=est_dd2$futime,status=est_dd2$fustat),nt=10,verbose = FALSE)
  fit2 <- plsRcox(est_dd2[,pre_var2],time=est_dd2$futime,event=est_dd2$fustat,nt=as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,type="lp",newdata=x[,-c(1,2)])))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('StepCox','[',direction,']',' + plsRcox')
  result <- rbind(result,cc)
}
##################################

#### 4.StepCox+Enet ####
##################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(as.numeric(futime,fustat))~.,est_dd),direction = direction)
  coef=coef(fit)
  lassoGene=names(coef)
  est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
  pre_var2 <- colnames(est_dd2)[-c(1:2)]
  x1 <- as.matrix(est_dd2[,pre_var2])
  x2 <- as.matrix(Surv(est_dd2$futime,est_dd2$fustat))
  for (alpha in seq(0,1,0.1)) {
    set.seed(seed)
    fit2 = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
    rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit2$lambda.min)))})
    cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
      rownames_to_column('id')
    cc$Model <- paste0('StepCox','[',direction,']','+Enet','[α=',alpha,']')
    result <- rbind(result,cc)
  }
}
##################################

#### 5.StepCox+GBM ####
##################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(as.numeric(futime,fustat))~.,est_dd),direction = direction)
  coef=coef(fit)
  lassoGene=names(coef)
  est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
  
  
  data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$futime,censoring.status=est_dd2$fustat,featurenames=colnames(est_dd2)[-c(1,2)])
  fit2 <- gbm(formula = Surv(futime,fustat)~.,data = est_dd2,distribution = 'coxph',
              n.trees = 10000,
              interaction.depth = 3,
              n.minobsinnode = 10,
              shrinkage = 0.001,
              cv.folds = 10,n.cores = 6)
  # find index for number trees with minimum CV error
  best <- which.min(fit2$cv.error)
  fit3 <- gbm(formula = Surv(futime,fustat)~.,data = est_dd2,distribution = 'coxph',
              n.trees = best,
              interaction.depth = 3,
              n.minobsinnode = 10,
              shrinkage = 0.001,
              cv.folds = 10,n.cores = 8)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit3,x,n.trees = best,type = 'link')))})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('StepCox','[',direction,']','+GBM')
  result <- rbind(result,cc)
}
##################################

#### 6.StepCox+CoxBoost ####
##################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(as.numeric(futime,fustat))~.,est_dd),direction = direction)
  coef=coef(fit)
  lassoGene=names(coef)
  est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
  pen <- optimCoxBoostPenalty(est_dd2[,'futime'],est_dd2[,'fustat'],as.matrix(est_dd2[,-c(1,2)]),
                              trace=TRUE,start.penalty=500,parallel = T)
  cv.res <- cv.CoxBoost(est_dd2[,'futime'],est_dd2[,'fustat'],as.matrix(est_dd2[,-c(1,2)]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  fit2 <- CoxBoost(est_dd2[,'futime'],est_dd2[,'fustat'],as.matrix(est_dd2[,-c(1,2)]),
                   stepno=cv.res$optimal.step,penalty=pen$penalty)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('StepCox','[',direction,']','+ CoxBoost')
  result <- rbind(result,cc)
}
##################################
write.table(result,file="StepCox1.txt",sep="\t",quote=F,col.names=T)   
#### 1.RSF+survival-SVM ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(futime,fustat)~.,data = est_dd,
             ntree = 1000,nodesize = rf_nodesize,  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance)}))%>%
  rownames_to_column('id')
imp= as.data.frame(fit$importance) %>% rownames_to_column("id")
imp= arrange(imp, -imp[,2])
imp= imp%>% head(20)
lassoGene=imp$id
est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
fit2 = survivalsvm(Surv(futime,fustat)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('RSF',' + survival-SVM')
result <- rbind(result,cc)
write.table(cc,file="RSF7.txt",sep="\t",quote=F,col.names=T)   
##################################
#### 2.RSF+SuperPC ####
##################################
data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$futime,censoring.status=est_dd2$fustat,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$futime,censoring.status=w$fustat,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('RSF',' + SuperPC')
result <- rbind(result,cc)
write.table(cc,file="RSF6.txt",sep="\t",quote=F,col.names=T)   
##################################
##################################
#### 3.RSF+plsRcox ####
##################################
pre_var2 <- colnames(est_dd2)[-c(1:2)]
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,pre_var2],time=est_dd2$futime,status=est_dd2$fustat),nt=10,verbose = FALSE)
fit2 <- plsRcox(est_dd2[,pre_var2],time=est_dd2$futime,event=est_dd2$fustat,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,type="lp",newdata=x[,-c(1,2)])))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('RSF',' + plsRcox')
result <- rbind(result,cc)
write.table(cc,file="RSF5.txt",sep="\t",quote=F,col.names=T)   
##################################
#### 4.RSF+Enet ####
##################################
x1 <- as.matrix(est_dd2[,pre_var2])
x2 <- as.matrix(Surv(est_dd2$futime,est_dd2$fustat))
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit2 = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit2$lambda.min)))})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('RSF',' + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
  write.table(cc,paste(alpha,".RSF.txt",sep=""),sep="\t",quote=F,col.names=T)   
}


##################################
#### 5.RSF+StepCox  ####
##################################
set.seed(seed)
for (direction in c("both", "backward", "forward")) {
  fit2 <- step(coxph(Surv(as.numeric(futime,fustat))~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit2,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('RSF',' +StepCox','[',direction,']')
  result <- rbind(result,cc)
  write.table(cc,paste(direction,".RSF.txt",sep=""),sep="\t",quote=F,col.names=T) 
}
##################################
#### 6.RSF+GBM  ####
##################################
fit2 <- gbm(formula = Surv(futime,fustat)~.,data = est_dd2,distribution = 'coxph',
            n.trees = 10000,
            interaction.depth = 3,
            n.minobsinnode = 10,
            shrinkage = 0.001,
            cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit2$cv.error)
fit2 <- gbm(formula = Surv(futime,fustat)~.,data = est_dd2,distribution = 'coxph',
            n.trees = best,
            interaction.depth = 3,
            n.minobsinnode = 10,
            shrinkage = 0.001,
            cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,x,n.trees = best,type = 'link')))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('RSF','+GBM')
result <- rbind(result,cc)
write.table(cc,file="RSF2.txt",sep="\t",quote=F,col.names=T)   
##################################
#### 7.RSF+CoxBoost ####
##################################

pen <- optimCoxBoostPenalty(est_dd2[,'futime'],est_dd2[,'fustat'],as.matrix(est_dd2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd2[,'futime'],est_dd2[,'fustat'],as.matrix(est_dd2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit2 <- CoxBoost(est_dd2[,'futime'],est_dd2[,'fustat'],as.matrix(est_dd2[,-c(1,2)]),
                 stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit2,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance[1])}))%>%
  rownames_to_column('id')
cc$Model <- paste0('RSF',' + CoxBoost')
result <- rbind(result,cc)
##################################
write.table(cc,file="RSF1.txt",sep="\t",quote=F,col.names=T)   
##################################


rm(result)
result <- data.frame()
#### 7.CoxBoost+RSF ####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[,'futime'],est_dd[,'fustat'],as.matrix(est_dd[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'futime'],est_dd[,'fustat'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'futime'],est_dd[,'fustat'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
coef=coef(fit)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=names(coef)[index]
est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
#est_dd2 <- est_data[,c('futime','fustat',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
pre_var2 <- colnames(est_dd2)[-c(1:2)]

fit <- rfsrc(Surv(futime,fustat)~.,data = est_dd2,
             ntree = 1000,nodesize = rf_nodesize,  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance)}))%>%
  rownames_to_column('id')
cc$Model <- paste0('CoxBoost','+RSF')
write.table(cc,file="CoxBoost+RSF.txt",sep="\t",quote=F,col.names=T)   

##################################

rm(result)
result <- data.frame()
#### 7.Lasso+RSF  ####
##################################
set.seed(seed)
x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$futime,est_dd$fustat))
fit = cv.glmnet(x1, x2,family = "cox",alpha=1,nfolds = 10)
coef=coef(fit)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
fit2 <- rfsrc(Surv(futime,fustat)~.,data = est_dd2,
              ntree = 1000,nodesize = rf_nodesize,  
              splitrule = 'logrank',
              importance = T,
              proximity = T,
              forest = T,
              seed = seed)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit2,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance)}))%>%
  rownames_to_column('id')
cc$Model <- paste0('Lasso','+RSF')
result <- rbind(result,cc)
write.table(cc,file="Lasso+RSF.txt",sep="\t",quote=F,col.names=T)   
##################################
rm(result)
result <- data.frame()
#### 7.StepCox+RSF ####
##################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(as.numeric(futime,fustat))~.,est_dd),direction = direction)
  coef=coef(fit)
  lassoGene=names(coef)
  est_dd2 <- est_data[,c('futime','fustat',lassoGene)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('futime','fustat',lassoGene)]})
  
  pre_var2 <- colnames(est_dd2)[-c(1:2)]
  set.seed(seed)
  fit2 <- rfsrc(Surv(futime,fustat)~.,data = est_dd2,
                ntree = 1000,nodesize = rf_nodesize,  
                splitrule = 'logrank',
                importance = T,
                proximity = T,
                forest = T,
                seed = seed)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit2,newdata = x)$predicted)})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(futime,fustat)~RS,x))$concordance)}))%>%
    rownames_to_column('id')
  cc$Model <- paste0('StepCox','[',direction,']','+RSF')
  result <- rbind(result,cc)
}
write.table(result,file="StepCox+RSF.txt",sep="\t",quote=F,col.names=T)   
##################################
