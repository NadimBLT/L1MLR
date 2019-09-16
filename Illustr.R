path_repository="~/L1MLR-master/"

library(Rcpp)
library(glmnet)
library(Matrix)

sourceCpp(paste0(path_repository,"MultinomLogist_StdLasso_rcpp.cpp"))
source(paste0(path_repository,"MultinomLogist_StdLasso.R"))
source(paste0(path_repository,"MultinomLogist_SymLasso.R"))

load(paste0(path_repository,"/DataEx.Rdata"))
# n=1000, p=20, K=7, n_k=(500,200,100,50,50,50,50)
X=Data$X
Y=Data$Y
True_Models=Data$True_Models

True_Models
head(X)
head(Y)

#RES_StdLasso       = MultinomLogist_StdLasso    (X,Y,intercept=TRUE, method=c("BIC","BIC-H"),version_cpp=FALSE)
#RES_SymLasso       = MultinomLogist_SymLasso    (X,Y,intercept=TRUE, method=c("BIC","BIC-H"))

load(paste0(path_repository,"ResEx.Rdata"))

# Comparing the estimated models with the true models

ALL = cbind(True_Models,RES_SymLasso$BIC_H[-1,],RES_StdLasso$BIC_H[-1,])
#ALL = cbind(True_Models,RES_SymLasso$BIC[-1,],RES_StdLasso$BIC[-1,])

#ALL[which(ALL>6)]=6
#ALL[which(ALL< -6)]=-6
RES_DataFrame = data.frame(DELTA=as.numeric(ALL), METH=rep(c("TRUE", "MultinomLogist_SymLasso", "MultinomLogist_StdLasso"), each=20*6), VAR=rep(paste0("Var_",1:20), 3*6), STRAT=  rep(rep(paste0("Cat_",1:6), each=20), 3))
RES_DataFrame$VAR=factor(RES_DataFrame$VAR,levels = paste0("Var_",1:20))
RES_DataFrame$METH=factor(RES_DataFrame$METH,levels = c("TRUE", "MultinomLogist_SymLasso", "MultinomLogist_StdLasso"),labels = c("TRUE","SymLasso","StdLasso"))


library(ggplot2)
ggplot(RES_DataFrame,aes(STRAT,VAR)) + geom_tile(aes(fill = as.numeric(as.character(DELTA))),colour = 'black') + 
  scale_fill_gradient2(high="black",mid="white",low="red",midpoint=0) + 
  facet_wrap(~ METH, ncol=5)+ 
  theme(legend.position="left",legend.title=element_blank(),panel.margin.y = unit(0.8, "lines"),
        strip.background=element_blank(),strip.text=element_text(face="bold",size=1.5),
        axis.ticks = element_line(size = 0.1), axis.ticks.length = unit(0.03, "cm"),
        axis.text=element_text(size=10,face = "bold"), 
        axis.title=element_blank(), axis.text.x=element_text(hjust = -0.4, vjust = 0.4, angle=-90),    
        strip.text.x = element_text(size = 12), legend.text=element_text(size=10,hjust = 1.2,vjust=0.8,face="bold"), legend.key.size = unit(0.5, 'lines'))

