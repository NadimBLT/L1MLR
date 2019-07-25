MultinomLogist_SymLasso <- function(X,Y,intercept=FALSE, method=c("BIC","BIC-H","CV")){
  K=max(Y)
  BIC=BIC_H=CV=NULL
  if("BIC" %in% method || "BIC-H" %in% method){  
    
    if("BIC-H" %in% method){
      RESBH=estimations_glmnet(X,Y,intercept=intercept, hybride=TRUE)
      BIC_H=matrix(RESBH$beta_hybride[, which.min(RESBH$BIC_hybride)],ncol=K)
      BIC_H= BIC_H[,-1]-BIC_H[,1]
      BIC=matrix(RESBH$beta[, which.min(RESBH$BIC)],ncol=K)
      BIC=BIC[,-1]-BIC[,1]
    }else{
      if("BIC" %in% method){
        RESB=estimations_glmnet(X,Y,intercept=intercept, hybride=FALSE)
        BIC=matrix(RESB$beta[, which.min(RESB$BIC)],ncol=K)
        BIC=BIC[,-1]-BIC[,1]
      }
    }
  }
  if("CV" %in% method){
    
    RESCV=cv_estimations_glmnet(X,Y,intercept=intercept, nfolds=10)
    CV=matrix(RESCV,ncol=K)
    CV=CV[,-1]-CV[,1]
  }
  
  return(list(BIC=BIC,BIC_H=BIC_H,CV=CV))
  
}

# Regression logistique multinomiale penalise par le package glmnet
## Entree
# X : matrice de taille n (nombre observations) * p (nombre de variables)
# Y : vecteur numerique ou facteur a K valeurs differentes
# hybride : nouvelle estimation des parametres non nuls sans la penalisation
## Sortie
# renvoit beta sous forme de matrice, avec p*K lignes (ou p+1 si intercept), chaque colonne correspond a un lambda
estimations_glmnet=function(X,Y,family="multinomial",intercept=FALSE,hybride=FALSE,lambda=NULL,...){
  t1=Sys.time()
  n=length(Y)
  K=length(unique(Y))
  niveaux=levels(factor(Y))
  
  modele=glmnet(X,Y,family=family,intercept=intercept,lambda=lambda)
  l=modele$lambda
  L=length(l)
  
  if (intercept){
    beta_glmnet=sapply(1:L,function(i){
      assign("v",eval(parse(text=paste0("c(",
                                        paste0("modele$a0[\"",niveaux,"\",",i,"],","modele$beta[[\"",niveaux,"\"]][,",i,"]",
                                               collapse=","),
                                        ")"))))
      return(v)
    })
    p=ncol(X)+1
  }
  else{
    beta_glmnet=sapply(1:L,function(i){
      assign("v",eval(parse(text=paste0("c(",
                                        paste0("modele$beta[[\"",niveaux,"\"]][,",i,"]",collapse=","),
                                        ")"))))
      return(v)
    })
    p=ncol(X)
  }
  
  dimnames(beta_glmnet)=list(as.vector(outer(paste0("beta",1:p),niveaux,FUN=paste,sep="_")),
                             paste0("lambda",1:L))
  
  dev=deviance(modele)
  nombre_parametre_non_nuls=colSums(modele$dfmat)
  BIC=dev+nombre_parametre_non_nuls*log(n)
  meilleur_lambda=n*l[which.min(BIC)]
  
  beta_glmnet_hybride=NULL
  BIC_hybride=NULL
  if (hybride){
    modele_hybride=lapply(1:L,function(i){
      ensemble=estimation_hybride(X,Y,beta_glmnet[,i],intercept=intercept,mu=TRUE,...)
      ensemble$beta_hybride[1:p]=-ensemble$beta_hybride[1:p]
      return(list(beta_hybride=ensemble$beta_hybride,BIC_hybride=ensemble$BIC_hybride))
    })
    beta_glmnet_hybride=sapply(1:L,function(i){
      modele_hybride[[i]]$beta_hybride
    })
    BIC_hybride=sapply(1:L,function(i){
      modele_hybride[[i]]$BIC_hybride
    })
  }
  t2=Sys.time()
  temps_calcul=as.numeric(t2-t1,units="secs")
  
  return(list(modele=modele,beta=beta_glmnet,beta_hybride=beta_glmnet_hybride,
              BIC=BIC,BIC_hybride=BIC_hybride,meilleur_lambda=meilleur_lambda,temps_calcul=temps_calcul))
}


# Regression logistique multinomiale penalise par le package glmnet par un critere de cross validation
# qui sélectionne le lambda qui donne la moyenne de vraisemblance la plus grande
## Entree
# X : matrice de taille n (nombre observations) * p (nombre de variables)
# Y : vecteur numerique ou facteur a K valeurs differentes
# hybride : nouvelle estimation des parametres non nuls sans la penalisation
## Sortie
# renvoit beta sous forme de vecteur, de taille  p*K (ou p+1 si intercept)
cv_estimations_glmnet=function(X,Y,family="multinomial",intercept=FALSE,nfolds=5,fold_id,random_number,...){
  
  n=length(Y)
  K=length(unique(Y))
  
  if (missing(fold_id)){
    fold_id=cv_partitions(n,nfolds,random_number)
  }
  else{
    if (length(fold_id)!=n){
      stop("fold_id must have the same length as Y")
    }
    nfolds=length(unique(fold_id))
  }
  
  modele_initial=estimations_glmnet(X,Y,family=family,intercept=intercept,...)
  sequence_lambda=modele_initial$modele$lambda
  L=length(sequence_lambda)
  
  matrice_log_vraisemblance=matrix(NA,nrow=L,ncol=nfolds)
  for (fold in 1:nfolds){
    app=which(fold_id!=fold)
    test=which(fold_id==fold)
    modele_apprentissage=estimations_glmnet(X[app,],Y[app],lambda=sequence_lambda,family=family,intercept=intercept,...)
    # Si non convergence pour tous les lambda
    Lapp=ncol(modele_apprentissage$beta)
    for (s in 1:Lapp){
      matrice_log_vraisemblance[s,fold]=calcul_BIC(X[test,],Y[test],modele_apprentissage$beta[,s],intercept=intercept,
                                                   reference=FALSE,mu=FALSE,...)$log_vraisemblance_classique
    }
  }
  rownames(matrice_log_vraisemblance)=colnames(modele_initial$beta)
  moyenne_vraisemblance=rowMeans(matrice_log_vraisemblance,na.rm=TRUE)
  
  indice_lambda_optimal=which.max(moyenne_vraisemblance)
  beta_final=modele_initial$beta[,indice_lambda_optimal]
  return(beta_final)
}

# Fonction qui utilise les beta renvoye des deux algorithmes precedents,
# elle prend une categorie en reference pour la mettre a 0
## Entree
# beta : vecteur beta donne par les deux algorithme precedents
# reference : un niveaux de factor(Y) qui va etre choisit comme reference, si non precise le premier niveau est pris
vecteur_beta_reference=function(beta,reference){
  
  if (is.matrix(beta)){
    
    categorie=sapply(rownames(beta),function(v){
      caracteres_separes=unlist(strsplit(v,split=""))
      indice=which(caracteres_separes=="_")
      return(substr(v,indice+1,nchar(v)))
    })
    niveaux=levels(factor(categorie))
    nombre_categorie=length(niveaux)
    if(missing(reference)){
      reference=niveaux[1]
    }
    indice_colonne=which(niveaux==reference)
    
    beta_ref=apply(beta,2,function(x){
      temp=matrix(x,ncol=nombre_categorie)
      temp=temp-temp[,indice_colonne]
      temp=temp[,-indice_colonne]
      return(as.vector(temp))
    })
    p=nrow(beta)/nombre_categorie
    rownames(beta_ref)=rownames(beta)[-c((p*(indice_colonne-1)+1):(p*indice_colonne))]
  }
  else{
    categorie=sapply(names(beta),function(v){
      caracteres_separes=unlist(strsplit(v,split=""))
      indice=which(caracteres_separes=="_")
      return(substr(v,indice+1,nchar(v)))
    })
    niveaux=levels(factor(categorie))
    nombre_categorie=length(niveaux)
    if(missing(reference)){
      reference=niveaux[1]
    }
    indice_colonne=which(niveaux==reference)
    
    temp=matrix(beta,ncol=nombre_categorie)
    temp=temp-temp[,indice_colonne]
    beta_ref=as.vector(temp[,-indice_colonne])
    p=length(beta)/nombre_categorie
    names(beta_ref)=names(beta)[-c((p*(indice_colonne-1)+1):(p*indice_colonne))]
  }
  return(beta_ref)
}
