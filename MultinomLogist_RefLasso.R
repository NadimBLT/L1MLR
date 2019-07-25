
MultinomLogist_RefLasso <- function(X,Y,intercept=FALSE,version_cpp=FALSE, method=c("BIC","BIC-H","CV")){
  K=max(Y)
  BIC=BIC_H=CV=NULL
  if("BIC" %in% method || "BIC-H" %in% method){  
    
    if("BIC-H" %in% method){
      
      RESBH=estimations_simple_penalisation(X,Y, intercept=intercept,version_cpp=version_cpp,hybride=TRUE)
      BIC_H=matrix(RESBH$beta_hybride[, which.min(RESBH$BIC_hybride)],ncol=K-1)
      BIC=matrix(RESBH$beta[, which.min(RESBH$BIC)],ncol=K-1)
    }else{
      if("BIC" %in% method){
        RESB=estimations_simple_penalisation(X,Y, intercept=intercept,version_cpp=version_cpp,hybride=FALSE)
        BIC=matrix(RESB$beta[, which.min(RESB$BIC)],ncol=K-1)
      }
    }
  }
  if("CV" %in% method){
    
    RESCV=cv_estimations_simple_penalisation(X,Y,intercept=intercept,version_cpp=version_cpp, nfolds=10)
    CV=matrix(RESCV,ncol=K-1)
  }
  
  return(list(BIC=BIC,BIC_H=BIC_H,CV=CV))
  
}




###################################################################################################
# Algorithmes avec beta sous forme de vecteur
###################################################################################################


# Algorithme qui effectue le cycle (a...b) du groupe symetrique Sn, necessaire pour la fonction suivante
cycle=function(sigma,a,b){
  N=length(sigma)
  
  if (a>b){
    stop("a must be inferior or equal to b")
  }
  
  if (a==b){
    s=1:N
  }
  
  else if ((1<a) & (b<N)){
    s=c(1:(a-1),b,a:(b-1),(b+1):N)
  }
  else if ((a==1) & (b<N)){
    s=c(b,1:(b-1),(b+1):N)
  }
  else if ((1<a) & (b==N)){
    s=c(1:(a-1),N,a:(N-1))
  }
  else{
    s=c(N,1:(N-1))
  }
  sigma=sigma[s]
  return(sigma)
}

# Modification de la structure de la matrice X pour realiser la regression logistique multinomiale en version vecteur, permet egalement
# de forcer des egalites de parametres
## Entree
# X : matrice de taille n (nombre observations) * p (nombre variables)
# nombre_classe : nombre de classe pour la regression logistique, avec ou sans la reference
# mu : si vrai, fait la decomposition adhequate pour faire la double penalisation
# penalise : enleve les colonnes specifie de la matrice finale de nombre de colonnes p*nombre_classe, afin de garantir d'avoir
# 0 aux composantes du vecteur beta enlevees pour la version hybride de l'algorithme general
matrice_bloc_multinomial_logistic=function(X,nombre_classe,mu=FALSE,tau,meme_influence=NULL,
                                           penalise=NULL,...){
  
  
  X0=vector("list",nombre_classe)
  X=as.matrix(X)
  p=ncol(X)
  
  if (mu){
    if (missing(tau)){
      tau=rep(1,nombre_classe)
    }
    else if (length(tau)!=(nombre_classe)){
      tau=rep(tau,length.out=nombre_classe)
    }
    
    X0=bdiag(lapply(1:nombre_classe,FUN=function(i){DX=X/tau[i]
    return(DX)}))
    L=do.call(rbind,lapply(1:nombre_classe,FUN=function(i){LX=X
    return(LX)}))
    X0=cbind(L,X0)
    # Dans ce cas de figure :
    # ncol(X0)=(nombre_classe+1)*p
  }
  
  else if (is.null(meme_influence)){
    X0=bdiag(lapply(1:nombre_classe,FUN=function(i){DX=X
    return(DX)}))
    # Dans ce cas de figure :
    # ncol(X0)=nombre_classe*p
  }
  
  else{
    L=NULL
    Xs=X[,-meme_influence]
    Xa=matrix(X[,meme_influence],ncol=length(meme_influence))
    
    X0=bdiag(lapply(1:nombre_classe,FUN=function(i){DX=Xs
    return(DX)}))
    L=do.call(rbind,lapply(1:nombre_classe,FUN=function(i){LX=Xa
    return(LX)}))
    X0=cbind(L,X0)
    sigma=1:p
    m=0
    for (i in sort(meme_influence)){
      m=m+1
      sigma=cycle(sigma,m,i)
    }
    ordre=order(sigma)
    X0[,1:p]=X0[,ordre]
    # Dans ce cas de figure :
    # ncol(X0)=nombre_classe*(p-length(meme_influence))+length(meme_influence)
  }
  
  if (is.null(penalise)==FALSE){
    X0=X0[,-penalise]
  }
  
  return(X0)
}

## Entree
# X : matrice de taille n (nombre observations) * p (nombre de variables)
# Y : vecteur numerique ou facteur a K valeurs differentes. Si reference=TRUE, la reference est choisie comme le premier niveau (c'est a dire levels(factor(Y))[1])
# lambda : vecteur lambda de penalisation a specifie, mis a 0.1 par defaut
# beta_initial : beta_initial pour l'algorithme iteratif
# reference : renvoyer beta sous K composantes si FALSE, K-1 si TRUE
# mu : si TRUE, change la structure de la matrice X pour effectuer la double penalisation, est utilise pour la fonction suivante
# tolerance : seul d'arret de l'algorithme iteratif pour la difference de beta_(m+1)-beta_m
# max_iter : nombre maximum d'iterations
# seed : graine pour le generateur de nombre aleatoire
# version_cpp : si TRUE, effectue l'algorithme iteratif avec une fonction faisant appel au langage Cpp. Cela rend la fonction plus rapide
## Sortie
# renvoit une liste donnant beta sous forme de vecteur de taille p0*K0, avec P0=p ou p+1 si intercept, K0=K ou K-1 si reference
# lambda_max donne le lambda maximal qu'on peut mettre dans l'entree pour avoir beta non nul (sans compter l'intercept)
multinomial_logistic_lasso_vecteur = function(X,
                                              Y,
                                              lambda,
                                              beta_initial,
                                              intercept = FALSE,
                                              reference=FALSE,
                                              mu=FALSE,
                                              meme_influence=NULL,
                                              tolerance = 1e-3,
                                              max_iter = 1e5,
                                              decay_rate = 0.25,
                                              decay_min = 0,
                                              seed = 1,
                                              version_cpp=FALSE,
                                              ...) {
  # Dimensions des objets
  if (intercept){
    X=cbind(rep(1,nrow(X)),X)
  }
  p=ncol(X)
  
  if (!is.factor(Y)) {
    Y = factor(Y)
  }
  K = nlevels(Y)
  n = length(Y)
  
  if (reference==FALSE){
    nom = levels(Y)
    epsilon=0
    X0=matrice_bloc_multinomial_logistic(X,nombre_classe = K,mu=mu,meme_influence=meme_influence,...)
  }
  else{
    nom=levels(Y)[-1]
    epsilon=1
    X0=matrice_bloc_multinomial_logistic(X,nombre_classe = K-1,mu=mu,meme_influence=meme_influence,...)
  }
  X0=as.matrix(X0) # a mettre ou pas, cela depend de la taille des donnees
  pp=ncol(X0)
  L=nrow(X0) # L est forcement un multiple de n
  m=L/n # vaut K-1 avec reference=TRUE, K sinon
  
  Y0=c(sapply(nom,FUN=function(i){
    a=rep(0,n)
    a[Y==i]=1
    return(a)
  }))
  
  if (missing(lambda)) {
    lambda = rep(0.1,pp)
  }
  else if (length(lambda) != pp) {
    lambda = rep(lambda, length.out = pp)
  }
  
  if (intercept){
    if (mu==FALSE){
      if (is.null(meme_influence)){
        indice_intercept=seq(1,pp,by=p)
        lambda[indice_intercept]=0
      }
      else if (1 %in% meme_influence){
        lambda[1]=0
      }
      else{
        indice_intercept=c(1,seq(p+1,pp,by=p-length(meme_influence)))
        lambda[indice_intercept]=0
      }
    }
    else{
      lambda[1]=0
    }
  }
  
  # Valeur initiale
  if (missing(beta_initial)){
    beta=rep(0, pp)
  }
  else{
    beta=rep(as.vector(beta_initial),length.out=pp)
  }
  
  # precomputation de certaines valeurs
  B = ((K - 1) / (2 * K)) * colSums(X0 ^ 2)
  softmax_delta = lambda / B
  XY=as.vector(Y0 %*% X0)
  NN=lapply(seq_len(pp),FUN=function(j){which(X0[,j]!=0)})
  
  # Probabilites d'appartenance a une classe
  Xbeta = as.vector(X0 %*% beta)
  E = exp(Xbeta)
  JJ=do.call(cbind,lapply(1:m,FUN=function(i){U=Diagonal(n)
  return(U)}))
  S=epsilon+as.vector(JJ %*% E)
  S0 = rep(S, m)
  P=E/S0
  
  beta_resamp = rep(1, pp)
  
  # Indicateurs de performance
  incr = .Machine$double.xmax
  
  set.seed(seed)
  
  if (version_cpp){
    modele_cpp=multinomial_logistic_lasso_vecteur_optimal(beta,
                                                          X0,
                                                          Y0,
                                                          pp,
                                                          L,
                                                          n,
                                                          P,
                                                          beta_resamp,
                                                          XY,
                                                          Xbeta,
                                                          E,
                                                          S,
                                                          S0,
                                                          B,
                                                          NN,
                                                          lambda,
                                                          epsilon,
                                                          max_iter,
                                                          reference,
                                                          tolerance,
                                                          decay_min,
                                                          decay_rate)
    
    beta=modele_cpp$beta
    Xbeta=modele_cpp$Xbeta
    S=modele_cpp$S
  }
  
  else{
    # Algorithme
    for (iter in 1:max_iter) {
      beta_prev = beta
      for (l in 1:pp) {
        beta_old = beta[l]
        
        if ((beta_old != 0) | (runif(1) <= beta_resamp[l])) {
          # grad = XY[l] - X0[,l] %*% P
          grad = XY[l] - sum(X0[,l]*P)
          beta_new = beta_old + grad / (B[l])
          
          # Calcul
          beta_new = sign(beta_new) * max(0, abs(beta_new) - softmax_delta[l])
          
          # Affectation si changement
          if (beta_new != beta_old) {
            Xbeta[NN[[l]]] = Xbeta[NN[[l]]] + X0[NN[[l]],l] * (beta_new - beta_old)
            # E_new_classe = exp(Xbeta[NN[[l]]])
            # E[NN[[l]]] = E_new_classe
            E[NN[[l]]]=exp(Xbeta[NN[[l]]])
            S=epsilon+(JJ %*% E)
            S0 = rep(as.vector(S), m)
            P=E/S0
            
            beta[l] = beta_new
          }
          
          # Changement des probabilites
          if ((beta_new != 0) & (beta_old == 0)) {
            beta_resamp[l] = 1
          }
          else if ((beta_new == 0) & (beta_old == 0)) {
            beta_resamp[l] = (beta_resamp[l] - decay_min) * decay_rate + decay_min
          }
        }
      }
      incr = norm(beta_prev - beta, "2") / (norm(beta_prev,"2") + .Machine$double.eps)
      
      if (incr < tolerance) {
        break
      }
    }
  }
  
  # log_vraisemblance=sum(Xbeta*Y0)-sum(log(S))-sum(lambda*abs(beta))
  log_vraisemblance_classique=sum(Xbeta*Y0)-sum(log(S))
  # nb_para_libre=length(unique(beta[beta!=0]))
  nb_para_libre=length(beta[beta!=0])
  BIC=-2*log_vraisemblance_classique+nb_para_libre*log(n)
  lambda_max=max(abs(XY-1/K*colSums(X0)))
  return(list(beta=beta,log_vraisemblance_classique=log_vraisemblance_classique,BIC=BIC,lambda_max=lambda_max))
}



# Fonction qui effectue la decomposition beta_j=mu+gamma_j, j=2..K
# Entree
# X Y et lambda pareil que fonction precedente
# tau sert a penaliser gamma, lambda_2=tau*lambda_1
mu_gamma_vecteur=function(X,Y,lambda,tau,intercept = FALSE,...){
  if (!is.factor(Y)) {
    Y = factor(Y)
  }
  K=nlevels(Y)
  noms=levels(Y)[-1]
  if (missing(tau)){
    tau=rep(1,K-1)
  }
  else if (length(tau)!=(K-1)){
    tau=rep(tau,length.out=K-1)
  }
  
  modele=multinomial_logistic_lasso_vecteur(X=X,Y=Y,lambda = lambda,reference=TRUE,mu=TRUE,tau=tau,intercept=intercept,...)
  beta_vecteur=modele$beta
  if (intercept){
    X=cbind(rep(1,nrow(X)),X)
  }
  normalise=matrix(rep(c(1,tau),each=ncol(X)),ncol=K)
  Omega=matrix(beta_vecteur/normalise,ncol=K,dimnames=list(colnames(X),c("mu",paste0("gamma_",noms))))
  beta=matrix(0,ncol=K-1,nrow=nrow(Omega),dimnames=list(colnames(X),paste0("beta_",noms)))
  for (i in 1:(K-1)){
    beta[,i]=Omega[,1]+Omega[,i+1]
  }
  
  # Mode <- function(x) {
  #   ux <- unique(x)
  #   ux[which.max(tabulate(match(x, ux)))]
  # }
  # 
  # M=sapply(1:ncol(X),FUN=function(i){
  #   return(Mode(beta[i,]))
  # })
  # 
  # M_beta=cbind(-M,beta-M)
  # nb_para_libre=length(M_beta[M_beta!=0])
  # BIC_M=-2*modele$log_vraisemblance_classique+nb_para_libre*log(length(Y))
  
  return(list(mu_gamma=Omega,beta=beta,modele=modele))
}


# Effectue algorithme mu_gamma_vecteur sur un ensemble de lambda et tau
## Entree
# ratio_lambda : rapport entre lambda_min et lambda_max
# sequence_lambda : si non specifie, l'algorithme choisit lui meme la sequence de lambda, sinon il utilise celle qui est donnee
# taille_obs : si TRUE, tau_j est donne comme sqrt(sum(Y==j))
# hybride : si TRUE, effectue l'estimation hybride, voir fonction estimation hybride
estimations_double_penalisation=function(X,Y,nombre_lambda=50,nombre_tau=50,
                                         ratio_lambda=0.001,ratio_tau=0.001,
                                         sequence_lambda=NULL,sequence_tau=NULL,
                                         lambda_max,
                                         taille_obs=FALSE,
                                         intercept=FALSE,hybride=FALSE,...){
  t1=Sys.time()
  if (!is.factor(Y)) {
    Y = factor(Y)
  }
  
  categories=levels(Y)
  K=nlevels(Y)
  
  
  if (taille_obs==FALSE){
    a=rep(1,K-1)
  }
  else{
    a=sqrt(as.vector(table(Y)))[-1]
  }
  
  # base=multinomial_logistic_lasso_vecteur(X,Y,lambda=0,intercept=intercept,...)
  base=mu_gamma_vecteur(X,Y,lambda=0,tau=a,intercept = intercept,...)
  if (missing(lambda_max)){
    C=base$modele$lambda_max
  }
  else{
    C=lambda_max
  }
  p=length(base$modele$beta)/K
  
  if (is.null(sequence_lambda)){
    lambda=exp(seq(from=log(ratio_lambda*C),to=log(C),length.out=nombre_lambda))
  }
  else{
    lambda=sequence_lambda
    nombre_lambda=length(lambda)
  }
  if(is.null(sequence_tau)==FALSE){
    tau=sequence_tau
    nombre_tau=length(tau)
  }
  
  giga_beta=matrix(0,nrow=p*(K-1),ncol=nombre_lambda*nombre_tau,
                                     dimnames=list(paste0(paste0("beta_",rep(categories[-1],each=p)),"_",rep(1:p,K-1)),
                                                   paste(rep(paste0("lambda_",1:nombre_lambda),each=nombre_tau),
                                                         rep(paste0("tau_",1:nombre_tau),nombre_lambda),sep="_")))
  
  
  giga_mu_gamma=matrix(0,nrow=p*K,ncol=nombre_lambda*nombre_tau,
                                             dimnames=list(c(paste0("mu_",1:p),paste0(paste0("gamma_",rep(categories[-1],each=p)),"_",rep(1:p,K-1))),
                                                           paste(rep(paste0("lambda_",1:nombre_lambda),each=nombre_tau),
                                                                 rep(paste0("tau_",1:nombre_tau),nombre_lambda),sep="_")))
  
  giga_BIC=matrice_tau=matrix(0,nrow=nombre_lambda,ncol=nombre_tau,dimnames=list(paste0("lambda_",1:nombre_lambda),paste0("tau_",1:nombre_tau)))
  # beta_initial=0
  
  if (hybride){
    giga_beta_hybride=giga_beta
    giga_mu_gamma_hybride=giga_mu_gamma
    giga_BIC_hybride=giga_BIC
  }
  
  else{
    giga_beta_hybride=giga_mu_gamma_hybride=giga_BIC_hybride=NULL
  }
  
  
  for (i in 1:nombre_lambda){
    # tau=seq(from=0.5,to=C/lambda[i],length.out=nombre_tau)
    if (is.null(sequence_tau)){
      tau=exp(seq(from=log(ratio_tau*C/(lambda[i])),to=log(C/lambda[i]),length.out=nombre_tau))
    }
    matrice_tau[i,]=tau
    beta_initial=0
    
    for (j in 1:nombre_tau){
      ensemble=mu_gamma_vecteur(X,Y,lambda=lambda[i],tau=a*tau[j],intercept=intercept,beta_initial=beta_initial,...)
      giga_beta[,nombre_tau*(i-1)+j]=as.vector(ensemble$beta)
      giga_mu_gamma[,nombre_tau*(i-1)+j]=as.vector(ensemble$mu_gamma)
      giga_BIC[i,j]=ensemble$modele$BIC
      beta_initial=as.vector(ensemble$modele$beta)
      
      if (hybride){
        ensemble_hybride=estimation_hybride(X,Y,as.vector(ensemble$mu_gamma),mu=TRUE,
                                            intercept=intercept,
                                            beta_initial=beta_initial,...)
        giga_mu_gamma_hybride[,nombre_tau*(i-1)+j]=ensemble_hybride$beta_hybride
        giga_beta_hybride[,nombre_tau*(i-1)+j]=giga_mu_gamma_hybride[,nombre_tau*(i-1)+j][-(1:p)]+
          rep(giga_mu_gamma_hybride[,nombre_tau*(i-1)+j][1:p],K-1)
        giga_BIC_hybride[i,j]=ensemble_hybride$BIC_hybride
      }
    }
  }
  meilleur_BIC=vector("list",2)
  names(meilleur_BIC)=c("indice_matrice_BIC", "indice_matrice_beta")
  meilleur_BIC[["indice_matrice_BIC"]]=which(giga_BIC==min(giga_BIC),arr.ind=TRUE)
  meilleur_BIC[["indice_matrice_beta"]]=sapply(nrow(meilleur_BIC[["indice_matrice_BIC"]]),function(l){
    return((meilleur_BIC[["indice_matrice_BIC"]][l,1]-1)*nombre_tau+meilleur_BIC[["indice_matrice_BIC"]][l,2])
  })
  # meilleur_BIC_hybride=which(giga_BIC_hybride==min(giga_BIC_hybride),arr.ind=TRUE)
  meilleur_lambda_tau=NULL
  for (l in 1:nrow(meilleur_BIC[["indice_matrice_BIC"]])){
    meilleur_lambda_tau=rbind(meilleur_lambda_tau,c(lambda[meilleur_BIC[["indice_matrice_BIC"]][l,1]],
                                                    matrice_tau[meilleur_BIC[["indice_matrice_BIC"]][l,1],
                                                                meilleur_BIC[["indice_matrice_BIC"]][l,2]]))
  }
  t2=Sys.time()
  temps_calcul=as.numeric(t2-t1,units="secs")
  
  return(list(beta=giga_beta,beta_hybride=giga_beta_hybride,
              mu_gamma=giga_mu_gamma,mu_gamma_hybride=giga_mu_gamma_hybride,
              BIC=giga_BIC,BIC_hybride=giga_BIC_hybride,
              meilleur_BIC=meilleur_BIC,lambda=lambda,tau=matrice_tau,
              meilleur_lambda_tau=meilleur_lambda_tau,temps_calcul=temps_calcul))
}


# Fonction similaire a la precedent mais avec un seul lambda
estimations_simple_penalisation=function(X,Y,nombre_lambda=50,ratio_lambda=0.001,sequence_lambda=NULL,
                                         reference=TRUE,intercept=FALSE,hybride=FALSE,...){
  t1=Sys.time()
  if (!is.factor(Y)) {
    Y = factor(Y)
  }
  K=nlevels(Y)
  if (reference){
    J=K-1
    noms=levels(Y)[-1]
  }
  else{
    J=K
    noms=levels(Y)
  }
  base=multinomial_logistic_lasso_vecteur(X,Y,lambda=0,reference=reference,intercept=intercept,...)
  C=base$lambda_max
  p=length(base$beta)/J
  
  if(is.null(sequence_lambda)){
    lambda=exp(seq(from=log(ratio_lambda*C),to=log(C),length.out=nombre_lambda))
  }
  else{
    lambda=sequence_lambda
    nombre_lambda=length(sequence_lambda)
  }
  
  giga_beta=matrix(NA,nrow=p*J,ncol=nombre_lambda,dimnames=list(paste0(paste0("beta_",rep(noms,each=p)),"_",rep(1:p,J)),
                                                               paste0("lambda_",1:nombre_lambda)))
  
  
  liste_BIC=rep(0,nombre_lambda)
  names(liste_BIC)=paste0("lambda_",1:nombre_lambda)
  beta_initial=0
  
  if (hybride){
    giga_beta_hybride=giga_beta
    liste_BIC_hybride=liste_BIC
  }
  else{
    giga_beta_hybride=liste_BIC_hybride=NULL
  }
  
  for (i in 1:nombre_lambda){
    ensemble=multinomial_logistic_lasso_vecteur(X,Y,lambda=lambda[i],mu=FALSE,reference=reference,intercept=intercept,beta_initial=beta_initial,...)
    giga_beta[,i]=ensemble$beta
    liste_BIC[i]=ensemble$BIC
    beta_initial=ensemble$beta
    
    if (hybride){
      ensemble_hybride=estimation_hybride(X,Y,ensemble$beta,intercept=intercept,reference=reference,mu=FALSE,...)
      giga_beta_hybride[,i]=ensemble_hybride$beta_hybride
      liste_BIC_hybride[i]=ensemble_hybride$BIC_hybride
    }
  }
  meilleur_BIC=which(liste_BIC==min(liste_BIC))
  meilleur_lambda=rep(NA,length(meilleur_BIC))
  for (l in 1:length(meilleur_BIC)){
    meilleur_lambda[l]=lambda[meilleur_BIC[l]]
  }
  t2=Sys.time()
  temps_calcul=as.numeric(t2-t1,units="secs")
  
  return(list(beta=giga_beta,beta_hybride=giga_beta_hybride,
              BIC=liste_BIC,BIC_hybride=liste_BIC_hybride,
              meilleur_BIC=meilleur_BIC,lambda=lambda,meilleur_lambda=meilleur_lambda,temps_calcul=temps_calcul))
}

# reeffectue fonction multinomial_logistic_lasso_vecteur sur les composantes non nulles de beta pour les reestimer sans la penalisation lambda
estimation_hybride=function(X,Y,beta,intercept=FALSE,mu=FALSE,reference=TRUE,...){
  
  L=length(beta)
  beta_hybride=rep(0,L)
  names(beta_hybride)=names(beta)
  K=length(unique(Y))
  BIC_hybride=2*length(Y)*log(K) # si tous les beta=0, BIC=2*n*log(K)
  
  
  if (mu==TRUE){
    beta_zero=NULL
    p=L/K
    Omega=matrix(beta,ncol=K,nrow=p)
    for (j in 1:p){
      sequence_indice=seq(from=j,to=p*K,by=p)
      if (all(Omega[j,]==0)){
        beta_zero=c(beta_zero,sequence_indice)
      }
      else if (prod(Omega[j,-1])!=0){
        beta_zero=c(beta_zero,j)
      }
      # else if (Omega[j,1]==0){
      #   indice_a_enlever=which(Omega[j,sequence_indice]==0)
      #   beta_zero=c(beta_zero,sequence_indice[indice_a_enlever])
      # }
      else {
        indice_a_enlever=which(Omega[j,]==0)
        beta_zero=c(beta_zero,sequence_indice[indice_a_enlever])
      }
    }
    beta_zero=sort(beta_zero)
    
    
    if (any(Omega!=0)){
      modele_hybride=multinomial_logistic_lasso_vecteur(X,Y,reference = TRUE,intercept=intercept,
                                                        lambda=0,mu=TRUE,tau=1,
                                                        penalise=beta_zero,...)
      beta_non_nul=modele_hybride$beta
      
      beta_hybride[setdiff(1:(p*K),beta_zero)]=beta_non_nul
      BIC_hybride=modele_hybride$BIC
    }
    
  }
  
  else{
    beta_zero=which(beta==0)
    if (length(beta_zero)==0){
      beta_zero=NULL
    }
    
    if (length(beta_zero)!=L){
      ensemble_hybride=multinomial_logistic_lasso_vecteur(X,Y,lambda=0,reference=reference,
                                                          intercept=intercept,penalise=beta_zero,...)
      beta_hybride[setdiff(1:L,beta_zero)]=ensemble_hybride$beta
      BIC_hybride=ensemble_hybride$BIC
    }
  }
  
  
  return(list(beta_hybride=beta_hybride,BIC_hybride=BIC_hybride))
}




calcul_BIC=function(X,Y,beta,intercept=FALSE,reference=FALSE,mu=FALSE,meme_influence=NULL,...){
  
  if (intercept){
    X=cbind(rep(1,nrow(X)),X)
  }
  p=ncol(X)
  
  if (!is.factor(Y)) {
    Y = factor(Y)
  }
  K = nlevels(Y)
  n = length(Y)
  
  if (reference==FALSE){
    nom = levels(Y)
    epsilon=0
    X0=matrice_bloc_multinomial_logistic(X,nombre_classe = K,mu=mu,meme_influence=meme_influence,...)
  }
  else{
    nom=levels(Y)[-1]
    epsilon=1
    X0=matrice_bloc_multinomial_logistic(X,nombre_classe = K-1,mu=mu,meme_influence=meme_influence,...)
  }
  X0=as.matrix(X0) # a mettre ou pas, cela depend de la taille des donnees
  pp=ncol(X0)
  L=nrow(X0) # L est forcement un multiple de n
  m=L/n # vaut K-1 avec reference=TRUE, K sinon
  
  Y0=c(sapply(nom,FUN=function(i){
    a=rep(0,n)
    a[Y==i]=1
    return(a)
  }))
  
  
  Xbeta = as.vector(X0 %*% beta)
  E = exp(Xbeta)
  JJ=do.call(cbind,lapply(1:m,FUN=function(i){U=Diagonal(n)
  return(U)}))
  S=epsilon+as.vector(JJ %*% E)
  
  log_vraisemblance_classique=sum(Xbeta*Y0)-sum(log(S))
  # nb_para_libre=length(unique(beta[beta!=0]))
  nb_para_libre=length(beta[beta!=0])
  BIC=-2*log_vraisemblance_classique+nb_para_libre*log(n)
  
  return(list(log_vraisemblance_classique=log_vraisemblance_classique,nb_para_libre=nb_para_libre,BIC=BIC))
}


# partition echantillon pour realiser cv
## Entree
# n : taille echantillon
# nfolds : nombre de sous echantillons
# random_number : graine (seed) optionnelle pour le generateur de nombre aleatoires
## Sortie
# vecteur de taille n avec des valeurs de 1 a nfolds
cv_partitions=function(n,nfolds,random_number){
  if (!missing(random_number)){
    set.seed(random_number)
  }
  tirage=sample(n,replace=FALSE)
  
  q=floor(n/nfolds)
  fold_id=rep(NA,n)
  
  if (q*nfolds<n){
    for (i in 1:(nfolds-1)){
      fold_id[tirage[((i-1)*q+1):(i*q)]]=i
    }
    fold_id[tirage[((nfolds-1)*q+1):n]]=nfolds
  }
  
  else {
    for (i in 1:nfolds){
      fold_id[tirage[((i-1)*q+1):(i*q)]]=i
    }
  }
  return(fold_id)
}


# Fonction qui choisit le lambda et tau optimal par un critere de cross validation 
# Entree 
# fold_id : mettre le vecteur fold_id de la fonction precedente, si non precise, la fonction va le creer
# a l'aide de la fonction cv_partitions precedente avec le nfolds specifie
## Sortie
# renvoit le vecteur Omega et beta du modele optimal choisit par cv
cv_estimations_double_penalisation=function(X,
                                            Y,
                                            nombre_lambda=50,
                                            nombre_tau=50,
                                            nfolds=5,
                                            fold_id,
                                            ratio_lambda=0.001,
                                            ratio_tau=0.001,
                                            sequence_lambda=NULL,
                                            sequence_tau=NULL,
                                            lambda_max,
                                            intercept=FALSE,
                                            random_number,
                                            ...) {
  n=length(Y)
  base=multinomial_logistic_lasso_vecteur(X,Y,lambda=0,intercept=intercept,...)
  if (missing(lambda_max)){
    C=base$lambda_max
  }
  else{
    C=lambda_max
  }
  
  if (missing(fold_id)){
    fold_id=cv_partitions(n,nfolds,random_number)
  }
  else{
    if (length(fold_id)!=n){
      stop("fold_id must have the same length as Y")
    }
    nfolds=length(unique(fold_id))
  }
  
  
  if (is.null(sequence_lambda)){
    sequence_lambda=exp(seq(from=log(ratio_lambda*C),to=log(C),length.out=nombre_lambda))
  }
  else{
    nombre_lambda=length(sequence_lambda)
  }
  
  if (!is.null(sequence_tau)){
    nombre_tau=length(sequence_tau)
  }
  
  matrice_log_vraisemblance=matrix(NA,nrow=nombre_lambda*nombre_tau,ncol=nfolds)
  for (fold in 1:nfolds){
    app=which(fold_id!=fold)
    test=which(fold_id==fold)
    modele_apprentissage=estimations_double_penalisation(X[app,],Y[app],sequence_lambda=sequence_lambda,
                                                         sequence_tau=sequence_tau,nombre_tau=nombre_tau,
                                                         lambda_max=C,intercept=intercept,...)
    
    for (s in 1:(nombre_lambda*nombre_tau)){
      matrice_log_vraisemblance[s,fold]=calcul_BIC(X[test,],Y[test],
                                                   modele_apprentissage$mu_gamma[,s],
                                                   intercept=intercept,
                                                   reference=TRUE,mu=TRUE,...)$log_vraisemblance_classique
    }
  }
  rownames(matrice_log_vraisemblance)=colnames(modele_apprentissage$mu_gamma)
  moyenne_vraisemblance=rowMeans(matrice_log_vraisemblance)
  meilleur_lambda_tau=which.max(moyenne_vraisemblance)
  
  indice_lambda_optimal=floor((meilleur_lambda_tau-1)/nombre_tau)+1
  lambda_optimal=sequence_lambda[indice_lambda_optimal]
  sequence_tau=modele_apprentissage$tau[indice_lambda_optimal,]
  indice_tau_optimal=((meilleur_lambda_tau)%%(nombre_tau))+nombre_tau*(((meilleur_lambda_tau)%%(nombre_tau))==0)
  tau_optimal=sequence_tau[indice_tau_optimal]
  modele_final=mu_gamma_vecteur(X,Y,intercept=intercept,lambda=lambda_optimal,tau=tau_optimal,...)
  return(list(mu_gamma=modele_final$mu_gamma,beta=modele_final$beta))
}


# Fonction similaire a la precedente avec juste un seul lambda qui penalise
cv_estimations_simple_penalisation=function(X,
                                            Y,
                                            nombre_lambda=50,
                                            nfolds=5,
                                            fold_id,
                                            ratio_lambda=0.001,
                                            sequence_lambda=NULL,
                                            intercept=FALSE,
                                            reference=TRUE,
                                            random_number,
                                            ...) {
  n=length(Y)
  base=multinomial_logistic_lasso_vecteur(X,Y,lambda=0,intercept=intercept,...)
  C=base$lambda_max
  
  if (missing(fold_id)){
    fold_id=cv_partitions(n,nfolds,random_number)
  }
  else{
    if (length(fold_id)!=n){
      stop("fold_id must have the same length as Y")
    }
    nfolds=length(unique(fold_id))
  }
  
  
  if (is.null(sequence_lambda)){
    sequence_lambda=exp(seq(from=log(ratio_lambda*C),to=log(C),length.out=nombre_lambda))
  }
  else{
    nombre_lambda=length(sequence_lambda)
  }
  
  matrice_log_vraisemblance=matrix(NA,nrow=nombre_lambda,ncol=nfolds)
  for (fold in 1:nfolds){
    app=which(fold_id!=fold)
    test=which(fold_id==fold)
    modele_apprentissage=estimations_simple_penalisation(X[app,],Y[app],sequence_lambda=sequence_lambda,
                                                         ratio_lambda=ratio_lambda,
                                                         intercept=intercept,reference=reference,...)
    
    for (s in 1:nombre_lambda){
      matrice_log_vraisemblance[s,fold]=calcul_BIC(X[test,],Y[test],
                                                   modele_apprentissage$beta[,s],
                                                   intercept=intercept,
                                                   reference=reference,mu=FALSE,...)$log_vraisemblance_classique
    }
  }
  rownames(matrice_log_vraisemblance)=colnames(modele_apprentissage$beta)
  moyenne_vraisemblance=rowMeans(matrice_log_vraisemblance)
  
  indice_lambda_optimal=which.max(moyenne_vraisemblance)
  lambda_optimal=sequence_lambda[indice_lambda_optimal]
  modele_final=multinomial_logistic_lasso_vecteur(X,Y,intercept=intercept,lambda=lambda_optimal,reference=reference,...)
  return(modele_final$beta)
}
