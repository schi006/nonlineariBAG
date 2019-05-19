#' iBag
#' @import MASS glmnet MCMCpack HyperbolicDist mgcv stats ggplot2 d3heatmap
#' @importFrom BART surv.bart
#' @importFrom glmnet cv.glmnet glmnet
NULL


#' Mechanistic model
#' @param meth methylation
#' @param mrna mrna
#' @param cnv copy number
#' @param dsurv survival outcome
#' @export
mechmodel<-function(meth,mrna,cnv,dsurv){
  OurSurvival<-dsurv[,2]
  event<-dsurv[,3]
  names(OurSurvival) <- dsurv[,1]
  OurMRNA<-data.frame(mrna[,2:length(mrna)],row.names = mrna[,1])
  OurMeth<-data.frame(meth[,2:length(meth)],row.names = meth[,1])
  OurCopyNumber<-data.frame(cnv[,2:length(cnv)],row.names = cnv[,1])
  OurGenes<-colnames(OurMRNA)
  cname<- colnames(OurMRNA)

  barcode<-rownames(OurMRNA)
  p<-length(OurGenes)
  n<-dim(OurMRNA)[1]
  k<-3
  num_scores_meth <- rep(NA,p)
  num_scores_CN <- rep(NA,p)
  SST <- rep(NA,p)
  SSM <- rep(NA,p)
  SSCN <- rep(NA,p)
  SSE <- rep(NA,p)
  X <- matrix(NA,nrow=n,ncol=p*k)
  rownames(X) <- barcode
  colnames(X) <- paste(rep(c("Meth","CN","Other"),each=p),rep(OurGenes,3),sep="_")
  for (i in 1:p) {
    ind_meth <- grep(OurGenes[i],colnames(OurMeth))
    if (length(ind_meth)==0) scores_meth <- rep(0,n) else {
      if (length(ind_meth)==1) {
        scores_meth <- as.matrix(OurMeth[,ind_meth])
        num_scores_meth[i] <- 1} else {  ## If only 1 data value, keep raw data (no PCA).
          PCA_meth <- princomp(OurMeth[,ind_meth])
          num_scores_meth[i] <- which(cumsum(PCA_meth$sdev^2/sum(PCA_meth$sdev^2))>=0.9)[1]
          scores_meth  <- PCA_meth$scores[,1:num_scores_meth[i]]
        } }
    # ========
    ind_CN   <- grep(OurGenes[i],colnames(OurCopyNumber))
    if (length(ind_CN)==0) scores_CN <- rep(0,n) else {
      if (length(ind_CN)==1) {
        scores_CN <- OurCopyNumber[,ind_CN]
        num_scores_CN[i] <- 1} else {
          PCA_CN <- princomp(OurCopyNumber[,ind_CN])
          num_scores_CN[i] <- which(cumsum(PCA_CN$sdev^2/sum(PCA_CN$sdev^2))>=0.9)[1]
          scores_CN  <- PCA_CN$scores[,1:num_scores_CN[i]]
        } }

    # ========  	USING GAM INSTEAD OF LEAST SQUARES
    # write formula
    if (length(scores_meth) == n) formula_meth <- "s(scores_meth)" else
    { formula_meth <- paste("s(scores_meth[,",paste(1:num_scores_meth[i], collapse="]) + s(scores_meth[,"),"])",sep='') }
    if (length(scores_CN) == n) formula_CN <- "s(scores_CN)" else
    { formula_CN <- paste("s(scores_CN[,",paste(1:num_scores_CN[i], collapse="]) + s(scores_CN[,"),"])",sep='') }
    formula_all <- paste("OurMRNA[,i] ~ ",formula_meth," + ",formula_CN)
    #
    ## !!! WARNING: I am using a kluge here b/c I know gene 49 is missing methylation.  If want to
    # do this with a new dataset, need to make this more general.

    gam.mRNA  <- gam(as.formula(formula_all))
    # If entire row is 0, coef is NA and scores%*%coef is NA.
    # Estimate pieces
    fit_meth <- as.matrix(predict.gam(gam.mRNA,type="terms")[,1:num_scores_meth[i]] )
    fit_CN <- as.matrix(predict.gam(gam.mRNA,type="terms")[,(num_scores_meth[i]+1):(num_scores_meth[i]+num_scores_CN[i])])

    M <- apply(fit_meth,1,sum)
    CN <- apply(fit_CN,1,sum)
    O <- gam.mRNA$residuals
    X[,paste("Meth",OurGenes[i],sep="_")]   <-  M
    X[,paste("CN",OurGenes[i],sep="_")] <- CN
    X[,paste("Other",OurGenes[i],sep="_")] <- O
    # Pseudo Sums of Squares (to use to find percentages of explained variance)
    SST[i] <- sum( (OurMRNA[,i] - mean(OurMRNA[,i]))^2 )
    SSM[i] <- sum( ( (coef(gam.mRNA)[1] + M) - mean(OurMRNA[,i]) )^2  )
    SSCN[i] <- sum( ( (coef(gam.mRNA)[1] + CN) - mean(OurMRNA[,i]) )^2  )
    SSE[i] <- SST[i] - SSM[i] - SSCN[i]

  }

  tmp<-data.frame(Methylation =SSM/SST*100,CopyNumber =SSCN/SST *100,Other =SSE/SST*100 ,row.names = cname )
  heatmap=d3heatmap(
    tmp,
    # scale(mtcars),
    dendrogram =  "row",
    xaxis_height = 200, yaxis_width = 80,
    xaxis_font_size = 25, yaxis_font_size = 20
  )
  return(list(X=X,
              OurSurvival=OurSurvival,
              event=event,
              num_scores_meth=num_scores_meth,
              num_scores_CN=num_scores_CN,
              OurMRNA=OurMRNA,
              SST=SST, SSM=SSM, SSCN=SSCN, SSE=SSE,heatmap=heatmap))
}


get_starting_values_NG <- function(S, p, k, n, X, Y, names_to_keep, my_seed=sample(99999,size=1)){


  set.seed(my_seed)

  PARAM <- matrix(nrow=S, ncol=sum(p)+1+k+k+sum(p))	# betas, sig_sq, lam, gam_n2, psi[j]s
  beta_names_mat <- matrix(NA,nrow=k,ncol=max(p))	# beta names in matrix form
  for (i in 1:k) {
    for (j in 1:max(p)) { beta_names_mat[i,j] <- paste("beta_",i,".",j,sep="")}
  }
  beta_names <- as.vector(t(beta_names_mat))	# beta names in vector form
  beta_names <- beta_names[ names_to_keep ]	# GET RID OF BETA FOR X COL W/ NO DATA
  psi_names_mat <- matrix(NA,nrow=k,ncol=max(p))
  for (i in 1:k) {
    for (j in 1:max(p)) {psi_names_mat[i,j] <- paste("psi_",i,".",j,sep="")}
  }
  psi_names <- as.vector(t(psi_names_mat))
  psi_names <- psi_names[ names_to_keep ]	# GET RID OF GAM^2 FOR X COL W/ NO DATA
  lam_names <- rep(NA,k)
  for (i in 1:k) { lam_names[i] <- paste("lam_",i,sep="")}
  gam_n2_names <- rep(NA,k)
  for (i in 1:k) { gam_n2_names[i] <- paste("gam_n2_",i,sep="")}
  colnames(PARAM) <- c(beta_names, "sig_sq", psi_names, lam_names, gam_n2_names)


  # ================
  # STARTING VALUES : These should be right in the high density areas of the posteriors
  # ================
  #
  # Get betas from frequentist lasso command glm.net().
  # Get sigma^2 by estimating var(residuals) but divide by n (so, the MLE).
  #

  lasso_out <- glmnet(X,Y)
  cv_lasso_out <- cv.glmnet(X,Y)
  PARAM[1,beta_names] <- coef(lasso_out,s=cv_lasso_out$lambda.min)[-1,]
  # [-1,] b/c don't want the intercept.
  # Using lambda.min instead of lambda.1se b/c lambda.1se returns the largest lambda -- all betas = 0.
  PARAM[1,"sig_sq"] <- sum( (Y-predict(lasso_out,newx=X, s=cv_lasso_out$lambda.min)) ^2 )/(n)
  PARAM[1,psi_names] <- 1
  PARAM[1,lam_names] <- 1
  PARAM[1,gam_n2_names] <- 1



  return( list(PARAM=PARAM,
               beta_names=beta_names,
               gam_n2_names=gam_n2_names,
               psi_names=psi_names,
               lam_names=lam_names,
               my_seed=my_seed) )
}

MC_samples_NG_no_sig_sq_in_beta_prior <- function(PARAM, X, Y, p, k, n, a, b, c, a_tilde, b_tilde, tune, beta_names, gam_n2_names, lam_names,
                                                  psi_names, my_seed=sample(999999,size=1) ){

  set.seed(my_seed) 		# NOTE: important to update initial betas first or run into problems with infinity.
  S <- nrow(PARAM)
  param <- t(as.matrix(PARAM[1,]))


  ## function to use when updating lambda[i]'s. xx is lambda; yy is gamma^(-2)
  prior <- function(xx,yy) {
    (1/xx)^a_tilde * exp(-b_tilde*yy/(2*xx) - c*xx)
  }
  accepted <- rep(0,k)

  for (s in 2:S) {

    ## Update betas ##
    our_mean <- solve(t(X)%*%X+param[1,"sig_sq"]*diag(1/param[1,psi_names]))%*%t(X)%*%Y
    our_cov  <- param[1,"sig_sq"]*solve(t(X)%*%X+param[1,"sig_sq"]*diag(1/param[1,psi_names]))
    param[1,beta_names] <- mvrnorm(n=1, mu=our_mean, Sigma=our_cov)   # n=1 gives one sample of length(mu).

    ## Update sig_sq ##
    our_a <- a+n/2
    our_b <- b+(1/2)*(t(Y-X%*%param[1,beta_names])%*%(Y-X%*%param[1,beta_names]) )
    # NOTE: don't need to use as.matrix for betas vector b/c R coerces it using as.matrix which
    # results in a column vector (which is what we need).
    param[1,"sig_sq"] <- rinvgamma(n=1, shape=our_a, scale=our_b)

    ## Update psi's ##
    for (jj in 1:k) {
      our_aa <- param[1,gam_n2_names[jj]]
      our_bb <- param[1,beta_names[ cumsum(c(1,p))[jj]:cumsum(p)[jj] ]]^2
      our_bb[our_bb<10^(-10)] <- 10^(-10)	# Kluge to keep rgig() happy.
      our_p <- param[1,lam_names[jj]] - (1/2)
      param[1,psi_names[ cumsum(c(1,p))[jj]:cumsum(p)[jj] ]] <-
        apply(as.matrix(our_bb), 1, function(t) rgig(n=1,Theta=c(lambda=our_p,chi=t,psi=our_aa)))
    }

    ## Update lam[i]'s ##
    for (jj in 1:k) {
      lam_star <- exp(tune*rnorm(1))*param[1,lam_names[jj]]
      log_r <- log(prior(lam_star,param[1,gam_n2_names[jj]])) - log(prior(param[1,lam_names[jj]],param[1,gam_n2_names[jj]])) +
        p[jj]*log(gamma(param[1,lam_names[jj]])) - p[jj]*log(gamma(lam_star)) +
        (lam_star-param[1,lam_names[jj]])*( -p[jj]*log(2) + p[jj]*log(param[1,gam_n2_names[jj]]) +
                                              sum(log(param[1,psi_names[cumsum(c(1,p))[jj]:cumsum(p)[jj]]])) )   +
        log(lam_star) - log(param[1,lam_names[jj]])  #This ratio should be here -- confirmed typo in G&B 2010.
      #But we didn't include it in gensips paper. (final significant effects the
      #same either way.)
      acc_prob <- min(log(1),log_r)
      param[1,lam_names[jj]] <- ifelse( log(runif(1))<acc_prob, lam_star, param[1,lam_names[jj]] )
      if (param[1,lam_names[jj]] == lam_star) accepted[jj]<-accepted[jj]+1
    }

    ## Update gam_n2[i]'s ##
    for (jj in 1:k) {
      our_a_tilde     <- p[jj]*param[1,lam_names[jj]] + a_tilde
      our_b_tilde 	<- (1/2)*(b_tilde/param[1,lam_names[jj]] + sum(param[1,psi_names[cumsum(c(1,p))[jj]:cumsum(p)[jj]]]) )
      param[1,gam_n2_names[jj]] <- rgamma(n=1, shape=our_a_tilde, rate=our_b_tilde)
    }

    PARAM[s,] <- param[1,]	# Store updated parameters.
    if (s%%100 == 0) print(paste(s," ",sep=""))	# Keep track of progress because it's SO SLOW.
  }

  return( list(PARAM=PARAM,
               my_seed=my_seed,
               accepted=accepted) )
}


#' BART clinical model
#' @param GBM_data fitted mechanistic model
#' @param nruns nruns
#' @param delta delta
#' @export
clinicalmodel_survbart<-function(GBM_data,nruns,delta){


  burn_in <- floor(0.05*nruns)+1
  n<-dim(GBM_data$OurMRNA)[1]
  p<-dim(GBM_data$OurMRNA)[2]

  oo <- matrix(nrow=3,ncol=p)
  for (i in 1:(p-1)){ oo[,i] <- c(i,i+p,i+2*p+1) }

  m_c=BART::surv.bart(x.train=GBM_data$X, time=GBM_data$OurSurvival,delta = GBM_data$event,x.test=GBM_data$X,ndpost=500)
  rmse=sqrt(mean((m_c$surv.test.mean-as.numeric(GBM_data$event))^2))

  delta_star <- c( log(1-delta), log(1+delta) )
  pos <- apply(m_c$varprob[-(1:burn_in),],2,function(t) mean(t>delta_star[2]))
  neg <- apply(m_c$varprob[-(1:burn_in),],2,function(t) mean(t<delta_star[1]))
  pgene<-which(pos>0.5)
  ngene<-which(neg>0.5)
  cname<- colnames(GBM_data$OurMRNA)
  p<-length(cname)
  pgene1<-pgene[pgene<=p]
  ngene1<-ngene[ngene<=p]
  pgene2<-pgene[(pgene<=2*p) & (pgene>=p+1)]-p
  ngene2<-ngene[(ngene<=2*p) & (ngene>=p+1)]-p
  pgene3<-pgene[(pgene<=3*p) & (pgene>=2*p+1)]-2*p
  ngene3<-ngene[(ngene<=3*p) & (ngene>=2*p+1)]-2*p
  tmp<-array(0,dim=p)
  tmp1<- tmp
  tmp1[pgene1]<-1
  tmp2<- tmp
  tmp2[ pgene2 ]<-1
  tmp3<- tmp
  tmp3[ pgene3]<-1


  t1<-data.frame(meth=tmp1,cnv=tmp2,other=tmp3,row.names = cname)
  maxl<- max(colSums(t1))
  v1<-array("",dim=c(maxl,3))
  c1<-cname[which(t1$meth>0)]
  if(length(c1)>0) v1[1:length(c1),1]<-c1
  c2<-cname[which(t1$cnv>0)]
  if(length(c2)>0) v1[1:length(c2),2]<-c2
  c3<-cname[which(t1$other>0)]
  if(length(c3)>0) v1[1:length(c3),3]<-c3
  colnames(v1)<-c("Methylation","Copy Number", "Other")

  tmp1<- tmp
  tmp1[ngene1]<-1
  tmp2<- tmp
  tmp2[ ngene2 ]<-1
  tmp3<- tmp
  tmp3[ ngene3]<-1


  t2<-data.frame(meth=tmp1,cnv=tmp2,other=tmp3,row.names = cname)
  maxl<- max(colSums(t2))
  v2<-array("",dim=c(maxl,3))
  c1<-cname[which(t2$meth>0)]
  if(length(c1)>0) v2[1:length(c1),1]<-c1
  c2<-cname[which(t2$cnv>0)]
  if(length(c2)>0) v2[1:length(c2),2]<-c2
  c3<-cname[which(t2$other>0)]
  if(length(c3)>0) v2[1:length(c3),3]<-c3
  colnames(v2)<-c("Methylation","Copy Number", "Other")

  pp <- c( sum(grepl("Meth_",colnames(X))), sum(grepl("CN_",colnames(X))),
           sum(grepl("Other_",colnames(X))))		# number of markers per platform
  k <- length(pp)		# number of platforms
  #aa<-data.frame(x=1:length(initial$beta_names),y=pos[oo],g=c(rep(1,to_gibbs$p[1]),rep(2,to_gibbs$p[2]),rep(3,to_gibbs$p[3]))[oo])
  aa<-data.frame(x=1:length(cname),y=pos[oo],g=c(rep(1,pp[1]),rep(2,pp[2]),rep(3,pp[3]))[oo])
  p1<-ggplot()+geom_bar(aes(x=x,y=y,fill = as.factor(g)), position = "dodge", stat="identity",data=aa)+scale_fill_discrete(labels=c("Methylation","Copy Number","Other"))+scale_x_continuous(breaks=c(seq(1,144,by=3)), labels=c(colnames(GBM_data$OurMRNA))  )+theme(axis.text.x = element_text(angle=90))+geom_hline(aes(yintercept=0.5))+theme(legend.title=element_blank())+ylab("Pr(beta > log(1+delta))")+ggtitle("Posterior Probabilities (Positive)")+xlab("Genes")


  aa<-data.frame(x=1:length(cname),y=neg[oo],g=c(rep(1,pp[1]),rep(2,pp[2]),rep(3,pp[3]))[oo])
  #aa<-data.frame(x=1:length(initial$beta_names),y=neg[oo],g=c(rep(1,to_gibbs$p[1]),rep(2,to_gibbs$p[2]),rep(3,to_gibbs$p[3]))[oo])
  p2<-ggplot()+geom_bar(aes(x=x,y=y,fill = as.factor(g)), position = "dodge", stat="identity",data=aa)+scale_fill_discrete(labels=c("Methylation","Copy Number","Other"))+scale_x_continuous(breaks=c(seq(1,144,by=3)), labels=c(colnames(GBM_data$OurMRNA))  )+theme(axis.text.x = element_text(angle=90))+geom_hline(aes(yintercept=0.5))+theme(legend.title=element_blank())+ylab("Pr(beta > log(1-delta))")+ggtitle("Posterior Probabilities (Negative)")+xlab("Genes")


  return(list(GBM_data=GBM_data,oo=oo,p1=p1,p2=p2,t1=v1,t2=v2,rmse=rmse,p=pp,k=k))

}

#' Linear clinical model
#' @param GBM_data fitted mechanistic model
#' @param nruns nruns
#' @param delta delta
#' @param take_log TRUE/FALSE of taking logrithm of clinical outcome
#' @export
clinicalmodel_linear<-function(GBM_data,nruns,delta,take_log){


  burn_in <- floor(0.05*nruns)+1
  n<-dim(GBM_data$OurMRNA)[1]
  p<-dim(GBM_data$OurMRNA)[2]

  oo <- matrix(nrow=3,ncol=p)
  for (i in 1:(p-1)){ oo[,i] <- c(i,i+p,i+2*p+1) }
  names_to_keep <- apply(GBM_data$X,2,function(t) sum(is.na(t))==0)
  X <- apply(GBM_data$X,2,function(t) (t-mean(t))/sd(t) )	# Standardize columns of X.
  if (sum(is.na(X))>0) X <- X[,-which(is.na(X[1,]))]

  if (take_log) {Y <- log(GBM_data$OurSurvival)} else {Y <- GBM_data$OurSurvival}
  Y <- Y-mean(Y)
  initial <- get_starting_values_NG(S=nruns, p=GBM_data$p, k=GBM_data$k, n=nrow(GBM_data$X), X=X, Y=Y,names_to_keep = names_to_keep)
  M <- mean(coef(lm(to_gibbs$Y~to_gibbs$X - 1))^2)	#b_tilde=M per Griffin & Brown (2009)
  pp <- c( sum(grepl("Meth_",colnames(X))), sum(grepl("CN_",colnames(X))),
           sum(grepl("Other_",colnames(X))))		# number of markers per platform
  k <- length(pp)		# number of platforms
  final <- MC_samples_NG_no_sig_sq_in_beta_prior(PARAM=initial$PARAM, X=X, Y=Y, p=pp, k=k, n=nrow(X),a=0.001, b=0.001, c=1, a_tilde=2, b_tilde=M, tune=0.6, beta_names=initial$beta_names,gam_n2_names=initial$gam_n2_names, lam_names=initial$lam_names, psi_names=initial$psi_names)

  post_means <- apply(final$PARAM[(burn_in+1):nrow(final$PARAM),],2,mean)

  delta_star <- c( log(1-delta), log(1+delta) )
  pos <- apply(final$PARAM[-(1:burn_in),initial$beta_names],2,function(t) mean(t>delta_star[2]))
  neg <- apply(final$PARAM[-(1:burn_in),initial$beta_names],2,function(t) mean(t<delta_star[1]))
  pgene<-which(pos>0.5)
  ngene<-which(neg>0.5)
  cname<- colnames(GBM_data$OurMRNA)
  p<-length(cname)
  pgene1<-pgene[pgene<=p]
  ngene1<-ngene[ngene<=p]
  pgene2<-pgene[(pgene<=2*p) & (pgene>=p+1)]-p
  ngene2<-ngene[(ngene<=2*p) & (ngene>=p+1)]-p
  pgene3<-pgene[(pgene<=3*p) & (pgene>=2*p+1)]-2*p
  ngene3<-ngene[(ngene<=3*p) & (ngene>=2*p+1)]-2*p
  tmp<-array(0,dim=p)
  tmp1<- tmp
  tmp1[pgene1]<-1
  tmp2<- tmp
  tmp2[ pgene2 ]<-1
  tmp3<- tmp
  tmp3[ pgene3]<-1


  t1<-data.frame(meth=tmp1,cnv=tmp2,other=tmp3,row.names = cname)
  maxl<- max(colSums(t1))
  v1<-array("",dim=c(maxl,3))
  c1<-cname[which(t1$meth>0)]
  if(length(c1)>0) v1[1:length(c1),1]<-c1
  c2<-cname[which(t1$cnv>0)]
  if(length(c2)>0) v1[1:length(c2),2]<-c2
  c3<-cname[which(t1$other>0)]
  if(length(c3)>0) v1[1:length(c3),3]<-c3
  colnames(v1)<-c("Methylation","Copy Number", "Other")

  tmp1<- tmp
  tmp1[ngene1]<-1
  tmp2<- tmp
  tmp2[ ngene2 ]<-1
  tmp3<- tmp
  tmp3[ ngene3]<-1


  t2<-data.frame(meth=tmp1,cnv=tmp2,other=tmp3,row.names = cname)
  maxl<- max(colSums(t2))
  v2<-array("",dim=c(maxl,3))
  c1<-cname[which(t2$meth>0)]
  if(length(c1)>0) v2[1:length(c1),1]<-c1
  c2<-cname[which(t2$cnv>0)]
  if(length(c2)>0) v2[1:length(c2),2]<-c2
  c3<-cname[which(t2$other>0)]
  if(length(c3)>0) v2[1:length(c3),3]<-c3
  colnames(v2)<-c("Methylation","Copy Number", "Other")

  pp <- c( sum(grepl("Meth_",colnames(X))), sum(grepl("CN_",colnames(X))),
           sum(grepl("Other_",colnames(X))))		# number of markers per platform
  k <- length(pp)		# number of platforms
  aa<-data.frame(x=1:length(initial$beta_names),y=pos[oo],g=c(rep(1,pp[1]),rep(2,pp[2]),rep(3,pp[3]))[oo])
  p1<-ggplot()+geom_bar(aes(x=x,y=y,fill = as.factor(g)), position = "dodge", stat="identity",data=aa)+scale_fill_discrete(labels=c("Methylation","Copy Number","Other"))+scale_x_continuous(breaks=c(seq(1,144,by=3)), labels=c(colnames(GBM_data$OurMRNA))  )+theme(axis.text.x = element_text(angle=90))+geom_hline(aes(yintercept=0.5))+theme(legend.title=element_blank())+ylab("Pr(beta > log(1+delta))")+ggtitle("Posterior Probabilities (Positive)")+xlab("Genes")



  aa<-data.frame(x=1:length(initial$beta_names),y=neg[oo],g=c(rep(1,pp[1]),rep(2,pp[2]),rep(3,pp[3]))[oo])
  p2<-ggplot()+geom_bar(aes(x=x,y=y,fill = as.factor(g)), position = "dodge", stat="identity",data=aa)+scale_fill_discrete(labels=c("Methylation","Copy Number","Other"))+scale_x_continuous(breaks=c(seq(1,144,by=3)), labels=c(colnames(GBM_data$OurMRNA))  )+theme(axis.text.x = element_text(angle=90))+geom_hline(aes(yintercept=0.5))+theme(legend.title=element_blank())+ylab("Pr(beta > log(1-delta))")+ggtitle("Posterior Probabilities (Negative)")+xlab("Genes")


  return(list(GBM_data=GBM_data,oo=oo,p1=p1,p2=p2,t1=v1,t2=v2,p=pp,k=k,post_means=post_means,final=final))

}
