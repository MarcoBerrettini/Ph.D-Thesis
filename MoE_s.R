#kn <- function(x, xl, xr, ndx, bdeg) {
#dx <- (xr - xl) / ndx
#knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
#knots
#}

bspline <- function(x, xl, xr, ndx, bdeg) {
dx <- (xr - xl) / ndx
knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
B <- spline.des(knots, x, bdeg + 1, 0 * x, outer.ok = TRUE)$design
B
}

MoE_s<-function (Y, X, G=3, delta = 1, a=1, b=0.005,v=0.01, nseg=20, bdeg=3, burn.in=0.1, maxiter = 5000, seed=303) {
	cat("Start Time =",date(),"\n")
	flush.console()
	starttime = proc.time()[3] # save the start time
	set.seed(seed) # set seed
   	N <- nrow(Y) # sample size
	X<-as.matrix(X) # concomitant covariates
	#X<-orthonormalization(X,basis=F,norm=F)
	p<-ncol(X) # number of concomitant covariates
	# re-define matrix Y using "dummy" variables 
	Y[,1]<-Y[,1]-min(Y[,1])+1 # redifine variable Y_1 in order to take values starting from 1
	Jq<-length(unique(Y[,1])) # number of possible values that the first manifest variable can take
	Y_d_old<-matrix(0,N,Jq)
	for (j in 1:N){
		Y_d_old[j,Y[j,1]]<-1 # 0-1 matrix, each row sums to one
	}
	# repeat the procedure for each other manifest variable, and bind the so-built matrices
	for(c in 2:ncol(Y)){
		Y[,c]<-Y[,c]-min(Y[,c])+1
		Jq<-c(Jq,length(unique(Y[,c])))
		Y_d<-matrix(0,N,Jq[c])
		for (j in 1:N){
			Y_d[j,Y[j,c]]<-1
		}
		Y_d_old<-cbind(Y_d_old,Y_d)
	}
	Y<-Y_d_old
	#
	if(G>1){
		# initialize the parameters of the multinomial distribution conditional on the group membership
		thetastore<-array(0,dim=c(G,ncol(Y),maxiter))
		theta <- thetastore[,,1]
		# draw from a Dirichlet, one manifest variable at a time
		theta[,1:(Jq[1])]<-rdirichlet(G,rep(delta,Jq[1])) 
		for (i in 2:length(Jq)){
			theta[,(sum(Jq[1:(i-1)])+1):(sum(Jq[1:i]))]<-rdirichlet(G,rep(delta,Jq[i]))
		}
		pr <- rdirichlet(1, rep(delta,G)) # initialize the mixing proportions drawing from a Dirichlet
		# initialize the indicator of group membership
		Sstore<-array(0, dim=c(N,G,maxiter))	
		S<-matrix(0,N,G) 
		W <- matrix(nrow = N, ncol = G) # initialize the posterior probability matrix 
		# starting allocation
		for(j in 1:N) {
     			for (g in 1:G){ 
				W[j,g] <- pr[g]*prod(theta[g,]^t(Y[j,])) #  compute posterior probabilities
			}
			S[j,sample(1:G,1,prob=W[j,]/sum(W[j,]))]<-1
		}	
		S<-S[,c(order(colSums(S))[-G],which.max(colSums(S)))]
		Z<-matrix(0,N,G) # initialize the differences in random utilities
		Lambda<-array(0,dim=c(N,N,G-1)) # initialize the weight matrix
 		for(g in 1:(G-1)){
			Lambda[,,g]<-diag(N)
		}
		for(j in 1:N){
			#repeat{
				for(g in 1:(G-1)){
					U<-runif(1)
					Z[j,g]<-log(U/(G-1)+S[j,g])-log(1-U+(1-S[j,g])/(G-1)) # draw the differences in random utilities from a truncated logistic distribution
				}
			#if(sum(Z[j,]>0)==1){break}
			#}
		}
		# the logistic distribution of the errors is approximated by a mixture of normale distributions with parameters drawn with fixed probabilities
		s<-c(0.84678,1.61,2.8904,5.0772,8.9109,15.923) # variances of the mixture of normal distributions
		omega<-c(5.8726,28.74,36.756,22.427,5.8701,0.33466)/100 # weights of the mixture of normal distributions
		w<-rep(0,6) # initialize a vector apredispongo un vettore per contenere i parametri della multinomiale
		d<-nseg+bdeg # sum the number of 
		# initialize the matrix of spline basis
		BX <- matrix(0,N,d*p) 
		for(h in 1:p){	
			# orthonormalization by Redd(2012), requires library orthogonalsplinebasis
			#nodi<-kn(x=X[,h], xl=min(X[,h]), xr=max(X[,h]), ndx=nseg, bdeg=bdeg)
			#nodi[bdeg+1]<-min(X[,h])
			#nodi[d+1]<-max(X[,h])
			#obase<-OBasis(nodi)	
			#BX[,(d*(h-1)+1):(h*d)] <- evaluate(obase,X[,h])
			# ...or, standard orthogonalization, not always appropriate. (requires library far)
			BX[,(d*(h-1)+1):(h*d)]<-bspline(x=X[,h], xl=min(X[,h]), xr=max(X[,h]), ndx=nseg, bdeg=bdeg)
		}
		# initialize the splines coefficients
		betastore<-array(0,dim=c(d*p,G,maxiter)) 
		beta<-betastore[,,1]
		# initialize the intercept (one per component)
		gammastore<-matrix(0,G,maxiter)  
		gamma<-gammastore[,1]
		Kj = INLA:::inla.rw1(d) # penalty matrix
		# initialize the smoothing parameters and the parameters of the inverse gamma distribution  
		taustore<-array(0,c(p,G-1,maxiter))
		tau2<-matrix(0,p,G-1)
		a_star<-(a+rankMatrix(Kj)/2)[1]
		b_star<-matrix(b,p,G-1)
		b<-b_star
		eta<-matrix(0,N,G) # initialize a matrix that stores the sum of the spline functions
		C<-matrix(log(G-1),N,G-1)
		m<-matrix(0,N,G-1)
		Cstore<-array(0,dim=c(N,G-1,maxiter))
		Zstore<-array(0,dim=c(N,G-1,maxiter))
		llik<-rep(0,maxiter) # initialize a vector to store the log-likelihood
		# start Gibbs sampling
		for (it in 1:maxiter) {
			for(g in 1:(G-1)){
				if(G>2){
					# vector containing log(sum(exp(predictors referred to all the groups but the g-th)))
					C[,g]<-rowLogSumExps(eta[,-g]+matrix(gamma[-g],N,G-1,byrow=TRUE))#
					#Cstore[,g,it]<-C[,g]
				}
				for(h in 1:p){
					tau2[h,g]<-rinvgamma(1,a_star,b_star[h,g]) # draw the (inverse )smoothing parameter
					taustore[h,g,it]<-tau2[h,g] # store the (inverse) smoothing parameter
					# posterior variance of the coefficients
					V<-chol2inv(chol(t(BX[,(1+d*(h-1)):(d*h)])%*%Lambda[,,g]%*%BX[,(1+d*(h-1)):(d*h)]+Kj/tau2[h,g]))#
					# posterior mean of the coefficients
					if(p==1){
						B<-V%*%t(BX[,(1+d*(h-1)):(d*h)])%*%Lambda[,,g]%*%(Z[,g]
							-gamma[g]
							+C[,g])
					}
					else{
						B<-V%*%t(BX[,(1+d*(h-1)):(d*h)])%*%Lambda[,,g]%*%(Z[,g]
							-gamma[g]
							-rowSums(BX[,-((1+d*(h-1)):(d*h))]%*%beta[-((1+d*(h-1)):(d*h)),g])
							+C[,g])
					}#
					L<-chol(V) # compute the choleski of the posterior variance
					TT<-rnorm(d) # draw from a standard distribution
					beta[(1+d*(h-1)):(d*h),g]<-as.vector(B+L%*%TT) # draw the coefficients of the splines
					# apply constrains in order to overcome identifiability issues in case p>1
					# see Appendix to Semiparametric Latent Variable Models with Bayesian P-splines 
					# Journal of Computational and Graphical Statistics, Xin-Yuan Song and Zhao-Hua Lu
			 		# page 4		
					if(p>1){
						beta[(1+d*(h-1)):(d*h),g]<-as.vector(beta[(1+d*(h-1)):(d*h),g]-
							V%*%colSums(BX[,(1+d*(h-1)):(d*h)])%*%chol2inv(chol(
							t(colSums(BX[,(1+d*(h-1)):(d*h)]))%*%V%*%colSums(BX[,(1+d*(h-1)):(d*h)])))%*%
							colSums(BX[,(1+d*(h-1)):(d*h)])%*%beta[(1+d*(h-1)):(d*h),g])
					}
					betastore[(1+d*(h-1)):(d*h),g,it]<-beta[(1+d*(h-1)):(d*h),g] # store the coefficients of the splines
					# update the parameters of the Inverse Gamma distribution (a_star is constant)
					b_star[h,g]<-(b[h,g]+(t(beta[(1+d*(h-1)):(d*h),g])%*%Kj%*%beta[(1+d*(h-1)):(d*h),g])/2)[1,1] 
				}
				# compute the predictor (without intercept)
				eta[,g]<-BX%*%beta[,g]
				if(p>1){	
					# repeat the (simplified) procedure for the intercept  			
					V<-1/(sum(diag(Lambda[,,g]))+v)
					B<-V*t(diag(Lambda[,,g]))%*%(Z[,g]-eta[,g]+C[,g])
					L<-chol(V)
					TT<-rnorm(1)
					gamma[g]<-B+L*TT
					gammastore[g,it]<-gamma[g]#	
				}
				m[,g]<-eta[,g]+gamma[g]		
			}
			if(G>2){
				for(g in 1:(G-1)){
					# vector containing log(sum(exp(predictors referred to all the groups but the g-th)))
					C[,g]<-rowLogSumExps(eta[,-g]+matrix(gamma[-g],N,G-1,byrow=TRUE))#
				}
			}
			#Cstore[,,it]<-C
			for(j in 1:N){
				#repeat{
					for(g in 1:(G-1)){
						U<-runif(1)
						# update the differences in random utilities 
						Z[j,g]<-log(exp(m[j,g]-C[j,g])*U+S[j,g])-log(1-U+exp(m[j,g]-C[j,g])*(1-S[j,g]))
					#}
				  #if(sum(Z[j,]>0)<2){break}
					#if(sum(Z[j,]>0)==1){break}
				#}
				#for(g in 1:(G-1)){	
					# update the weights of multinomial distribution
					for(h in 1:6){
						w[h]<-omega[h]*dnorm(Z[j,g]+C[j,g],m[j,g],sqrt(s[h])) 
					}
					# draw an indicator from a multinomial distribution and
					# update the weights of the mixture of normal distribution
					Lambda[j,j,g]<-1/s[sample(1:6,1,prob=w)]	
					#}
				}
			}
			#Zstore[,,it]<-Z
			# update the conditional probabilities of the manifest variables depending on group membership  
			for(g in 1:G){
				if(sum(S[,g])>1){
					theta[g,1:(Jq[1])]<-rdirichlet(1,delta+colSums(as.matrix(Y[which(S[,g]==1),1:(Jq[1])])))
					for (i in 2:length(Jq)){
						theta[g,(sum(Jq[1:(i-1)])+1):(sum(Jq[1:i]))]<-rdirichlet(1,delta+
							colSums(as.matrix(Y[which(S[,g]==1),(sum(Jq[1:(i-1)])+1):(sum(Jq[1:i]))])))
					}
				}
  				if(sum(S[,g])==0){
   					theta[g,1:(Jq[1])]<-rdirichlet(1,rep(delta,Jq[1]))
   					for (i in 2:length(Jq)){
    						theta[g,(sum(Jq[1:(i-1)])+1):(sum(Jq[1:i]))]<-rdirichlet(1,rep(delta,Jq[i]))
   					}
  				}
			}
			thetastore[,,it]<-theta # store the conditional probabilities
			# initialize the indicator of group membership
			S<-matrix(0,N,G)
			for(j in 1:N){
				# compute the posterior probabilities for the j-th unit
				for(g in 1:G){
					W[j,g] <- exp(eta[j,g]+gamma[g])/ 
						sum(exp(eta[j,]+gamma))*
						prod(theta[g,]^t(Y[j,]))
				}
				S[j,sample(1:G,1,prob=W[j,]/sum(W[j,]))]<-1	# allocate the j-th unit
				# compute the log-likelihood for the first j units
				llik[it]<-llik[it]+log(prod(theta[which(S[j,]==1),]^t(Y[j,]))) 
			}
			#
			Sstore[,,it]<-S # store the allocations
			if(is.element(it, c(1:5, 10, 20, 50, 100, 200)) | it%%500==0){
				if((proc.time()[3] - starttime)<60){
					cat("sim =",
					it, 
					"/ duration of iter proc so far:", 
                  		round(diff <- proc.time()[3] - starttime, 2),
                  		"sec.\n")
					flush.console()
				}
				else{
					cat("sim =",
						it, 
						"/ duration of iter proc so far:", 
                  			round(diff <- ((proc.time()[3] - starttime)/60), 2),
                  			"min.\n")
						flush.console()
				}
			}
			#if (it>29){
			#	if(sum(effectiveSize(as.mcmc(t(gammastore[-G,1:it])))>minESS(1))==G-1){
			#		break
			#		maxiter<-it
			#	}
			#}
		}
		burn.in<-round(burn.in*maxiter)
		group<-rep(0,N)
		# apply MAP rule
		S<-apply(Sstore[,,(burn.in+1):maxiter],c(1,2),mean)
		for(j in 1:N){
			group[j]<-which.max(S[j,])
		}
		if(p==1 & G==2){mean(taustore[,,(burn.in+1):maxiter])}
		else{tau2<-apply(taustore[,,(burn.in+1):maxiter],c(1,2),mean)}
		finish = proc.time()[3]
		duration = finish - starttime
		output<-list(beta=apply(betastore[,,(burn.in+1):maxiter],c(1,2),mean),
			gamma=rowMeans(gammastore[,(burn.in+1):maxiter]),
			tau2=tau2,
			theta=apply(thetastore[,,(burn.in+1):maxiter],c(1,2),mean),
			#AICM=aicm(llik[(burn.in+1):maxiter]),
			AICM=2*(mean(llik[(burn.in+1):maxiter])-var(llik[(burn.in+1):maxiter])),
			group=group,
			betastore=betastore,
			gammastore=gammastore,
			taustore=taustore,
			thetastore=thetastore,
			llik=llik,
			BX=BX,
			#Sstore=Sstore,
			#Zstore=Zstore,
			#Cstore=Cstore,
			duration=duration
		)
		return(output)
	}
	# special case: G=1
	if(G==1){
		thetastore<-matrix(0,maxiter,ncol(Y))
		theta <- thetastore[1,]
		theta[1:(Jq[1])]<-rdirichlet(1,rep(delta,Jq[1]))
		for (i in 2:length(Jq)){
			theta[(sum(Jq[1:(i-1)])+1):(sum(Jq[1:i]))]<-rdirichlet(1,rep(delta,Jq[i]))
		}
		llik<-rep(0,maxiter)
		W<-rep(0,N)
		for (it in 1:maxiter) {
			theta[1:(Jq[1])]<-rdirichlet(1,delta+colSums(Y[,1:(Jq[1])]))
			for (i in 2:length(Jq)){
				theta[(sum(Jq[1:(i-1)])+1):(sum(Jq[1:i]))]<-rdirichlet(1,delta+
					colSums(Y[,(sum(Jq[1:(i-1)])+1):(sum(Jq[1:i]))]))
			}
			thetastore[it,]<-theta
			for(j in 1:N){
				W[j] <- prod(theta^t(Y[j,]))
			}
			llik[it]<-sum(log(W))
			if(is.element(it, c(1:5, 10, 20, 50, 100, 200)) | it%%500==0){
				if((proc.time()[3] - starttime)<60){
					cat("sim =",
					it, 
					"/ duration of iter proc so far:", 
                  		round(diff <- proc.time()[3] - starttime, 2),
                  		"sec.\n")
					flush.console()
				}
				else{
					cat("sim =",
					it, 
					"/ duration of iter proc so far:", 
                  		round(diff <- ((proc.time()[3] - starttime)/60), 2),
                  		"min.\n")
					flush.console()
				}
			}
		}
		finish = proc.time()[3]
		duration = finish - starttime
		output<-list(theta=colSums(thetastore),
			#AICM=aicm(llik[(burn.in+1):maxiter]),
			AICM=2*(mean(llik[(burn.in+1):maxiter])-var(llik[(burn.in+1):maxiter])),
			thetastore=thetastore,
			llik=llik,
			duration=duration)
		return(output)
	}
}

