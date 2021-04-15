bspline <- function(x, xl, xr, ndx, bdeg) {
dx <- (xr - xl) / ndx
knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
B <- spline.des(knots, x, bdeg + 1, 0 * x, outer.ok = TRUE)$design
B
}

MoE_sN<-function (Y, X, G=3, delta = 1, a=1, b=0.005, a.n=2.5, b.n=0.5, v=0.01, nseg=20, bdeg=3, burn.in=0.1, maxiter = 5000, seed=303) {
	cat("Start Time =",date(),"\n")
	flush.console()
	starttime = proc.time()[3] 
	set.seed(seed)
	X<-as.matrix(X) 
	#X<-orthonormalization(X,basis=F,norm=F)
	p<-ncol(X) 
	Y<-as.matrix(Y)
	Y<-as.matrix(Y);J<-ncol(Y);N <- nrow(Y);mu0=colMeans(Y) 
	if(G>1){
		S <- matrix(0,nrow = N, ncol = G) 
		Sstore<-array(0, dim=c(N,G,maxiter))
		prstore<-matrix(0,G,maxiter)
		pr <- rdirichlet(1, rep(delta,G)) 
		if(J==1){
			v.n<-1#2.6/(diff(range(Y)))^2
			sigma2store<-matrix(0,G,maxiter)
			mustore<-matrix(0,G,maxiter)
			sigma2<-rep(var(Y)/G^2,G)
			mu<-rnorm(G,rep(mu0,G),sqrt(sigma2))
			W<-matrix(0,N,G)
			for(i in 1:N) {
				S[i,sample(1:G,size=1,prob=pr*dnorm(Y[i],mu,sqrt(sigma2)))]<-1	 
			}
			S<-S[,c(order(colSums(S))[-G],which.max(colSums(S)))]
		}
		if(J>1){
			v.n<-1#2.6/(diff(range(Y)))^2
			sigma2store<-array(0,dim=c(J,J,G,maxiter))
			sigma2<-sigma2store[,,,1]
			for(g in 1:G){sigma2[,,g]<-var(Y)/G^2}
			mustore<-array(0,dim=c(J,G,maxiter))
			mu<-t(rmnorm(G,mu0,sigma2[,,1]))
			W<-matrix(0,N,G)
			for(i in 1:N) {
     				for (g in 1:G){ 
					W[i,g] <- pr[g]*dmnorm(Y[i,],mu[,g],sigma2[,,g]) 
				}
				S[i,sample(1:G,size=1,prob=W[i,])]<-1	
			}
			S<-S[,c(order(colSums(S))[-G],which.max(colSums(S)))]
		} 
		Z<-matrix(0,N,G-1) 
		Lambda<-array(0,dim=c(N,N,G-1))
 		for(g in 1:(G-1)){
			Lambda[,,g]<-diag(N)
			for(j in 1:N){
				U<-runif(1)
				Z[j,g]<-log(U/(G-1)+S[j,g])-log(1-U+(1-S[j,g])/(G-1)) 
			}
		}
		s<-c(0.84678,1.61,2.8904,5.0772,8.9109,15.923) 
		omega<-c(5.8726,28.74,36.756,22.427,5.8701,0.33466)/100 
		#s<-c(0.68159,1.2419,2.2388,4.0724,7.4371,13.772) 
		#omega<-c(1.8446,17.268,37.393,31.697,10.89,0.90745)/100 
		#s<-c(1.2131,2.9955,7.5458,16.257) 
		#omega<-c(25.22,58.523,16.257)/100 
		w<-rep(0,length(omega)) 
		d<-nseg+bdeg 
		BX <- matrix(0,N,d*p) 
		for(h in 1:p){	
			BX[,(d*(h-1)+1):(h*d)]<-bspline(x=X[,h], xl=min(X[,h]), xr=max(X[,h]), ndx=nseg, bdeg=bdeg)
		}
		betastore<-array(0,dim=c(d*p,G,maxiter)) 
		beta<-betastore[,,1]
		gammastore<-matrix(0,G,maxiter)  
		gamma<-gammastore[,1]
		Kj = INLA:::inla.rw1(d)+diag(1e-3,d)
		taustore<-array(0,c(p,G-1,maxiter))
		tau2<-matrix(0,p,G-1)
		a_star<-(a+rankMatrix(Kj)/2)[1]
		b_star<-matrix(b,p,G-1)
		C<-matrix(log(G-1),N,G-1)
		m<-matrix(0,N,G)
		llik<-rep(0,maxiter) 
		for (it in 1:maxiter) {
			for(g in 1:(G-1)){
				if(G>2){
					C[,g]<-rowLogSumExps(m[,-g])
				}
				for(h in 1:p){
					tau2<-rgamma(1,shape=a_star,rate=b_star[h,g]) 
					taustore[h,g,it]<-tau2
					V<-chol2inv(chol(t(BX[,(1+d*(h-1)):(d*h)])%*%Lambda[,,g]%*%BX[,(1+d*(h-1)):(d*h)]+Kj*tau2))
					if(p==1){
						B<-V%*%t(BX[,(1+d*(h-1)):(d*h)])%*%Lambda[,,g]%*%(Z[,g]
							-gamma[g]
							+C[,g])
					}
					else{
						B<-V%*%t(BX[,(1+d*(h-1)):(d*h)])%*%Lambda[,,g]%*%(Z[,g]
							-gamma[g]
							-BX[,-((1+d*(h-1)):(d*h))]%*%beta[-((1+d*(h-1)):(d*h)),g]
							+C[,g])
					}
					L<-chol(V) 
					TT<-rnorm(d) 
					beta[(1+d*(h-1)):(d*h),g]<-as.vector(B+L%*%TT) 
					if(p>1){
						beta[(1+d*(h-1)):(d*h),g]<-as.vector(beta[(1+d*(h-1)):(d*h),g]-
							V%*%colSums(BX[,(1+d*(h-1)):(d*h)])%*%chol2inv(chol(
							t(colSums(BX[,(1+d*(h-1)):(d*h)]))%*%V%*%colSums(BX[,(1+d*(h-1)):(d*h)])))%*%
							colSums(BX[,(1+d*(h-1)):(d*h)])%*%beta[(1+d*(h-1)):(d*h),g])
					}
					betastore[(1+d*(h-1)):(d*h),g,it]<-beta[(1+d*(h-1)):(d*h),g] 
					b_star[h,g]<-(b+(t(beta[(1+d*(h-1)):(d*h),g])%*%Kj%*%beta[(1+d*(h-1)):(d*h),g])/2)[1,1] 
				}
				if(p>1){
					V<-1/(sum(diag(Lambda[,,g]))+v)
					B<-V*t(diag(Lambda[,,g]))%*%(Z[,g]-BX%*%beta[,g]+C[,g])
					L<-chol(V)
					TT<-rnorm(1)
					gamma[g]<-B+L*TT
					gammastore[g,it]<-gamma[g]
				}
				m[,g]<-BX%*%beta[,g]+gamma[g]	
			}
			if(G>2){
				for(g in 1:(G-2)){
					C[,g]<-rowLogSumExps(m[,-g])
				}
			}
			for(j in 1:N){
				for(g in 1:(G-1)){
					U<-runif(1)
					Z[j,g]<-log(exp(m[j,g]-C[j,g])*U+S[j,g])-log(1-U+exp(m[j,g]-C[j,g])*(1-S[j,g]))
					for(h in 1:length(omega)){
						w[h]<-omega[h]*dnorm(Z[j,g]+C[j,g],m[j,g],sqrt(s[h])) 
					}
					Lambda[j,j,g]<-1/s[sample(1:length(omega),1,prob=w)]	
				}
			}
			# update the conditional probabilities of the manifest variables depending on group membership  
			if(J==1){
				repeat{
					mu<-rnorm(G,
						(mu0*v.n+colSums(Y*S))/(v.n+colSums(S)),
						sqrt(sigma2/(v.n+colSums(S)))
					)
					if (max(mu)<max(Y) & min(mu)>min(Y)){break}
				}
				sigma2<-1/rgamma(G,
					shape=a.n+(colSums(S)+1)/2,
					rate=	
						b.n*var(Y)
						+0.5*colSums((Y-matrix(mu,N,G,byrow=TRUE))^2*S)
						+0.5*v.n*(mu-mu0)^2
					)
				sigma2store[,it]<-sigma2 
				mustore[,it]<-mu
				S<-matrix(0,N,G)
				# compute the posterior probabilities and allocate units into groups drawing from a multinomial
				for(i in 1:N){
					for(g in 1:G){
						W[i,g] <- exp(m[i,g])/ 
							sum(exp(m[i,]))*
							dnorm(Y[i],mu[g],sqrt(sigma2[g]))
					}
					S[i,sample(1:G,1,prob=W[i,])]<-1
					llik[it]<-llik[it]+log(dnorm(as.numeric(Y[i]),mu[which(S[i,]==1)],sqrt(sigma2[which(S[i,]==1)])))
				}
			}
			if(J>1){
				for(g in 1:G){
					repeat{
						mu[,g]<-rmnorm(1,
							(mu0*v.n+colSums(Y*matrix(S[,g],N,J)))/(v.n+(sum(S[,g]))),
							sigma2[,,g]/(v.n+(sum(S[,g])))
						)
						if (sum(mu[,g]<apply(Y,2,max))==J & sum(mu[,g]>apply(Y,2,min))==J){break}
					}
					#if(sum(S[,g])==0){
					#	sigma2[,,g]<-riwish(J*a.n+1,b.n*var(Y)+v.n*(mu0-mu[,g])%*%t(mu0-mu[,g]))
					#}
					#if(sum(S[,g])==1){
					#	sigma2[,,g]<-riwish(J*a.n+2,
					#		J*b.n*var(Y)+v.n/(v.n+1)*(Y[which(S[,g]==1),]-mu0)%*%t(Y[which(S[,g]==1),]-mu0))
					#}
					if(sum(S[,g])>1){
						sigma2[,,g]<-riwish(
							a.n+0.5*(sum(S[,g])+1),
							#chol2inv(chol(
								b.n*var(Y)
								#+t(Y*matrix(S[,g],N,J)-matrix(colMeans(Y[which(S[,g]==1),]),N,J,byrow=TRUE)*matrix(S[,g],N,J))%*%
								#(Y*matrix(S[,g],N,J)-matrix(colMeans(Y[which(S[,g]==1),]),N,J,byrow=TRUE)*matrix(S[,g],N,J))
								#+v.n*sum(S[,g])/(v.n+sum(S[,g]))*(colMeans(Y[which(S[,g]==1),])-mu0)%*%t(colMeans(Y[which(S[,g]==1),])-mu0)
								+0.5*t(Y*matrix(S[,g],N,J)-matrix(mu[,g],N,J,byrow=TRUE)*matrix(S[,g],N,J))%*%
								(Y*matrix(S[,g],N,J)-matrix(mu[,g],N,J,byrow=TRUE)*matrix(S[,g],N,J))
								+0.5*v.n*(mu0-mu[,g])%*%t(mu0-mu[,g])#))
						)
						#sigma2[,,g]<-chol2inv(chol(sigma2[,,g]))
					}
				}	
				sigma2store[,,,it]<-sigma2
				mustore[,,it]<-mu
				S<-matrix(0,N,G)
				for(i in 1:N){
					for(g in 1:G){
						W[i,g] <- exp(m[i,g])/ 
							sum(exp(m[i,]))*
							dmnorm(as.numeric(Y[i,]),mu[,g],sigma2[,,g])

					}
					S[i,sample(1:G,1,prob=W[i,])]<-1
					llik[it]<-llik[it]+log(dmnorm(as.numeric(Y[i,]),mu[,which(S[i,]==1)],sigma2[,,which(S[i,]==1)]))	
				}
			}
			Sstore[,,it]<-S 
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
		burn.in<-round(burn.in*maxiter)
		group<-rep(0,N)
		S<-apply(Sstore[,,(burn.in+1):maxiter],c(1,2),mean)
		for(j in 1:N){
			group[j]<-which.max(S[j,])
		}
		finish = proc.time()[3]
		duration = finish - starttime
		output<-list(AICM=-2*(mean(llik[(burn.in+1):maxiter])-var(llik[(burn.in+1):maxiter])),
			group=group,
			betastore=betastore,
			gammastore=gammastore,
			taustore=taustore,
			sigma2store=sigma2store,
			mustore=mustore,
			llik=llik,
			BX=BX,
			duration=duration
		)
		return(output)
	}
	if(G==1){
		if(J==1){
			v.n<-2.6/(diff(range(Y)))^2
			sigma2store<-rep(0,maxiter)
			mustore<-rep(0,maxiter)
			sigma2<-var(Y)
			mu<-rnorm(1,mu0,sqrt(sigma2))
		}
		if(J>1){
			v.n<-2.6/(diff(range(Y)))^2
			sigma2store<-array(0,dim=c(J,J,maxiter))
			sigma2<-var(Y)
			mustore<-matrix(0,J,maxiter)
			mu<-t(rmnorm(1,mu0,sigma2))
		}
		llik=rep(0,maxiter) 
		for(it in 1:maxiter){
			if(J==1){
				repeat{
					mu<-rnorm(1,(mu0*v.n+colSums(Y))/(v.n+N),
						sqrt(sigma2/(v.n+N))
					)
				}
				if (mu<max(Y) & mu>min(Y)){break}
				sigma2<-1/rgamma(1,
					shape=a.n+(N+1)/2,
					rate=b.n*var(Y)+0.5*colSums((Y-mu)^2)+0.5*v.n*(mu-mu0)^2
				)
				sigma2store[it]<-sigma2 
				mustore[it]<-mu
				llik[it]<-sum(log(dnorm(Y,mu,sqrt(sigma2))))
			}
			if(J>1){
				repeat{
					mu<-rmnorm(1,colSums(mu0*v.n+Y)/(v.n+N),sigma2/(v.n+N))
					if (sum(mu<apply(Y,2,max))==J & sum(mu>apply(Y,2,min))==J){break}
				}
				sigma2<-riwish(
						J*a.n+N+1,
						J*b.n*var(Y)
						+t(Y-matrix(mu,N,J,byrow=TRUE))%*%
						(Y-matrix(mu,N,J,byrow=TRUE))
						+v.n*t(mu-mu0)%*%(mu0-mu)
						#+t(Y-matrix(colMeans(Y),N,J,byrow=TRUE))%*%
						#(Y-matrix(colMeans(Y),N,J,byrow=TRUE))
						#+v.n*N/(v.n+N)*(colMeans(Y)-mu0)%*%t(colMeans(Y)-mu0)
				)
				sigma2store[,,it]<-sigma2
				mustore[,it]<-mu
				S<-matrix(0,N,G)
				for(i in 1:N){llik[it]<-llik[it]+log(dmnorm(Y[i,],t(mu),sigma2))}	
			}
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
		burn.in<-round(burn.in*maxiter)
		finish = proc.time()[3]
		duration = finish - starttime
		output<-list(AICM=-2*(mean(llik[(burn.in+1):maxiter])-var(llik[(burn.in+1):maxiter])),
			sigma2store=sigma2store,
			mustore=mustore,
			llik=llik,
			duration=duration
		)
		return(output)
	}
}


mix_N<-function (Y, G=3, delta = 1, a.n=2.5, b.n=0.5, burn.in=0.1, maxiter = 5000, seed=303) {
	cat("Start Time =",date(),"\n")
	flush.console()
	starttime = proc.time()[3] 
	set.seed(seed)
	Y<-as.matrix(Y)
	Y<-as.matrix(Y);J<-ncol(Y);N <- nrow(Y);mu0=colMeans(Y) 
	if(G>1){
		S <- matrix(0,nrow = N, ncol = G) 
		Sstore<-array(0, dim=c(N,G,maxiter))
		prstore<-matrix(0,G,maxiter)
		pr <- rdirichlet(1, rep(delta,G)) 
		if(J==1){
			v.n<-1
			sigma2store<-matrix(0,G,maxiter)
			mustore<-matrix(0,G,maxiter)
			sigma2<-rep(var(Y)/G^2,G)
			mu<-rnorm(G,rep(mu0,G),sqrt(sigma2))
			W<-matrix(0,N,G)
			for(i in 1:N) {
				S[i,sample(1:G,size=1,prob=pr*dnorm(Y[i],mu,sqrt(sigma2)))]<-1	 
			}
			S<-S[,c(order(colSums(S))[-G],which.max(colSums(S)))]
		}
		if(J>1){
			v.n<-1#2.6/(diff(range(Y)))^2
			sigma2store<-array(0,dim=c(J,J,G,maxiter))
			sigma2<-sigma2store[,,,1]
			for(g in 1:G){sigma2[,,g]<-var(Y)/G^2}
			mustore<-array(0,dim=c(J,G,maxiter))
			mu<-t(rmnorm(G,mu0,sigma2[,,1]))
			W<-matrix(0,N,G)
			for(i in 1:N) {
     				for (g in 1:G){ 
					W[i,g] <- pr[g]*dmnorm(Y[i,],mu[,g],sigma2[,,g]) 
				}
				S[i,sample(1:G,size=1,prob=W[i,])]<-1	
			}
			S<-S[,c(order(colSums(S))[-G],which.max(colSums(S)))]
		} 
		llik<-rep(0,maxiter) 
		for (it in 1:maxiter) {
			pr <- rdirichlet(1, delta+colSums(S))
			# update the conditional probabilities of the manifest variables depending on group membership  
			if(J==1){
				repeat{
					mu<-rnorm(G,
						(mu0*v.n+colSums(Y*S))/(v.n+colSums(S)),
						sqrt(sigma2/(v.n+colSums(S)))
					)
					if (max(mu)<max(Y) & min(mu)>min(Y)){break}
				}
				sigma2<-1/rgamma(G,
					shape=a.n+(colSums(S)+1)/2,
					rate=	
						b.n*var(Y)
						+0.5*colSums((Y-matrix(mu,N,G,byrow=TRUE))^2*S)
						+0.5*v.n*(mu-mu0)^2
					)
				sigma2store[,it]<-sigma2 
				mustore[,it]<-mu
				S<-matrix(0,N,G)
				# compute the posterior probabilities and allocate units into groups drawing from a multinomial
				for(i in 1:N){
					for(g in 1:G){
						W[i,g] <- pr[g]*
							dnorm(Y[i],mu[g],sqrt(sigma2[g]))
					}
					S[i,sample(1:G,1,prob=W[i,])]<-1
					llik[it]<-llik[it]+log(dnorm(as.numeric(Y[i]),mu[which(S[i,]==1)],sqrt(sigma2[which(S[i,]==1)])))
				}
			}
			if(J>1){
				for(g in 1:G){
					repeat{
						mu[,g]<-rmnorm(1,
							(mu0*v.n+colSums(Y*matrix(S[,g],N,J)))/(v.n+(sum(S[,g]))),
							sigma2[,,g]/(v.n+(sum(S[,g])))
						)
						if (sum(mu[,g]<apply(Y,2,max))==J & sum(mu[,g]>apply(Y,2,min))==J){break}
					}
					#if(sum(S[,g])==0){
					#	sigma2[,,g]<-riwish(J*a.n+1,b.n*var(Y)+v.n*(mu0-mu[,g])%*%t(mu0-mu[,g]))
					#}
					#if(sum(S[,g])==1){
					#	sigma2[,,g]<-riwish(J*a.n+2,
					#		J*b.n*var(Y)+v.n/(v.n+1)*(Y[which(S[,g]==1),]-mu0)%*%t(Y[which(S[,g]==1),]-mu0))
					#}
					if(sum(S[,g])>1){
						sigma2[,,g]<-riwish(
							a.n+0.5*(sum(S[,g])+1),
							#chol2inv(chol(
								b.n*var(Y)
								#+t(Y*matrix(S[,g],N,J)-matrix(colMeans(Y[which(S[,g]==1),]),N,J,byrow=TRUE)*matrix(S[,g],N,J))%*%
								#(Y*matrix(S[,g],N,J)-matrix(colMeans(Y[which(S[,g]==1),]),N,J,byrow=TRUE)*matrix(S[,g],N,J))
								#+v.n*sum(S[,g])/(v.n+sum(S[,g]))*(colMeans(Y[which(S[,g]==1),])-mu0)%*%t(colMeans(Y[which(S[,g]==1),])-mu0)
								+0.5*t(Y*matrix(S[,g],N,J)-matrix(mu[,g],N,J,byrow=TRUE)*matrix(S[,g],N,J))%*%
								(Y*matrix(S[,g],N,J)-matrix(mu[,g],N,J,byrow=TRUE)*matrix(S[,g],N,J))
								+0.5*v.n*(mu0-mu[,g])%*%t(mu0-mu[,g])#))
						)
						#sigma2[,,g]<-chol2inv(chol(sigma2[,,g]))
					}
				}	
				sigma2store[,,,it]<-sigma2
				mustore[,,it]<-mu
				S<-matrix(0,N,G)
				for(i in 1:N){
					for(g in 1:G){
						W[i,g] <-pr[g]*
							dmnorm(as.numeric(Y[i,]),mu[,g],sigma2[,,g])

					}
					S[i,sample(1:G,1,prob=W[i,])]<-1
					llik[it]<-llik[it]+log(dmnorm(as.numeric(Y[i,]),mu[,which(S[i,]==1)],sigma2[,,which(S[i,]==1)]))	
				}
			}
			Sstore[,,it]<-S 
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
		burn.in<-round(burn.in*maxiter)
		group<-rep(0,N)
		S<-apply(Sstore[,,(burn.in+1):maxiter],c(1,2),mean)
		for(j in 1:N){
			group[j]<-which.max(S[j,])
		}
		finish = proc.time()[3]
		duration = finish - starttime
		output<-list(AICM=-2*(mean(llik[(burn.in+1):maxiter])-var(llik[(burn.in+1):maxiter])),
			group=group,
			sigma2store=sigma2store,
			mustore=mustore,
			llik=llik,
			duration=duration
		)
		return(output)
	}
	if(G==1){
		if(J==1){
			v.n<-1#2.6/(diff(range(Y)))^2
			sigma2store<-rep(0,maxiter)
			mustore<-rep(0,maxiter)
			sigma2<-var(Y)
			#mu<-rnorm(1,mu0,sqrt(sigma2))
		}
		if(J>1){
			v.n<-1#2.6/(diff(range(Y)))^2
			sigma2store<-array(0,dim=c(J,J,maxiter))
			sigma2<-var(Y)
			mustore<-matrix(0,J,maxiter)
			#mu<-t(rmnorm(1,mu0,sigma2))
		}
		llik=rep(0,maxiter) 
		for(it in 1:maxiter){
			if(J==1){
				repeat{
					mu<-rnorm(1,(mu0*v.n+colSums(Y))/(v.n+N),
						sqrt(sigma2/(v.n+N))
					)
				if (mu<max(Y) & mu>min(Y)){break}
				}
				sigma2<-1/rgamma(1,
					shape=a.n+(N+1)/2,
					rate=b.n*var(Y)+0.5*colSums((Y-mu)^2)+0.5*v.n*(mu-mu0)^2
				)
				sigma2store[it]<-sigma2 
				mustore[it]<-mu
				llik[it]<-sum(log(dnorm(Y,mu,sqrt(sigma2))))
			}
			if(J>1){
				repeat{
					mu<-rmnorm(1,(mu0*v.n+colSums(Y))/(v.n+N),sigma2/(v.n+N))
					if (sum(mu<apply(Y,2,max))==J & sum(mu>apply(Y,2,min))==J){break}
				}
				sigma2<-riwish(
						a.n+(N+1)/2,
						b.n*var(Y)
						+0.5*t(Y-matrix(mu,N,J,byrow=TRUE))%*%
						(Y-matrix(mu,N,J,byrow=TRUE))
						+0.5*v.n*t(mu-mu0)%*%(mu0-mu)
						#+t(Y-matrix(colMeans(Y),N,J,byrow=TRUE))%*%
						#(Y-matrix(colMeans(Y),N,J,byrow=TRUE))
						#+v.n*N/(v.n+N)*(colMeans(Y)-mu0)%*%t(colMeans(Y)-mu0)
				)
				sigma2store[,,it]<-sigma2
				mustore[,it]<-t(mu)
				S<-matrix(0,N,G)
				for(i in 1:N){llik[it]<-llik[it]+log(dmnorm(Y[i,],t(mu),sigma2))}	
			}
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
		burn.in<-round(burn.in*maxiter)
		finish = proc.time()[3]
		duration = finish - starttime
		output<-list(AICM=-2*(mean(llik[(burn.in+1):maxiter])-var(llik[(burn.in+1):maxiter])),
			sigma2store=sigma2store,
			mustore=mustore,
			llik=llik,
			duration=duration
		)
		return(output)
	}
}
