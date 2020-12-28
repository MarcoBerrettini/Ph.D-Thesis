

bspline <- function(x, xl, xr, ndx, bdeg) {
dx <- (xr - xl) / ndx
knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
B <- spline.des(knots, x, bdeg + 1, 0 * x, outer.ok = TRUE)$design
B
}

full_nl<-function (Y, X, G=3, delta = 1, a=1, b=0.005, a.n=1.28, b.n=0.36, v=0.01, nseg=20, bdeg=3, eps=0, burn.in=0.5, maxiter = 8000, seed=303) {
	cat("Start Time =",date(),"\n")
	flush.console()
	starttime = proc.time()[3] 
	set.seed(seed) 
	maxiter0<-maxiter
	X<-as.matrix(X) 
	p<-ncol(X) 
	Y<-as.numeric(Y)
	N<-length(Y)
	mu0=mean(Y)
	d<-nseg+bdeg  
	BX <- matrix(0,N,d*p) 
	for(h in 1:p){	
		BX[,(d*(h-1)+1):(h*d)]<-bspline(x=X[,h], xl=min(X[,h]), xr=max(X[,h]), ndx=nseg, bdeg=bdeg)
	}
	Kj = INLA:::inla.rw1(d) 
	Kj2<-Kj<-Kj+diag(eps,d)
	if(G>1){
		S <- matrix(0,nrow = N, ncol = G) 
		Sstore<-array(0, dim=c(N,G,maxiter+maxiter-(burn.in*maxiter)))
		prstore<-matrix(0,G,maxiter+maxiter-(burn.in*maxiter))
		pr <- rdirichlet(1, rep(delta,G)) 
		v.n<-2.6/(diff(range(Y)))^2
		sigma2store<-matrix(0,G,maxiter+maxiter-(burn.in*maxiter))
		mustore<-array(0,dim=c(N,G,maxiter+maxiter-(burn.in*maxiter)))
		sigma2<-rep(var(Y)/G^2,G)
		mu<-matrix(rnorm(G,rep(mu0,G),sqrt(sigma2)),N,G,byrow=T)
		W<-matrix(0,N,G)
		for(i in 1:N) {
			S[i,sample(1:G,size=1,prob=pr*dnorm(Y[i],mu[i,],sqrt(sigma2)))]<-1	 
		}
		S<-S[,c(order(colSums(S))[-G],which.max(colSums(S)))]
		betamustore<-array(0,dim=c(d*p,G,maxiter+maxiter-(burn.in*maxiter))) 
		betamu<-betamustore[,,1]
		gammamustore<-matrix(0,G,maxiter+maxiter-(burn.in*maxiter))  
		gammamu<-gammamustore[,1]
		taumustore<-array(0,c(p,G,maxiter+maxiter-(burn.in*maxiter)))
		taumu2<-matrix(rgamma(p*G,shape=a,rate=b),p,G)
		bmu_star<-matrix(b,p,G)
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
		w<-rep(0,6) 
		betastore<-array(0,dim=c(d*p,G,maxiter+maxiter-(burn.in*maxiter))) 
		beta<-betastore[,,1]
		gammastore<-matrix(0,G,maxiter+maxiter-(burn.in*maxiter))  
		gamma<-gammastore[,1]
		taustore<-array(0,c(p,G-1,maxiter+maxiter-(burn.in*maxiter)))
		tau2<-matrix(rgamma(p*(G-1),shape=a,rate=b),p,G-1)
		a_star<-(a+rankMatrix(Kj)/2)[1]
		b_star<-matrix(b,p,G-1)
		eta<-matrix(0,N,G) 
		C<-matrix(log(G-1),N,G-1)
		m<-matrix(0,N,G-1)
		llik<-rep(0,maxiter+maxiter-(burn.in*maxiter)) 
		it<-0
		while (it <= maxiter) {
			it<-it+1
			for(g in 1:(G-1)){
				if(G>2){
					C[,g]<-rowLogSumExps(eta[,-g]+matrix(gamma[-g],N,G-1,byrow=TRUE))
				}
				for(h in 1:p){
					V<-chol2inv(chol(t(BX[,(1+d*(h-1)):(d*h)])%*%Lambda[,,g]%*%BX[,(1+d*(h-1)):(d*h)]+Kj*tau2[h,g]))
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
					tau2[h,g]<-rgamma(1,shape=a_star,rate=b_star[h,g]) 
					taustore[h,g,it]<-tau2[h,g] 
				}
				eta[,g]<-BX%*%beta[,g]
				if(p>1){
					V<-1/(sum(diag(Lambda[,,g]))+v)
					B<-V*t(diag(Lambda[,,g]))%*%(Z[,g]-eta[,g]+C[,g])
					L<-chol(V)
					TT<-rnorm(1)
					gamma[g]<-B+L*TT
					gammastore[g,it]<-gamma[g]
				}
				m[,g]<-eta[,g]+gamma[g]	
			}
			if(G>2){
				for(g in 1:(G-2)){
					C[,g]<-rowLogSumExps(eta[,-g]+matrix(gamma[-g],N,G-1,byrow=TRUE))
				}
			}
			for(j in 1:N){
				for(g in 1:(G-1)){
					U<-runif(1)
					Z[j,g]<-log(exp(m[j,g]-C[j,g])*U+S[j,g])-log(1-U+exp(m[j,g]-C[j,g])*(1-S[j,g]))
					for(h in 1:6){
						w[h]<-omega[h]*dnorm(Z[j,g]+C[j,g],m[j,g],sqrt(s[h])) 
					}
					Lambda[j,j,g]<-1/s[sample(1:6,1,prob=w)]	
				}
			}
			for(g in 1:G){	
				if(sum(S[,g])>1){	
					for(h in 1:p){
						V<-chol2inv(chol(t(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]*sigma2[g]+Kj2*taumu2[h,g]))
						if(p==1){
							B<-V%*%t(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%(Y[which(S[,g]==1)]-gammamu[g])*sigma2[g]
						}
						else{
							B<-V%*%t(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%(Y[which(S[,g]==1)]-gammamu[g]-
								BX[which(S[,g]==1),-((1+d*(h-1)):(d*h))]%*%betamu[-((1+d*(h-1)):(d*h)),g])*sigma2[g]
						}
						L<-chol(V) 
						T<-rnorm(d)
						betamu[(1+d*(h-1)):(d*h),g]<-as.vector(B+L%*%T)
						if(p>1){ 
							betamu[(1+d*(h-1)):(d*h),g]<-as.vector(betamu[(1+d*(h-1)):(d*h),g]-
								V%*%colSums(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%chol2inv(chol(
								t(colSums(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]))%*%V%*%colSums(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])))%*%
								colSums(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%betamu[(1+d*(h-1)):(d*h),g])
						}
						betamustore[(1+d*(h-1)):(d*h),g,it]<-betamu[(1+d*(h-1)):(d*h),g] 
						bmu_star[h,g]<-(b+(t(betamu[(1+d*(h-1)):(d*h),g])%*%Kj2%*%betamu[(1+d*(h-1)):(d*h),g])/2)[1,1]
						taumu2[h,g]<-rgamma(1,shape=a_star,rate=bmu_star[h,g]) 
						taumustore[h,g,it]<-taumu2[h,g] 
					}
					if(p>1){
						V<-1/(sum(S[,g])*sigma2[g]+v)
						B<-V*sum(Y[which(S[,g]==1)]-BX[which(S[,g]==1),]%*%betamu[,g])*sigma2[g]
						L<-sqrt(V)
						T<-rnorm(1)
						gammamu[g]<-B+L*T
						gammamustore[g,it]<-gammamu[g]
					}
				}
				if(sum(S[,g])==1){
					for(h in 1:p){
						V<-ginv(as.matrix(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]%*%t(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])*sigma2[g]+Kj2*taumu2[h,g]))
						if(p==1){
							B<-V%*%(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])*(Y[which(S[,g]==1)]-gammamu[g])*sigma2[g]
						}
						else{
							B<-V%*%(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%(Y[which(S[,g]==1)]-gammamu[g]-
								BX[which(S[,g]==1),-((1+d*(h-1)):(d*h))]%*%betamu[-((1+d*(h-1)):(d*h)),g])*sigma2[g]
						}
						L<-chol(V,pivot=TRUE) 
						T<-rnorm(d)
						betamu[(1+d*(h-1)):(d*h),g]<-as.vector(B+L%*%T)
						if(p>1){ 
							betamu[(1+d*(h-1)):(d*h),g]<-as.vector(betamu[(1+d*(h-1)):(d*h),g]-
								V%*%BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]%*%chol2inv(chol(
								t(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%V%*%BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]))%*%
								BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]%*%betamu[(1+d*(h-1)):(d*h),g])
						}
						betamustore[(1+d*(h-1)):(d*h),g,it]<-betamu[(1+d*(h-1)):(d*h),g] # store the coefficients of the splines
						bmu_star[h,g]<-(b+(t(betamu[(1+d*(h-1)):(d*h),g])%*%Kj2%*%betamu[(1+d*(h-1)):(d*h),g])/2)[1,1] 
						taumu2[h,g]<-rgamma(1,a_star,bmu_star[h,g]) 
						taumustore[h,g,it]<-taumu2[h,g]
					}
					if(p>1){
						V<-1/(sigma2[g]+v)
						B<-V*(Y[which(S[,g]==1)]-t(BX[which(S[,g]==1),])%*%betamu[,g])*sigma2[g]
						L<-sqrt(V)
						T<-rnorm(1)
						gammamu[g]<-B+L*T
						gammamustore[g,it]<-gammamu[g]
					}
				}
				mu[,g]<-gammamu[g]+BX%*%betamu[,g]
				sigma2[g]<-rgamma(1,a+sum(S[,g])/2,b+0.5*sum((Y[which(S[,g]==1)]-gammamu[g]-BX[which(S[,g]==1),]%*%betamu[,g])^2)) 				
			}
			sigma2store[,it]<-sigma2 
			mustore[,,it]<-mu
			S<-matrix(0,N,G)
			for(i in 1:N){
				for(g in 1:G){
					W[i,g] <- exp(eta[i,g]+gamma[g])/sum(exp(eta[i,]+gamma))*dnorm(Y[i],mu[i,g],1/sqrt(sigma2[g]))
				}
				S[i,sample(1:G,1,prob=W[i,])]<-1
				llik[it]<-llik[it]+log(dnorm(as.numeric(Y[i]),mu[i,which(S[i,]==1)],1/sqrt(sigma2[which(S[i,]==1)])))
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
			#if(it==maxiter0){
			#	albero<-rpart(llik[1:maxiter]~c(1:maxiter))
			#	nodi<-albero$where
			#	b.new<-which.max(nodi==nodi[maxiter])
			#	if(b.new>maxiter*burn.in){
			#		diff.bin<-b.new-maxiter*burn.in
			#		maxiter<-maxiter+diff.bin
			#		burn.in<-b.new/maxiter			
			#	}
			#}
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
			betastore=betastore[,,1:maxiter],
			betamustore=betamustore[,,1:maxiter],
			gammastore=gammastore[,1:maxiter],
			gammamustore=gammamustore[,1:maxiter],
			taustore=taustore[,,1:maxiter],
			taumustore=taumustore[,,1:maxiter],
			sigma2store=sigma2store[,1:maxiter],
			llik=llik[1:maxiter],
			BX=BX,
			S=S,
			duration=duration
		)
		return(output)
	}
	if(G==1){
		v.n<-2.6/(diff(range(Y)))^2
		sigma2store<-rep(0,maxiter)
		mustore<-matrix(0,N,maxiter)
		sigma2<-var(Y)
		mu<-rnorm(N,mu0,sqrt(sigma2))
		betamustore<-matrix(0,d*p,maxiter) 
		betamu<-betamustore[,1]
		gammamustore<-rep(0,maxiter)  
		gammamu<-0
		taumustore<-matrix(0,p,maxiter)
		taumu2<-rgamma(p,shape=a,rate=b)
		a_star<-(a+rankMatrix(Kj)/2)[1]
		bmu_star<-rep(b,p)
		llik<-rep(0,maxiter)
		for (it in 1:maxiter){		
			for(h in 1:p){
				V<-chol2inv(chol(t(BX[,(1+d*(h-1)):(d*h)])%*%BX[,(1+d*(h-1)):(d*h)]*sigma2+Kj*taumu2[h]))
				if(p==1){
					B<-V%*%t(BX[,(1+d*(h-1)):(d*h)])%*%(Y-gammamu)*sigma2
				}
				else{
					B<-V%*%t(BX[,(1+d*(h-1)):(d*h)])%*%(Y-gammamu-
						BX[,-((1+d*(h-1)):(d*h))]%*%betamu[-((1+d*(h-1)):(d*h))])*sigma2
				}
				L<-chol(V) 
				T<-rnorm(d)
				betamu[(1+d*(h-1)):(d*h)]<-as.vector(B+L%*%T) 
				if(p>1){
					betamu[(1+d*(h-1)):(d*h)]<-as.vector(betamu[(1+d*(h-1)):(d*h)]-
						V%*%colSums(BX[,(1+d*(h-1)):(d*h)])%*%chol2inv(chol(
						t(colSums(BX[,(1+d*(h-1)):(d*h)]))%*%V%*%colSums(BX[,(1+d*(h-1)):(d*h)])))%*%
						colSums(BX[,(1+d*(h-1)):(d*h)])%*%betamu[(1+d*(h-1)):(d*h)])
				}
				betamustore[(1+d*(h-1)):(d*h),it]<-betamu[(1+d*(h-1)):(d*h)] 
				bmu_star[h]<-(b+(t(betamu[(1+d*(h-1)):(d*h)])%*%Kj%*%betamu[(1+d*(h-1)):(d*h)])/2)[1,1] 
				taumu2[h]<-rgamma(1,shape=a_star,rate=bmu_star[h]) 
				taumustore[h,it]<-taumu2[h]
			}
			if(p>1){
				V<-1/(N*sigma2+v)
				B<-V*sum(Y-BX%*%betamu)*sigma2
				L<-sqrt(V)
				T<-rnorm(1)
				gammamu<-B+L*T
				gammamustore[it]<-gammamu
			}
			mu<-gammamu+BX%*%betamu
			sigma2<-rgamma(1,a+N/2,b+0.5*sum((Y-gammamu-BX%*%betamu)^2)) 
			sigma2store[it]<-sigma2 
			mustore[,it]<-mu
			for(i in 1:N){
				llik[it]<-llik[it]+log(dnorm(as.numeric(Y[i]),mu[i],1/sqrt(sigma2)))
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
		finish = proc.time()[3]
		duration = finish - starttime
		output<-list(
			AICM=-2*(mean(llik[(burn.in+1):maxiter])-var(llik[(burn.in+1):maxiter])),
			betamustore=betamustore,
			gammamustore=gammamustore,
			taumustore=taumustore,
			BX=BX,
			sigma2store=sigma2store,
			llik=llik,
			duration=duration)
		return(output)
	}
}

full_l<-function (Y, X, G=3, delta = 1, a=1, b=0.005, a.n=1.28, b.n=0.36, v=0.01, burn.in=0.5, maxiter = 8000, seed=303) {
	cat("Start Time =",date(),"\n")
	flush.console()
	starttime = proc.time()[3] 
	set.seed(seed) 
	maxiter0<-maxiter
	X<-as.matrix(cbind(1,X)) 
	X<-orthonormalization(X,basis=F,norm=F)
	p<-ncol(X) 
	Y<-as.numeric(Y)
	N<-length(Y)
	mu0=mean(Y)
	if(G>1){
		S <- matrix(0,nrow = N, ncol = G) 
		Sstore<-array(0, dim=c(N,G,maxiter+maxiter-(burn.in*maxiter)))
		prstore<-matrix(0,G,maxiter+maxiter-(burn.in*maxiter))
		pr <- rdirichlet(1, rep(delta,G)) 
		v.n<-2.6/(diff(range(Y)))^2
		sigma2store<-matrix(0,G,maxiter+maxiter-(burn.in*maxiter))
		mustore<-array(0,dim=c(N,G,maxiter+maxiter-(burn.in*maxiter)))
		sigma2<-rep(var(Y)/G^2,G)
		mu<-matrix(rnorm(G,rep(mu0,G),sqrt(sigma2)),N,G,byrow=T)
		W<-matrix(0,N,G)
		for(i in 1:N) {
			S[i,sample(1:G,size=1,prob=pr*dnorm(Y[i],mu[i,],sqrt(sigma2)))]<-1	 
		}
		S<-S[,c(order(colSums(S))[-G],which.max(colSums(S)))]
		betamustore<-array(0,dim=c(p,G,maxiter+maxiter-(burn.in*maxiter))) 
		betamu<-betamustore[,,1]
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
		w<-rep(0,6) 
		betastore<-array(0,dim=c(p,G,maxiter+maxiter-(burn.in*maxiter))) 
		beta<-betastore[,,1]
		eta<-matrix(0,N,G) 
		C<-matrix(log(G-1),N,G-1)
		m<-matrix(0,N,G-1)
		llik<-rep(0,maxiter+maxiter-(burn.in*maxiter)) 
		it<-0
		while (it <= maxiter) {
			it<-it+1
			for(g in 1:(G-1)){
				if(G>2){
					C[,g]<-rowLogSumExps(eta[,-g])
				}
				V<-chol2inv(chol(t(X)%*%Lambda[,,g]%*%X+v*diag(p)))
				B<-V%*%t(X)%*%Lambda[,,g]%*%(Z[,g]+C[,g])
				L<-chol(V) 
				TT<-rnorm(p) 
				beta[,g]<-as.vector(B+L%*%TT) 
				betastore[,g,it]<-beta[,g] 
				eta[,g]<-X%*%beta[,g]
				m[,g]<-eta[,g]	
			}
			if(G>2){
				for(g in 1:(G-2)){
					C[,g]<-rowLogSumExps(eta[,-g])
				}
			}
			for(j in 1:N){
				for(g in 1:(G-1)){
					U<-runif(1)
					Z[j,g]<-log(exp(m[j,g]-C[j,g])*U+S[j,g])-log(1-U+exp(m[j,g]-C[j,g])*(1-S[j,g]))
					for(h in 1:6){
						w[h]<-omega[h]*dnorm(Z[j,g]+C[j,g],m[j,g],sqrt(s[h])) 
					}
					Lambda[j,j,g]<-1/s[sample(1:6,1,prob=w)]	
				}
			}
			for(g in 1:G){	
				if(sum(S[,g])>1){	
					V<-chol2inv(chol(t(X[which(S[,g]==1),])%*%X[which(S[,g]==1),]*sigma2[g]+v*diag(p)))
					B<-V%*%t(X[which(S[,g]==1),])%*%(Y[which(S[,g]==1)])*sigma2[g]
					L<-chol(V) 
					T<-rnorm(p)
					betamu[,g]<-as.vector(B+L%*%T)
					betamustore[,g,it]<-betamu[,g]
					mu[,g]<-X%*%betamu[,g]
					sigma2[g]<-rgamma(1,a+sum(S[,g])/2,b+0.5*sum((Y[which(S[,g]==1)]-X[which(S[,g]==1),]%*%betamu[,g])^2)) 
				}
				if(sum(S[,g])==1){
					V<-chol2inv(chol(X[which(S[,g]==1),]%*%t(X[which(S[,g]==1),])*sigma2[g]+v*diag(p)))
					B<-V%*%(X[which(S[,g]==1),])*(Y[which(S[,g]==1)])*sigma2[g]
					L<-chol(V) 
					T<-rnorm(p)
					betamu[,g]<-as.vector(B+L%*%T)
					betamustore[,g,it]<-betamu[,g] 
					mu[,g]<-X%*%betamu[,g]
					sigma2[g]<-rgamma(1,a+sum(S[,g])/2,b+0.5*sum((Y[which(S[,g]==1)]-X[which(S[,g]==1),]%*%betamu[,g])^2)) 				
				}

			}
			sigma2store[,it]<-sigma2 
			mustore[,,it]<-mu
			S<-matrix(0,N,G)
			for(i in 1:N){
				for(g in 1:G){
					W[i,g] <- exp(eta[i,g])/ 
							sum(exp(eta[i,]))*
							dnorm(Y[i],mu[i,g],1/sqrt(sigma2[g]))
					}
				S[i,sample(1:G,1,prob=W[i,])]<-1
				llik[it]<-llik[it]+log(dnorm(as.numeric(Y[i]),mu[i,which(S[i,]==1)],1/sqrt(sigma2[which(S[i,]==1)])))
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
			#if(it==maxiter0){
			#	albero<-rpart(llik[1:maxiter]~c(1:maxiter))
			#	nodi<-albero$where
			#	b.new<-which.max(nodi==nodi[maxiter])
			#	if(b.new>maxiter*burn.in){
			#		diff.bin<-b.new-maxiter*burn.in
			#		maxiter<-maxiter+diff.bin
			#		b.new<-max(b.new,(round(maxiter*burn.in)+1))
			#		burn.in<-b.new/maxiter			
			#	}
			#}
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
			betastore=betastore[,,1:maxiter],
			betamustore=betamustore[,,1:maxiter],
			sigma2store=sigma2store[,1:maxiter],
			X=X,
			llik=llik[1:maxiter],
			S=S,
			duration=duration
		)
		return(output)
	}
	if(G==1){
		v.n<-2.6/(diff(range(Y)))^2
		sigma2store<-rep(0,maxiter)
		mustore<-matrix(0,N,maxiter)			
		sigma2<-var(Y)
		mu<-rnorm(N,mu0,sqrt(sigma2))
		betamustore<-matrix(0,p,maxiter) 
		betamu<-betamustore[,1]
		llik<-rep(0,maxiter)
		for (it in 1:maxiter){		
			V<-chol2inv(chol(t(X)%*%X*sigma2+v*diag(p)))
			B<-V%*%t(X)%*%Y*sigma2
			L<-chol(V) 
			T<-rnorm(p)
			betamu<-as.vector(B+L%*%T) 
			betamustore[,it]<-betamu
			mu<-X%*%betamu
			sigma2<-rgamma(1,a+N/2,b+0.5*sum((Y-X%*%betamu)^2)) 
			sigma2store[it]<-sigma2 
			mustore[,it]<-mu
			for(i in 1:N){
				llik[it]<-llik[it]+log(dnorm(as.numeric(Y[i]),mu[i],1/sqrt(sigma2)))
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
		finish = proc.time()[3]
		duration = finish - starttime
		output<-list(
			AICM=-2*(mean(llik[(burn.in+1):maxiter])-var(llik[(burn.in+1):maxiter])),
			betamustore=betamustore,	
			sigma2store=sigma2store,
			X=X,
			llik=llik,
			duration=duration)
		return(output)
	}
}

mix_nl<-function (Y, X, G=3, delta = 1, a=1, b=0.005, a.n=1.28, b.n=0.36, v=0.01, nseg=20, bdeg=3, eps=0, burn.in=0.5, maxiter = 8000, seed=303) {
	cat("Start Time =",date(),"\n")
	flush.console()
	starttime = proc.time()[3] 
	set.seed(seed) 
	maxiter0<-maxiter
	X<-as.matrix(X) 
	p<-ncol(X) 
	Y<-as.numeric(Y)
	N<-length(Y)
	mu0=mean(Y)
	d<-nseg+bdeg  
	BX <- matrix(0,N,d*p) 
	for(h in 1:p){	
		BX[,(d*(h-1)+1):(h*d)]<-bspline(x=X[,h], xl=min(X[,h]), xr=max(X[,h]), ndx=nseg, bdeg=bdeg)
	}
	Kj = INLA:::inla.rw1(d) 
	Kj2<-Kj<-Kj+diag(eps,d)
	if(G>1){
		S <- matrix(0,nrow = N, ncol = G) 
		Sstore<-array(0, dim=c(N,G,maxiter+maxiter-(burn.in*maxiter)))
		prstore<-matrix(0,G,maxiter+maxiter-(burn.in*maxiter))
		pr <- rdirichlet(1, rep(delta,G)) 
		v.n<-2.6/(diff(range(Y)))^2
		sigma2store<-matrix(0,G,maxiter+maxiter-(burn.in*maxiter))
		mustore<-array(0,dim=c(N,G,maxiter+maxiter-(burn.in*maxiter)))
		sigma2<-rep(var(Y)/G^2,G)
		mu<-matrix(rnorm(G,rep(mu0,G),sqrt(sigma2)),N,G,byrow=T)
		W<-matrix(0,N,G)
		for(i in 1:N) {
			S[i,sample(1:G,size=1,prob=pr*dnorm(Y[i],mu[i,],sqrt(sigma2)))]<-1	 
		}
		S<-S[,c(order(colSums(S))[-G],which.max(colSums(S)))]
		betamustore<-array(0,dim=c(d*p,G,maxiter+maxiter-(burn.in*maxiter))) 
		betamu<-betamustore[,,1]
		gammamustore<-matrix(0,G,maxiter+maxiter-(burn.in*maxiter))  
		gammamu<-gammamustore[,1]
		taumustore<-array(0,c(p,G,maxiter+maxiter-(burn.in*maxiter)))
		taumu2<-matrix(rgamma(p*G,shape=a,rate=b),p,G)
		bmu_star<-matrix(b,p,G)
		a_star<-(a+rankMatrix(Kj)/2)[1]
		llik<-rep(0,maxiter+maxiter-(burn.in*maxiter)) 
		it<-0
		while (it <= maxiter) {
			it<-it+1
			pr<-rdirichlet(1,delta+colSums(S))
			prstore[,it]<-pr
			for(g in 1:G){	
				if(sum(S[,g])>1){	
					for(h in 1:p){
						V<-chol2inv(chol(t(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]*sigma2[g]+Kj2*taumu2[h,g]))
						if(p==1){
							B<-V%*%t(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%(Y[which(S[,g]==1)]-gammamu[g])*sigma2[g]
						}
						else{
							B<-V%*%t(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%(Y[which(S[,g]==1)]-gammamu[g]-
								BX[which(S[,g]==1),-((1+d*(h-1)):(d*h))]%*%betamu[-((1+d*(h-1)):(d*h)),g])*sigma2[g]
						}
						L<-chol(V) 
						T<-rnorm(d)
						betamu[(1+d*(h-1)):(d*h),g]<-as.vector(B+L%*%T)
						if(p>1){ 
							betamu[(1+d*(h-1)):(d*h),g]<-as.vector(betamu[(1+d*(h-1)):(d*h),g]-
								V%*%colSums(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%chol2inv(chol(
								t(colSums(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]))%*%V%*%colSums(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])))%*%
								colSums(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%betamu[(1+d*(h-1)):(d*h),g])
						}
						betamustore[(1+d*(h-1)):(d*h),g,it]<-betamu[(1+d*(h-1)):(d*h),g] 
						bmu_star[h,g]<-(b+(t(betamu[(1+d*(h-1)):(d*h),g])%*%Kj2%*%betamu[(1+d*(h-1)):(d*h),g])/2)[1,1] 
						taumu2[h,g]<-rgamma(1,shape=a_star,rate=bmu_star[h,g]) 
						taumustore[h,g,it]<-taumu2[h,g]
					}
					if(p>1){
						V<-1/(sum(S[,g])*sigma2[g]+v)
						B<-V*sum(Y[which(S[,g]==1)]-BX[which(S[,g]==1),]%*%betamu[,g])*sigma2[g]
						L<-sqrt(V)
						T<-rnorm(1)
						gammamu[g]<-B+L*T
						gammamustore[g,it]<-gammamu[g]
					}
				}
				if(sum(S[,g])==1){
					for(h in 1:p){
						V<-ginv(as.matrix(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]%*%t(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])*sigma2[g]+Kj2*taumu2[h,g]))
						if(p==1){
							B<-V%*%(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])*(Y[which(S[,g]==1)]-gammamu[g])*sigma2[g]
						}
						else{
							B<-V%*%(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%(Y[which(S[,g]==1)]-gammamu[g]-
								BX[which(S[,g]==1),-((1+d*(h-1)):(d*h))]%*%betamu[-((1+d*(h-1)):(d*h)),g])*sigma2[g]
						}
						L<-chol(V,pivot=TRUE) 
						T<-rnorm(d)
						betamu[(1+d*(h-1)):(d*h),g]<-as.vector(B+L%*%T)
						if(p>1){ 
							betamu[(1+d*(h-1)):(d*h),g]<-as.vector(betamu[(1+d*(h-1)):(d*h),g]-
								V%*%BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]%*%chol2inv(chol(
								t(BX[which(S[,g]==1),(1+d*(h-1)):(d*h)])%*%V%*%BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]))%*%
								BX[which(S[,g]==1),(1+d*(h-1)):(d*h)]%*%betamu[(1+d*(h-1)):(d*h),g])
						}
						betamustore[(1+d*(h-1)):(d*h),g,it]<-betamu[(1+d*(h-1)):(d*h),g] # store the coefficients of the splines
						bmu_star[h,g]<-(b+(t(betamu[(1+d*(h-1)):(d*h),g])%*%Kj2%*%betamu[(1+d*(h-1)):(d*h),g])/2)[1,1] 
						taumu2[h,g]<-rgamma(1,shape=a_star,rate=bmu_star[h,g]) 
						taumustore[h,g,it]<-taumu2[h,g]
					}
					if(p>1){
						V<-1/(sigma2[g]+v)
						B<-V*(Y[which(S[,g]==1)]-t(BX[which(S[,g]==1),])%*%betamu[,g])*sigma2[g]
						L<-sqrt(V)
						T<-rnorm(1)
						gammamu[g]<-B+L*T
						gammamustore[g,it]<-gammamu[g]
					}
				}
				mu[,g]<-gammamu[g]+BX%*%betamu[,g]
				sigma2[g]<-rgamma(1,a+sum(S[,g])/2,b+0.5*sum((Y[which(S[,g]==1)]-gammamu[g]-BX[which(S[,g]==1),]%*%betamu[,g])^2)) 				
			}
			sigma2store[,it]<-sigma2 
			mustore[,,it]<-mu
			S<-matrix(0,N,G)
			for(i in 1:N){
				for(g in 1:G){
					W[i,g] <- pr[g]*dnorm(Y[i],mu[i,g],1/sqrt(sigma2[g]))
				}
				S[i,sample(1:G,1,prob=W[i,])]<-1
				llik[it]<-llik[it]+log(dnorm(as.numeric(Y[i]),mu[i,which(S[i,]==1)],1/sqrt(sigma2[which(S[i,]==1)])))
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
			#if(it==maxiter0){
			#	albero<-rpart(llik[1:maxiter]~c(1:maxiter))
			#	nodi<-albero$where
			#	b.new<-which.max(nodi==nodi[maxiter])
			#	if(b.new>maxiter*burn.in){
			#		diff.bin<-b.new-maxiter*burn.in
			#		maxiter<-maxiter+diff.bin
			#		b.new<-max(b.new,(round(maxiter*burn.in)+1))
			#		burn.in<-b.new/maxiter			
			#	}
			#}
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
			prstore=prstore[,1:maxiter],
			betamustore=betamustore[,,1:maxiter],
			gammamustore=gammamustore[,1:maxiter],
			taumustore=taumustore[,,1:maxiter],
			sigma2store=sigma2store[,1:maxiter],
			llik=llik[1:maxiter],
			BX=BX,
			S=S,
			duration=duration
		)
		return(output)
	}
	if(G==1){
		v.n<-2.6/(diff(range(Y)))^2
		sigma2store<-rep(0,maxiter)
		mustore<-matrix(0,N,maxiter)
		sigma2<-var(Y)
		mu<-rnorm(N,mu0,sqrt(sigma2))
		betamustore<-matrix(0,d*p,maxiter) 
		betamu<-betamustore[,1]
		gammamustore<-rep(0,maxiter)  
		gammamu<-0
		taumustore<-matrix(0,p,maxiter)
		taumu2<-rgamma(p,shape=a,rate=b)
		a_star<-(a+rankMatrix(Kj)/2)[1]
		bmu_star<-rep(b,p)
		llik<-rep(0,maxiter)
		for (it in 1:maxiter){		
			for(h in 1:p){
				V<-chol2inv(chol(t(BX[,(1+d*(h-1)):(d*h)])%*%BX[,(1+d*(h-1)):(d*h)]*sigma2+Kj*taumu2[h]))
				if(p==1){
					B<-V%*%t(BX[,(1+d*(h-1)):(d*h)])%*%(Y-gammamu)*sigma2
				}
				else{
					B<-V%*%t(BX[,(1+d*(h-1)):(d*h)])%*%(Y-gammamu-
						BX[,-((1+d*(h-1)):(d*h))]%*%betamu[-((1+d*(h-1)):(d*h))])*sigma2
				}
				L<-chol(V) 
				T<-rnorm(d)
				betamu[(1+d*(h-1)):(d*h)]<-as.vector(B+L%*%T) 
				if(p>1){
					betamu[(1+d*(h-1)):(d*h)]<-as.vector(betamu[(1+d*(h-1)):(d*h)]-
						V%*%colSums(BX[,(1+d*(h-1)):(d*h)])%*%chol2inv(chol(
						t(colSums(BX[,(1+d*(h-1)):(d*h)]))%*%V%*%colSums(BX[,(1+d*(h-1)):(d*h)])))%*%
						colSums(BX[,(1+d*(h-1)):(d*h)])%*%betamu[(1+d*(h-1)):(d*h)])
				}
				betamustore[(1+d*(h-1)):(d*h),it]<-betamu[(1+d*(h-1)):(d*h)] # store the coefficients of the splines
				bmu_star[h]<-(b+(t(betamu[(1+d*(h-1)):(d*h)])%*%Kj%*%betamu[(1+d*(h-1)):(d*h)])/2)[1,1] 
				taumu2[h]<-rgamma(1,shape=a_star,rate=bmu_star[h]) 
				taumustore[h,it]<-taumu2[h]
			}
			if(p>1){
				V<-1/(N*sigma2+v)
				B<-V*sum(Y-BX%*%betamu)*sigma2
				L<-sqrt(V)
				T<-rnorm(1)
				gammamu<-B+L*T
				gammamustore[it]<-gammamu
			}
			mu<-gammamu+BX%*%betamu
			sigma2<-rgamma(1,a+N/2,b+0.5*sum((Y-gammamu-BX%*%betamu)^2)) 
			sigma2store[it]<-sigma2 
			mustore[,it]<-mu
			for(i in 1:N){
				llik[it]<-llik[it]+log(dnorm(as.numeric(Y[i]),mu[i],1/sqrt(sigma2)))
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
		finish = proc.time()[3]
		duration = finish - starttime
		output<-list(
			AICM=-2*(mean(llik[(burn.in+1):maxiter])-var(llik[(burn.in+1):maxiter])),
			betamustore=betamustore,
			gammamustore=gammamustore,
			taumustore=taumustore,
			BX=BX,
			sigma2store=sigma2store,
			llik=llik,
			duration=duration)
		return(output)
	}
}

mix_l<-function (Y, X, G=3, delta = 1, a=1, b=0.005, a.n=1.28, b.n=0.36, v=0.01, burn.in=0.5, maxiter = 8000, seed=303) {
	cat("Start Time =",date(),"\n")
	flush.console()
	starttime = proc.time()[3] 
	set.seed(seed) 
	maxiter0<-maxiter
	X<-as.matrix(cbind(1,X)) 
	p<-ncol(X)
	X<-orthonormalization(X,basis=F,norm=F)
	Y<-as.numeric(Y)
	N<-length(Y)
	mu0=mean(Y)
	if(G>1){
		S <- matrix(0,nrow = N, ncol = G) 
		Sstore<-array(0, dim=c(N,G,maxiter+maxiter-(burn.in*maxiter)))
		prstore<-matrix(0,G,maxiter+maxiter-(burn.in*maxiter))
		pr <- rdirichlet(1, rep(delta,G)) 
		v.n<-2.6/(diff(range(Y)))^2
		sigma2store<-matrix(0,G,maxiter+maxiter-(burn.in*maxiter))
		mustore<-array(0,dim=c(N,G,maxiter+maxiter-(burn.in*maxiter)))
		sigma2<-rep(var(Y)/G^2,G)
		mu<-matrix(rnorm(G,rep(mu0,G),sqrt(sigma2)),N,G,byrow=T)
		W<-matrix(0,N,G)
		for(i in 1:N) {
			S[i,sample(1:G,size=1,prob=pr*dnorm(Y[i],mu[i,],sqrt(sigma2)))]<-1	 
		}
		S<-S[,c(order(colSums(S))[-G],which.max(colSums(S)))]
		betamustore<-array(0,dim=c(p,G,maxiter+maxiter-(burn.in*maxiter))) 
		betamu<-betamustore[,,1]		
 		llik<-rep(0,maxiter+maxiter-(burn.in*maxiter)) 
		it<-0
		while (it <= maxiter) {
			it<-it+1
			pr<-rdirichlet(1,delta+colSums(S))
			prstore[,it]<-pr
			for(g in 1:G){	
				if(sum(S[,g])>1){	
					V<-chol2inv(chol(t(X[which(S[,g]==1),])%*%X[which(S[,g]==1),]*sigma2[g]+v*diag(p)))
					B<-V%*%t(X[which(S[,g]==1),])%*%(Y[which(S[,g]==1)])*sigma2[g]
					L<-chol(V) 
					T<-rnorm(p)
					betamu[,g]<-as.vector(B+L%*%T)
					betamustore[,g,it]<-betamu[,g]
					mu[,g]<-X%*%betamu[,g]
					sigma2[g]<-rgamma(1,a+sum(S[,g])/2,b+0.5*sum((Y[which(S[,g]==1)]-X[which(S[,g]==1),]%*%betamu[,g])^2)) 
				}
				if(sum(S[,g])==1){
					V<-chol2inv(chol(X[which(S[,g]==1),]%*%t(X[which(S[,g]==1),])*sigma2[g]+v*diag(p)))
					B<-V%*%(X[which(S[,g]==1),])*(Y[which(S[,g]==1)])*sigma2[g]
					L<-chol(V) 
					T<-rnorm(p)
					betamu[,g]<-as.vector(B+L%*%T)
					betamustore[,g,it]<-betamu[,g] 
					mu[,g]<-X%*%betamu[,g]
					sigma2[g]<-rgamma(1,a+sum(S[,g])/2,b+0.5*sum((Y[which(S[,g]==1)]-X[which(S[,g]==1),]%*%betamu[,g])^2)) 				
				}
			}
			sigma2store[,it]<-sigma2 
			mustore[,,it]<-mu
			S<-matrix(0,N,G)
			for(i in 1:N){
				for(g in 1:G){
					W[i,g] <- pr[g]*dnorm(Y[i],mu[i,g],1/sqrt(sigma2[g]))
				}
				S[i,sample(1:G,1,prob=W[i,])]<-1
				llik[it]<-llik[it]+log(dnorm(as.numeric(Y[i]),mu[i,which(S[i,]==1)],1/sqrt(sigma2[which(S[i,]==1)])))
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
			#if(it==maxiter0){
			#	albero<-rpart(llik[1:maxiter]~c(1:maxiter))
			#	nodi<-albero$where
			#	b.new<-which.max(nodi==nodi[maxiter])
			#	if(b.new>maxiter*burn.in){
			#		diff.bin<-b.new-maxiter*burn.in
			#		maxiter<-maxiter+diff.bin
			#		b.new<-max(b.new,(round(maxiter*burn.in)+1))
			#		burn.in<-b.new/maxiter			
			#	}
			#}
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
			prstore=prstore[,1:maxiter],
			betamustore=betamustore[,,1:maxiter],
			sigma2store=sigma2store[,1:maxiter],
			mustore=mustore[,,1:maxiter],
			X=X,
			llik=llik[1:maxiter],
			S=S,
			duration=duration
		)
		return(output)
	}
	if(G==1){
		v.n<-2.6/(diff(range(Y)))^2
		sigma2store<-rep(0,maxiter)
		mustore<-matrix(0,N,maxiter)			
		sigma2<-var(Y)
		mu<-rnorm(N,mu0,sqrt(sigma2))
		betamustore<-matrix(0,p,maxiter) 
		betamu<-betamustore[,1]
		llik<-rep(0,maxiter)
		for (it in 1:maxiter){		
			V<-chol2inv(chol(t(X)%*%X*sigma2+v*diag(p)))
			B<-V%*%t(X)%*%Y*sigma2
			L<-chol(V) 
			T<-rnorm(p)
			betamu<-as.vector(B+L%*%T) 
			betamustore[,it]<-betamu
			mu<-X%*%betamu
			sigma2<-rgamma(1,a+N/2,b+0.5*sum((Y-X%*%betamu)^2)) 
			sigma2store[it]<-sigma2 
			mustore[,it]<-mu
			for(i in 1:N){
				llik[it]<-llik[it]+log(dnorm(as.numeric(Y[i]),mu[i],1/sqrt(sigma2)))
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
		finish = proc.time()[3]
		duration = finish - starttime
		output<-list(
			AICM=-2*(mean(llik[(burn.in+1):maxiter])-var(llik[(burn.in+1):maxiter])),
			betamustore=betamustore,	
			sigma2store=sigma2store,
			X=X,
			llik=llik,
			duration=duration)
		return(output)
	}
}
