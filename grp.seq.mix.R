# function to get arm size for sar statistic
# packages: rmutil, mvtnorm
grp.seq.mix <- function(num.stages=2,alpha,beta,rho,theta,K=1,mu=0,sigma=1,f="normal",lambda=1,
                            B=10000,method="sar",norm.approx=TRUE){
  m.start = 1
  if(method=="sr"){method<-"rerank"}
  f <- ifelse(f=="t3","t",f)
  k <- num.stages
  ratio <- lambda
  std <- TRUE
  pi1 <- c(alpha*(1/k)^rho,diff(alpha*((1:k)/k)^rho))
  pi2 <- c(beta*(1/k)^rho,diff(beta*((1:k)/k)^rho))
  d <- f
  arg1 <- mu; arg2 <- sigma
  if(f=="normal"){
    f <- function(n,arg1=0,arg2=1){rnorm(n=n,mean=arg1,sd=arg2)}
    ppsi <- function(x,arg1=0,arg2=1){pnorm(x,mean=arg1,sd=arg2)}
    dpsi <- function(x,arg1=0,arg2=1,lg=FALSE){dnorm(x,mean=arg1,sd=arg2,log=lg)}
    qpsi <- function(q,arg1=0,arg2=1) qnorm(q)*arg2 + arg1
  } else if(f=="logistic"){
    f <- function(n,arg1=0,arg2){rlogis(n=n,location=arg1,arg2*sqrt(3)/pi)}
    ppsi <- function(x,arg1=0,arg2=1){plogis(x,location=arg1,scale=arg2*sqrt(3)/pi)}
    dpsi <- function(x,arg1=0,arg2=1,lg=FALSE){
      dlogis(x,location=arg1,scale=arg2*sqrt(3)/pi,log=lg)}
    qpsi <- function(q,arg1=0,arg2=1){ qlogis(q,location=arg1,scale=arg2*sqrt(3)/pi) }
  } else if(f=="laplace"){
    # uses package rmutil
    f <- function(n,arg1=0,arg2){rmutil::rlaplace(n=n,m=arg1,s=arg2/sqrt(2))}
    ppsi <- function(x,arg1=0,arg2=1){rmutil::plaplace(x,m=arg1,s=arg2/sqrt(2))}
    dpsi <- function(x,arg1=0,arg2=1,lg=FALSE){
      rmutil::dlaplace(x,m=arg1,s=arg2/sqrt(2),log=lg)}
    qpsi <- function(q,arg1=0,arg2=1){ rmutil::qlaplace(q,m=arg1,s=arg2/sqrt(2)) }
  } else if(f=="t"){
    f <- function(n,arg1=0,arg2){arg1 + arg2*rt(n=n,df=3,ncp=0)/sqrt(3)}
    ppsi <- function(x,arg1=0,arg2=1){pt((x-arg1)*sqrt(3)/arg2,df=3,ncp=0)}
    dpsi <- function(x,arg1=0,arg2=1,lg=FALSE){
      if(lg==FALSE){ return(dt((x-arg1)*sqrt(3)/arg2,df=3,ncp=0)*sqrt(3)/arg2)
      } else { dt((x-arg1)*sqrt(3)/arg2,df=3,ncp=0,log=T)+log(sqrt(3))-log(arg2) }
    }
    qpsi <- function(q,arg1=0,arg2=1) qt(q,df=3)/sqrt(3)*arg2 + arg1
  }
  # normal approx for one stage
  delta <- K*sigma
  h1<-function(y,theta,K) {ppsi(y)*((1-theta)*dpsi(y) + theta*dpsi(y-K))}
  gamval<-integrate(h1,-Inf,Inf,theta,K)$value
  h2<-function(y,theta,K) {((1-(1-theta)*ppsi(y)-theta*ppsi(y-K))^2)*dpsi(y)} 
  xi1<-integrate(h2,-Inf,Inf,theta,K)$value-gamval^2
  h3<-function(y,theta,K) {(ppsi(y))^2*((1-theta)*dpsi(y) + theta*dpsi(y-K))}
  xi2<-integrate(h3,-Inf,Inf,theta,K)$value-gamval^2
  lambda <- 1/2
  # P(X<Y)
  pxy <- integrate(function(y){
    ppsi(y)*((1-theta)*dpsi(y)+theta*dpsi(y-K))},-Inf,Inf)$value
  # P(X1 < Y1, X1 < Y2)
  pxyy <- integrate(function(y){
    (1-(1-theta)*ppsi(y)-theta*ppsi(y-K))^2*dpsi(y)},-Inf,Inf)$value
  # P(X1 < Y1, X2 < Y2)
  pxxy <- integrate(function(y){
    (ppsi(y))^2*((1-theta)*dpsi(y)+theta*dpsi(y-K))},-Inf,Inf)$value
  # jeske yao formula for fixed sample arm size equal arm sizes
  m.approx <- ((qnorm(1-alpha)/sqrt(6) + qnorm(1-beta)*sqrt(xi1+xi2))/(gamval-1/2))^2
  # jeske yao formula unequal arm sizes
  # ratio <- 0.8
  m.approx2 <- (((qnorm(1-alpha)*sqrt((ratio+1)/(12*ratio))) + qnorm(1-beta)*sqrt(xi1+xi2/ratio))/(gamval-1/2))^2
  
  if(norm.approx==FALSE){
    if(k==1){ # simulate just one stage
      delta <- K*sigma
      p <- 1
      if((m.approx-10)>m.start){m.start <- floor(m.approx - 10)}
      inc <- function(incr = 1, m = m.start){
        while(p>=pi2){
          ifelse(m<15,m<-m+1,m <- m+incr)
          n <- ceiling(ratio*m)
          w <- vector("numeric",B)
          w.alt <- vector("numeric",B)
          for(i in 1:B){
            # obs <- matrix(f(n=k*2*m,arg2=arg2),ncol=2,nrow=(k*m),byrow=F)
            x.null <- f(n=k*m,arg2=arg2); y.null <- f(n=k*n,arg2=arg2)
            # obs.alt <- NULL
            # y <- sapply(runif(m),FUN = function(i){
            #   if(i<theta){f(n=1,arg1=delta,arg2=arg2)}else{f(n=1,arg2=arg2)}})
            # obs.alt <- rbind(obs.alt,matrix(c(f(n=m,arg2=sigma),y),ncol=2,nrow=m,byrow=F))
            x.alt <- f(n=m,arg2=sigma)
            y.alt <- sapply(runif(n),FUN = function(i){
              if(i<theta){f(n=1,arg1=delta,arg2=arg2)}else{f(n=1,arg2=arg2)}})
            # calculate rank sum value
            # w[i] <- sum(rank(obs)[(m+1):(m+m)])
            # w.alt[i] <- sum(rank(obs.alt)[(m+1):(m+m)])
            w[i] <- sum(rank(c(x.null,y.null))[(m+1):(m+n)])
            w.alt[i] <- sum(rank(c(x.alt,y.alt))[(m+1):(m+n)])
          }
          
          p.r <- 0
          # index.max <- m*m + m*(m+1)/2
          index.max <- m*n + n*(n+1)/2
          u <- sort(unique(w),decreasing = TRUE)
          index <- index.max + 1
          while(p.r <= pi1){
            p.rej <- p.r
            index <- index - 1
            p.r <- sum(w>=(index))/B
            if(index==1){break}
          }
          if(index==index.max & sum(w>=index)/B>pi1[1]){
            r <- Inf
          } else if(sum(w>=index)/B <= pi1[1]) {
            r <- index
          } else {
            r <- index+1
          }
          if(std==TRUE){mu <- n*((m+n)+1)/2; sig <- sqrt(m*n*((m+n)+1)/12)
          r <- (r-mu)/sig
          w.alt <- (w.alt - mu)/sig; w <- (w-mu)/sig}
          
          p <- sum(w.alt<r)/B
        }
        p.rej <- sum(w>=r)/B
        p.acc <- sum(w.alt<r)/B
        # out <- list(r=r,m=m,theta=theta,K=K,B=B,distribution=d,p.rej=p.rej,p.acc=p.acc)
        out <- list("Design Alternative"=data.frame(Stages=k,f=d,theta=theta,K=K,row.names = ""),
                    "Arm Size"=data.frame(m=m,n=n,row.names = ""),
                    "Critical Value"=round(r,4),
                    "Type I Error Rate"=round(p.rej,4),
                    "Type II Error Rate"=round(p.acc,4))
        return(out)
      }
      incr.size <- 10
      res <- inc(incr=incr.size)
      if(res$`Arm Size`$m>=15){
        res <- inc(incr = 1, m = res$`Arm Size`$m-incr.size-1)
        return(res)
      } else {return(res)}
    } else { # simulate 2+ stages
      delta <- K*sigma
      p.start <- 1
      # use normal approx for one stage to get close initially
      if((m.approx/k-10)>m.start){m.start <- floor(m.approx/k - 10)}
      incr.size <- 10
      incr.arm.size <- function(incr = 1, p = 1, m = m.start){
        while(p>=pi2[k]){
          p.previous <- p
          ifelse(m<15,m<-m+1,m <- m+incr)
          n <- ceiling(ratio*m)
          r <- vector("numeric",k)
          a <- vector("numeric",k)
          if(method=="rerank"){w <- matrix(0,ncol = k,nrow=B)
            w.alt <- matrix(0,ncol = k,nrow=B)
          }else{sar <- matrix(0,ncol = k,nrow=B); sar.alt <- matrix(0,ncol = k,nrow=B)}
          
          ### simulate the observations and obtain ranks
          for(i in 1:B){
            # obs <- matrix(f(n=k*2*m,arg2=sigma),ncol=2,nrow=(k*m),byrow=F)
            # ind <- sample(1:2,prob=c(1-theta,theta),size=k*m,replace = TRUE)
            # mus <- c(0,delta)
            # y <- f(n=k*m,arg1=mus[ind],arg2=sigma)
            # obs.alt <- matrix(c(f(n=k*m,arg2=sigma),y),ncol=2,nrow=(k*m),byrow=F)
            # for(i1 in 1:k){
            #   if(method=="rerank"){
            #     w[i,i1] <- sum(rank(obs[1:(i1*m),])[(i1*m+1):(i1*m+i1*m)])
            #     w.alt[i,i1] <- sum(rank(obs.alt[1:(i1*m),])[(i1*m+1):(i1*m+i1*m)])
            #   }
            #   else{
            #     sar[i,i1] <- sum(rank(obs[((i1-1)*m+1):(i1*m),])[(n+1):(n+m)])
            #     sar.alt[i,i1] <- sum(rank(obs.alt[((i1-1)*m+1):(i1*m),])[(n+1):(n+m)])
            #   }
            # }
            x.null <- matrix(f(k*m,arg2=arg2),ncol=k,nrow=m); y.null <- matrix(f(k*n,arg2=arg2),ncol=k,nrow=n)
            x.alt <- matrix(f(k*m,arg2=arg2),ncol=k,nrow=m)
            ind <- sample(1:2,prob=c(1-theta,theta),size=k*n,replace = TRUE)
            mus <- c(0,delta)
            y.alt <- matrix(f(n=k*n,arg1=mus[ind],arg2=sigma),ncol=k,nrow=n)
            if(method=="rerank"){
              w[i,] <- sapply(1:k,FUN=function(z){sum(rank(c(x.null[,1:z],y.null[,1:z]))[(z*m+1):(z*(m+n))])})
              w.alt[i,] <- sapply(1:k,FUN=function(z){sum(rank(c(x.alt[,1:z],y.alt[,1:z]))[(z*m+1):(z*(m+n))])})
            } else {
              sar[i,] <- sapply(1:k,FUN=function(z){sum(rank(c(x.null[,z],y.null[,z]))[(m+1):(m+n)])})
              sar.alt[i,] <- sapply(1:k,FUN=function(z){sum(rank(c(x.alt[,z],y.alt[,z]))[(m+1):(m+n)])})
            }
          }
          # ts = test statistic
          if(method=="rerank" & std==TRUE){
            # mu <- (1:k)*m*((1:k)*(m+m)+1)/2
            # sig <- sqrt((1:k)^2*m*m*((1:k)*(m+m)+1)/12)
            mu <- (1:k)*n*((1:k)*m*0.5 + ((1:k)*n+1)/2)
            sig <- sqrt((1:k)^2*m*n*(1/2 - 0.5^2*((1:k)*(m+n)-1) + ((1:k)*n-1)*1/3 + ((1:k)*m-1)*1/3))
            ts <- t(apply(w,1,FUN=function(x) (x - mu)/(sig)))
            ts.alt <- t(apply(w.alt,1,FUN=function(x) (x - mu)/(sig)))
          }
          else if(method=="rerank" & std==FALSE){
            ts <- w; ts.alt <- w.alt
          }
          if(method!="rerank" & std==TRUE){ # sar standardized
            # mu <- m*((m+m)+1)/2
            # sig <- sqrt(m*m*((m+m)+1)/12)
            mu <- n*(m*0.5 + (n+1)/2)
            sig <- sqrt(m*n*(1/2 - 0.5^2*((m+n)-1) + (n-1)*1/3 + (m-1)*1/3))
            ts <- t(apply(sar,1,FUN=function(x){
              (cumsum(x)/seq_along(x) - mu)/(sig/sqrt(seq_along(x)))}))
            ts.alt <- t(apply(sar.alt,1,FUN=function(x){
              (cumsum(x)/seq_along(x) - mu)/(sig/sqrt(seq_along(x)))}))
          }
          if(method!="rerank" & std==FALSE){ # sar not standardized
            ts <- t(apply(sar,1,FUN=function(x) (cumsum(x)/seq_along(x))))
            ts.alt <- t(apply(sar.alt,1,FUN=function(x) (cumsum(x)/seq_along(x))))
          }
          
          p.r <- vector("numeric",length=k)
          p.a <- vector("numeric",length=k)
          p <- 0
          index <- 0
          u <- sort(unique(ts[,1]),decreasing = T)
          while(p <= pi1[1]){
            p.rej <- p
            index <- index + 1
            p <- sum(ts[,1]>=u[index])/B
            if(index==length(u)){break}
          }
          if(index==1 & sum(ts[,1]>=u[index])/B>pi1[1]){
            r[1] <- Inf
          } else if(sum(ts[,1]>=u[index])/B <= pi1[1]) {
            r[1] <- u[index]
          } else {
            r[1] <- u[index-1]
          }
          p.r[1] <- sum(ts[,1]>=r[1])/B
          
          p.acc <- 0
          u.alt <- sort(unique(ts.alt[,1]))
          index <- 0
          while(p.acc <= pi2[1]){
            index <- index + 1
            p.acc <- sum(ts.alt[,1]<=u.alt[index])/B
          }
          if(index==1 & sum(ts.alt[,1]<=u.alt[index])/B>pi2[1]){
            a[1] <- -Inf
          } else if(sum(ts.alt[,1]<=u.alt[index])/B <= pi2[1]) {
            a[1] <- u.alt[index]
          } else {
            a[1] <- u.alt[index-1]
          }
          p.a[1] <- sum(ts.alt[,1]<=a[1])/B
          
          s <- ts[ts[,1]<r[1] & ts[,1]>a[1],]
          s.alt <- ts.alt[ts.alt[,1]<r[1] & ts.alt[,1]>a[1],]
          for(i in 2:k){
            if(dim(s)[1]==0 || dim(s.alt)[1]==0){break}
            p.rej <- 0
            u <- sort(unique(s[,i]),decreasing = TRUE)
            index <- 0
            while(p.rej <= pi1[i]){
              index <- index + 1
              p.rej <- sum(s[,i]>=(u[index]))/B
              if(index==length(u)){break}
            }
            if(index==1 & sum(s[,i]>=u[index])/B>pi1[i]){
              r[i] <- Inf
            } else if(sum(s[,i]>=u[index])/B <= pi1[i]) {
              r[i] <- u[index]
            } else {
              r[i] <- u[index-1]
            }
            p.r[i] <- sum(s[,i]>=r[i])/B
            
            p.acc <- 0
            u.alt <- sort(unique(s.alt[,i]))
            index <- 0
            while(p.acc <= pi2[i]){
              index <- index + 1
              p.acc <- sum(s.alt[,i]<=u.alt[index])/B
              if(index==length(u.alt)){break}
            }
            if(index==1 & sum(s.alt[,i]<=u.alt[index])/B>pi2[i]){
              a[i] <- -Inf
            } else if(sum(s.alt[,i]<=u.alt[index])/B <= pi2[i]) {
              a[i] <- u.alt[index]
            } else {
              a[i] <- u.alt[index-1]
            }
            p.a[i] <- sum(s.alt[,i]<=a[i])/B
            if(i!=k){
              s <- s[s[,i]<r[i] & s[,i]>a[i],]
              s.alt <- s.alt[s.alt[,i]<r[i] & s.alt[,i]>a[i],]
            }
          }
          p <- sum(s.alt[,k]<r[k])/B
          p.out <- p
        }
        if(dim(s)[1]==0 || dim(s.alt)[1]==0){return("Unable to achieve design.")
        } else{
          # return(list(r=r,a=a,m=m,p.r=p.r,p.a=p.a,p.out=p.out,K=K,theta=theta))
          out <- list("Design Alternative"=data.frame(Stages=k,f=d,theta=theta,K=K,row.names = ""),
                      "Arm Size"=data.frame(m=m,n=n),
                      "Upper Critical Values"=round(r,4),
                      "Lower Critical Values"=round(c(a[1:(k-1)],r[k]),4),
                      "Type I Error Rate"=round(c(p.r,sum(p.r)),4),
                      "Type II Error Rate"=round(c(c(p.a[1:(k-1)],p.out),
                                              sum(c(p.a[1:(k-1)],p.out))),4))
          names(out$`Type I Error Rate`) <- c(paste("Stage",1:k),"Overall")
          names(out$`Type II Error Rate`) <- c(paste("Stage",1:k),"Overall")
          return(out)
        }
      }
      res1 <- incr.arm.size(incr=incr.size,p=p.start,m=m.start)
      if(is.list(res1)==FALSE){return(res1)
      } else if(res1$`Arm Size`$m>=15){
        # res2 <- incr.arm.size(incr=1,p=1, m = res1$m-incr.size-1)
        res2 <- incr.arm.size(incr=1,p=1, m = res1$`Arm Size`$m-incr.size-1)
        return(res2)
      } else{return(res1)}
    }
  } else { # norm.approx == TRUE
    if(k==1){ # fixed sample
      delta <- K*sigma
      p <- 1
      if((m.approx-10)>m.start){m <- floor(m.approx - 10)}else{m <- 0}
      while(p>=pi2){
        m <- m+1
        n <- ceiling(ratio*m)
        mu <- n*(m+n+1)/2
        v <- m*n*(m+n+1)/12
        r <- qnorm(1-pi1)
        # p <- pnorm(r,mean=(n*(m*pxy + (n+1)/2) - mu)/sqrt(v),
        #            sd=m^2*(pxy*(1-pxy)+(m-1)*(pxyy-pxy^2+pxxy-pxy^2))/v)
        p <- pnorm(r,mean=(n*(m*pxy + (n+1)/2) - mu)/sqrt(v),
                   sd=sqrt((m*n*pxy*(1-pxy) + m*n*(n-1)*(pxyy - pxy^2) + m*(m-1)*n*(pxxy-pxy^2)))/sqrt(v))
        # w = u + n*(n+1)/2
        # u = w - n*(n+1)/2
      }
      
      # return(list(r=r,m=m,p=p,K=K,theta=theta))
      out <- list("Design Alternative"=data.frame(Stages=k,f=d,theta=theta,K=K,row.names = ""),
                  "Arm Size"=data.frame(m=m,n=n,row.names = ""),
                  "Critical Value"=round(r,4),
                  "Type I Error Rate"=round(pnorm(r,lower.tail = FALSE),4),
                  "Type II Error Rate"=round(p,4))
      return(out)
      
    }else{ # multiple stages
      delta <- K*sigma
      p <- 1
      # if((m.approx/k-10)>m.start){m <- floor(m.approx/k - 10)}else{m <- 0}
      m <- 0
      while(p>=pi2[k]){
        m <- m+1
        n <- ceiling(ratio*m)
        r <- vector("numeric",k)
        a <- vector("numeric",k)
        
        if(method=="sar" & std=="TRUE"){
          # mu <- m*(m+m+1)/2
          # v <- m*m*(m+m+1)/12
          mu <- n*(m+n+1)/2
          v <- m*n*(m+n+1)/12
          null.mu <- rep(0,k)
          null.var <- v/(1:k)
          null.cov <- outer(1:k,1:k,Vectorize(function(x,y){sqrt(x*y)/max(x,y)}))
          # alt.mu <- (m*(m*gamval + (m+1)/2) - mu)/(sqrt(v)/sqrt(1:k))
          # alt.v <- m^2*(pxy*(1-pxy)+(m-1)*(pxyy-pxy^2+pxxy-pxy^2))/v
          # alt.cov <- outer(1:k,1:k,Vectorize(function(x,y){sqrt(x*y)/max(x,y)}))*alt.v
          alt.mu <- (m*n*pxy + n*(n+1)/2 - mu)/(sqrt(v)/sqrt(1:k))
          alt.v <- (m*n*pxy*(1-pxy) + m*n*(n-1)*(pxyy - pxy^2) + m*(m-1)*n*(pxxy-pxy^2))/v
          alt.cov <- outer(1:k,1:k,Vectorize(function(x,y){sqrt(x*y)/max(x,y)}))*alt.v
        }
        if(method=="rerank" & std=="TRUE"){
          # mu <- (1:k)*m*((1:k)*(m+m)+1)/2
          # v <- (1:k)^2*m*m*((1:k)*(m+m)+1)/12
          mu <- (1:k)*n*((1:k)*(m+n)+1)/2
          v <- (1:k)^2*m*n*((1:k)*(m+n)+1)/12
          null.mu <- rep(0,k)
          null.var <- rep(1,k)
          null.cov <- outer(1:k,1:k,Vectorize(function(x,y){
            (((min(x,y)*m)^2*(2*max(x,y)*m+1))/12)/(sqrt(v[x]*v[y]))}))
          null.cov <- outer(1:k,1:k,Vectorize(function(x,y){
            s <- min(x,y)
            sp <- max(x,y)
            ((s^2*m*n)*(1/2 - (sp*(m+n)-1)*0.5^2 + (sp*n-1)*1/3 + (sp*m-1)*1/3))/sqrt(v[x]*v[y])
          }))
          # alt.mu <- ((1:k)*m*((1:k)*m*gamval + ((1:k)*m+1)/2) - mu)/sqrt(v)
          # alt.var <- (((1:k)*m)^2*(pxy*(1-pxy) + ((1:k)*m-1)*(pxyy-pxy^2) + 
          #                            ((1:k)*m-1)*(pxxy-pxy^2)))/v
          # alt.cov <- outer(1:k,1:k,Vectorize(function(x,y){
          #   ((min(x,y)*m)^2*
          #      (pxy*(1-pxy) + (max(x,y)*m-1)*(pxyy-pxy^2 + pxxy-pxy^2)))/sqrt(v[x]*v[y])
          # }))
          alt.mu <- ((1:k)*n*((1:k)*m*pxy + ((1:k)*n+1)/2) - mu)/sqrt(v)
          alt.var <- (((1:k)*m)*((1:k)*n)*(pxy - pxy^2*(((1:k)*m)+((1:k)*n)-1) + (((1:k)*n)-1)*pxyy + (((1:k)*m)-1)*pxxy))/v
          alt.cov <- outer(1:k,1:k,Vectorize(function(x,y){
            s <- min(x,y)
            sp <- max(x,y)
            # sig.s <- sqrt((s)^2*m*n*(1/2 - 0.5^2*((s)*(m+n)-1) + ((s)*n-1)*1/3 + ((s)*m-1)*1/3))
            # sig.sp <- sqrt((sp)^2*m*n*(1/2 - 0.5^2*((sp)*(m+n)-1) + ((sp)*n-1)*1/3 + ((sp)*m-1)*1/3))
            # (1/(sig.s*sig.sp))*((s^2*m*n)*(pxy - (sp*(m+n)-1)*pxy^2 + (sp*n-1)*pxyy + (sp*m-1)*pxxy))
            ((s^2*m*n)*(pxy - (sp*(m+n)-1)*pxy^2 + (sp*n-1)*pxyy + (sp*m-1)*pxxy))/sqrt(v[x]*v[y])
          }))
        }
        
        
        
        # critical values
        r[1] <- qnorm(1-pi1[1],mean = null.mu[1],sd=sqrt(null.cov[1,1]))
        a[1] <- qnorm(pi2[1],mean=alt.mu[1],sd=sqrt(alt.cov[1,1]))
        for(i in 2:k){ # uses mvtnorm package
          if(m>3){
            r[i] <- uniroot(f=function(x){
              as.numeric(mvtnorm::pmvnorm(lower=c(a[1:(i-1)],x),
                                          upper=c(r[1:(i-1)],Inf),
                                          mean=null.mu[1:i],
                                          sigma=null.cov[1:i,1:i])) - pi1[i]},
                            lower=0,upper=4)$root
          } else {
            r[i] <- optimize(f=function(x){
              (as.numeric(mvtnorm::pmvnorm(lower=c(a[1:(i-1)],x),
                                           upper=c(r[1:(i-1)],Inf),
                                           mean=null.mu[1:i],
                                           sigma=null.cov[1:i,1:i])) - pi1[i])^2},
                             lower=0,upper=4)$minimum
          }
          if(i!=k){
            if(m>3){
              a[i] <- uniroot(f=function(x){
                as.numeric(mvtnorm::pmvnorm(lower=c(a[1:(i-1)],-Inf),
                                            upper=c(r[1:(i-1)],x),
                                            mean=alt.mu[1:i],
                                            sigma=alt.cov[1:i,1:i])) - pi2[i]},
                              lower=-4,upper=4)$root
            } else {
              a[i] <- optimize(f=function(x){
                (as.numeric(mvtnorm::pmvnorm(lower=c(a[1:(i-1)],-Inf),
                                             upper=c(r[1:(i-1)],x),
                                             mean=alt.mu[1:i],
                                             sigma=alt.cov[1:i,1:i])) - pi2[i])^2},
                               lower=-4,upper=4)$minimum
            }
          }
          rm(i)
        }
        # uses mvtnorm package
        p <- as.numeric(mvtnorm::pmvnorm(lower=c(a[1:(k-1)],-Inf),
                                         upper=c(r[1:(k-1)],r[k]),
                                         mean=alt.mu[1:k],
                                         sigma=alt.cov[1:k,1:k]))
        
        
        p.r <- p.a <- vector("numeric",k)
        p.r[1] <- pnorm(r[1],mean=null.mu[1],sd=sqrt(null.cov[1,1]),lower.tail = FALSE)
        p.a[1] <- pnorm(a[1],mean=alt.mu[1],sd=sqrt(alt.cov[1,1]))
        for(i in 2:k){
          if(i!=k){ # uses mvtnorm package
            p.r[i] <- mvtnorm::pmvnorm(lower=c(a[1:(i-1)],r[i]),upper=c(r[1:(i-1)],Inf),
                                       mean=null.mu[1:i],sigma=null.cov[1:i,1:i])
            p.a[i] <- mvtnorm::pmvnorm(lower=c(a[1:(i-1)],-Inf),upper=c(r[1:(i-1)],a[i]),
                                       mean=alt.mu[1:i],sigma=alt.cov[1:i,1:i])
          } else {
            p.r[i] <- mvtnorm::pmvnorm(lower=c(a[1:(i-1)],r[i]),upper=c(r[1:(i-1)],Inf),
                                       mean=null.mu[1:i],sigma=null.cov[1:i,1:i])
            p.a[i] <- mvtnorm::pmvnorm(lower=c(a[1:(i-1)],-Inf),upper=c(r[1:(i-1)],r[i]),
                                       mean=alt.mu[1:i],sigma=alt.cov[1:i,1:i])
          }
          rm(i)
        }
      }
      # return(list(r=r,a=c(a[1:(k-1)],r[k]),m=m,p.r=p.r,p.a=p.a,K=K,theta=theta))
      out <- list("Design Alternative"=data.frame(Stages=k,f=d,theta=theta,K=K,row.names = ""),
                  "Arm Size"=data.frame(m=m,n=n,row.names = ""),
                  "Upper Critical Values"=round(r,4),
                  "Lower Critical Values"=round(c(a[1:(k-1)],r[k]),4),
                  "Type I Error Rate"=round(c(p.r,sum(p.r)),4),
                  "Type II Error Rate"=round(c(p.a,sum(p.a)),4))
      names(out$`Type I Error Rate`) <- c(paste("Stage",1:k),"Overall")
      names(out$`Type II Error Rate`) <- c(paste("Stage",1:k),"Overall")
      return(out)
    }
  }
}
