library(dr)

step.cop<-function(x,y,H,alpha.in,alpha.out,my.range,k){
  x=as.matrix(x)		
  p=NCOL(x)
  n=nrow(x)
  if(k==1){
    lambdar=0
    aa=NULL
    for(j in 1:p){
      slice.1<-sapply(split(scale(x[,j]),as.factor(dr.slices(y,nslices=H)[[1]])),mean,simplify=TRUE)
      aa[j]<-var(slice.1)
    }
    lambdaf<-max(aa)
    id<-which.max(aa)
    cop<-n*(lambdaf-lambdar)/(1-lambdaf)
    if(cop>=qchisq(alpha.in,1)){
      my.current.sel=id
    }else{
      stop("There is no significant predictor!")
    }
  }else{
    my.current.sel=1:k
  }
  my.step=1			
  my.forward="conti"
  chi.in=qchisq(alpha.in,k)
  chi.out=qchisq(alpha.out,k)
  while(my.forward=="conti"&my.step<=my.range){		
    set.all<-1:p
    set.redundant<-setdiff(set.all,my.current.sel)
    pp=length(set.redundant)        
    if(length(my.current.sel)!=1){	   
      lambdar=dr(y~x[,my.current.sel],nslices=H)$evalues[1:k]
    }else{
      lambdar=lambdaf	
    }
    cop=lambdaf=NULL	  	
    for(j in 1:pp){
      ind1<-c(my.current.sel,set.redundant[j])
      xnew=x[,ind1]
      temp<-dr(y~xnew,nslices=H)
      lambdaf<-temp$evalues[1:k]
      cop[j]<-sum(n*(lambdaf-lambdar)/(1-lambdaf)) 	     	
    }
    cop.stata=max(cop[!is.na(cop)])	
    sel<-which(cop==cop.stata)[1]	
    if(cop.stata>=chi.in){
      my.forward="conti"
      my.current.sel<-c(my.current.sel,set.redundant[sel])
      my.backward="conti"
      while(my.backward=="conti"&length(my.current.sel)>2){
        pp=length(my.current.sel)
        cop=NULL
        for(l in 1:pp){
          ind1<-my.current.sel[-l]
          xfull=x[,my.current.sel]
          xreduce=x[,ind1]	
          temp1<-dr(y~scale(xfull),nslices=H)
          temp2<-dr(y~scale(xreduce),nslices=H)
          cop[l]<-sum(n*(temp1$evalues[1:k]-temp2$evalues[1:k])/(1-temp1$evalues[1:k]))
        }
        cop.statd=min(cop[!is.na(cop)])	
        sel<-which(cop==cop.statd)[1]
        if(cop.statd<=chi.out){
          my.backward="conti"
          my.current.sel<-my.current.sel[-sel]
        }else{
          my.backward="stop"
        }
      }
    }else{
      my.current.sel<-my.current.sel
      my.forward="stop"
    }
    my.step=length(my.current.sel)
  }	
  return(my.current.sel=my.current.sel)
}




scalar <- function(mat, m){
  n = nrow(mat)
  scalar_responses <- matrix(NA, n, m)
  # Generate random direction vectors and project Y along these directions
  for (i in 1:m) {
    # Generate a random direction vector with unit length
    p = ncol(mat)
    set.seed(sample(100,1)+p)
    direction <- rnorm(p, mean=0, sd=1)
    direction <- direction / sqrt(sum(direction^2))
    # Project Y along the direction vector
    scalar_responses[, i] <- mat %*% direction
  }
  return(scalar_responses)
}

step.multicop.x <- function(i_y, X, scalar.Y, alpha.in.list, alpha.out.list, k0=10){
  tryCatch({
    x = X
    y = scalar.Y[,i_y]
    GIC_score = NULL; kk = c()
    
    for (i in 1:length(alpha.in.list)){
      my.d = NULL
      for (j in 1:min(k0-2,max(nrow(x),ncol(x)))){
        my.cop.sel = step.cop(x,y,H = 5,alpha.in = alpha.in.list[i],alpha.out = alpha.out.list[i],
                              my.range=my.range,k=j+1)
        my.d[j]=GIC(x,y,my.sel = my.cop.sel, KK = j)
      }
      K = which.min(my.d)
      GIC_score[i] = my.d[K]
      kk = c(kk, K)
    }
    # if (is.na(GIC_score)){print("no K is selected")}
    alpha = which.min(GIC_score)
    alpha.in = alpha.in.list[alpha]; alpha.out = alpha.out.list[alpha]; K = kk[alpha]
    x_sel = step.cop(x,y,H = 5,alpha.in = alpha.in,alpha.out = alpha.out,my.range=100,k=K+1)
    return(x_sel)
  }, error = function(e) {
  })
}

select.idx <- function(result.X){
  x_sel_vote  <- Filter(function(x) !is.null(x), result.X)
  num_elements <- sapply(x_sel_vote, length)
  if (length(num_elements)==0){
    X_sel_final = c()
  }else{
    K_X <- as.integer(names(sort(table(num_elements), decreasing = TRUE)[1]))
    freq_table <- table(unlist(x_sel_vote))
    sorted_freq_table <- sort(freq_table, decreasing = TRUE)
    X_sel_final <- names(sorted_freq_table[1:K_X])
  }
  return(X_sel_final)
}

step.multicop.y <- function(i_x, Y, scalar.X, alpha.in.list, alpha.out.list, k0=8){
  tryCatch({
    x = Y
    y = scalar.X[,i_x]
    GIC_score = NULL; kk = c()
    for (i in 1:length(alpha.in.list)){
      # Select K, the number of principal profile correlation directions
      my.d=NULL
      for(j in 1:min(k0-2,max(nrow(x),ncol(x)))){
        my.cop.sel = step.cop(x,y,H = 5,alpha.in = alpha.in.list[i], alpha.out = alpha.out.list[i],
                              my.range=100,k=j+1)
        my.d[j]=GIC(x,y,my.sel = my.cop.sel, KK = j)
      }
      K = which.min(my.d)
      GIC_score[i] = my.d[K]
      kk = c(kk, K)
    }
    alpha = which.min(GIC_score)
    alpha.in = alpha.in.list[alpha]; alpha.out = alpha.out.list[alpha]; K = kk[alpha]
    x_sel = step.cop(x,y,H = 5,alpha.in = alpha.in, alpha.out = alpha.out,my.range=100,k=K+1)
    return(x_sel)
  }, error = function(e) {
    # cat("skip ",i_x, "\n")
  })
}


GIC<-function(x,y,my.sel,KK){
  x1=x[,my.sel]
  p=ncol(x1)
  n=nrow(x1)
  phi=dr(y~x[,my.sel])$M
  omega=phi+diag(1,p)
  tao=length(eigen(omega)$values>1)
  ss=min(tao,KK)
  theta = eigen(omega)$values[(1+ss):p]
  logL = 50*n/2*sum(log(theta)+1-theta)
  Gk = -(logL-log(n)*KK*(2*p-KK+1)/2)
  return(Gk)
}






