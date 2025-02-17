library(sads)
library(Seurat)
library(emdbook)
library(parallel)

# Load in training data and set initial values
load('training_data.rda')
a.init=2240958
mu.init=c(-12.73153,-8.301349)
sigma.on.init=1.936757
sigma.off.init=2.025042
rho.init=0.5

# Train new data with old data
# data_list = list of counts matrices (genes by cells)
# labels_list = list of cell type/cluster labels in the same order 
# file = path to save params to
# mc.cores = number of cores to use for parallelization
new_params <- function(data_list,labels_list,file,mc.cores) {
  common_genes <- intersect(rownames(panglao),
                            Reduce(intersect,lapply(data_list,function(data)
                              rownames(data))))
  processed_data <- lapply(1:length(data_list),function(x) {
    data <- data_list[[x]]
    labels <- as.character(labels_list[[x]])
    grouped_indices <- split(1:ncol(data),labels)
    data <- sapply(unique(labels), function(id) {
      cols <- grouped_indices[[id]]
      rs <- rowSums(data[, cols, drop = FALSE])
      N <- sum(rs)
      rs[rs==0] <- 1
      rs[common_genes]/N
    })
  })
  all_data <- cbind(panglao[common_genes,],Reduce(cbind,processed_data))
  params <- train_all_genes(all_data,mc.cores)
  save(params,file=file)
}

# Train all genes
train_all_genes <- function(y,mc.cores) {
  fits <- mclapply(1:nrow(y),function(j) train_one_gene(y[j,,drop=F]),
                   mc.cores=mc.cores)
  a.all.g <- sapply(fits,function(j) j[[1]])
  sigma.all.g <- t(sapply(fits,function(j) j[[3]]))
  g.on.all.g <- sapply(fits,function(j) j[[4]])
  g.off.all.g <- sapply(fits,function(j) j[[5]])
  pi.all.g <- t(sapply(fits,function(j) j[[6]]))
  names(a.all.g) <- rownames(sigma.all.g) <- names(g.on.all.g) <-
    names(g.off.all.g) <- rownames(pi.all.g) <- rownames(y)
  return(list(pi.all.g,g.on.all.g,g.off.all.g,a.all.g,sigma.all.g))
}

# Train one gene
train_one_gene <- function(y,a=a.init,mu=mu.init,sigma.on=sigma.on.init,
                  sigma.off=sigma.off.init,rho=rho.init) {
  K <- ncol(y)
  J <- nrow(y)
  
  pi <- c(1/3,1/3)
  g.on <- rep(0,J)
  g.off <- rep(0,J)
  gamma <- list(array(1/3,dim=c(J,K)),array(1/3,dim=c(J,K)))
  sigma <- c(1,1)
  a.init <- a
  y <- matrix(as.numeric(y),nrow=J,ncol=K)
  loglik <- loglikelihood(y,pi,a,mu,sigma,gamma,g.on,g.off,sigma.on,sigma.off,rho)
  diff <- 1000
  
  while(diff>0.1) {
    gamma <- E_step(y,pi,a,mu,sigma,gamma,g.on,g.off)
    a <- update_a(y,gamma,a.init)
    sigma <- update_sigma(y,mu,gamma,sigma,g.on,g.off)
    g.on <- update_gon(y,mu,sigma,gamma,g.off,sigma.on,sigma.off,rho)
    g.off <- update_goff(y,mu,sigma,gamma,g.on,sigma.off,sigma.on,rho)
    pi <- update_pi(gamma)
    loglik.new <- loglikelihood(y,pi,a,mu,sigma,gamma,g.on,g.off,sigma.on,sigma.off,rho)
    diff <- abs(loglik.new-loglik)
    loglik <- loglik.new
  }
  
  return(list(a,mu,sigma,g.on,g.off,pi))
}

# Helper functions for training

E_step <- function(y,pi,a,mu,sigma,gamma,g.on,g.off) {
  J <- nrow(y)
  K <- ncol(y)
  gamma.new <- list(array(0,dim=c(J,K)),array(0,dim=c(J,K)))
  gamma.new[[1]] <- pi[1]*dexp(as.matrix(y),a)/
    (pi[1]*dexp(as.matrix(y),a)+pi[2]*dlnorm(as.matrix(y),mu[1]+g.off,sigma[1])+
       (1-pi[1]-pi[2])*dlnorm(as.matrix(y),mu[2]+g.on,sigma[2]))
  gamma.new[[2]] <- pi[2]*dlnorm(as.matrix(y),mu[1]+g.off,sigma[1])/
    (pi[1]*dexp(as.matrix(y),a)+pi[2]*dlnorm(as.matrix(y),mu[1]+g.off,sigma[1])+
       (1-pi[1]-pi[2])*dlnorm(as.matrix(y),mu[2]+g.on,sigma[2])) 
  return(gamma.new)
}

update_gon <- function(y,mu,sigma,gamma,g.off,sigma.on,sigma.off,rho) {
  J <- nrow(y)
  K <- ncol(y)
  g <- max(((rho*g.off/((1-rho^2)*sigma.off*sigma.on))+sum((1-gamma[[1]]-gamma[[2]])*(log(y)-mu[2])/sigma[2]^2))/
             ((1/((1-rho^2)*sigma.on^2))+sum((1-gamma[[1]]-gamma[[2]])/sigma[2]^2)),g.off+mu[1]-mu[2])
  return(g)
}

update_goff <- function(y,mu,sigma,gamma,g.on,sigma.off,sigma.on,rho) {
  J <- nrow(y)
  K <- ncol(y)
  g <- min(((rho*g.on/((1-rho^2)*sigma.off*sigma.on))+sum(gamma[[2]]*(log(y)-mu[1])/sigma[1]^2))/
             ((1/((1-rho^2)*sigma.off^2))+sum(gamma[[2]]/sigma[1]^2)),g.on+mu[2]-mu[1])
  return(g)
}

update_sigma <- function(y,mu,gamma,sigma,g.on,g.off) {
  K <- ncol(y)
  sigma.new <- rep(0,2)
  sigma.new[1] <- sqrt(max(0.5^2,(sum(gamma[[2]]*(log(y)-mu[1]-g.off)^2)/(sum(gamma[[2]])))))
  sigma.new[2] <- sqrt(max(0.5^2,sum((1-gamma[[1]]-gamma[[2]])*(log(y)-mu[2]-g.on)^2)/sum(1-gamma[[1]]-gamma[[2]])))
  return(sigma.new)
}

update_pi <- function(gamma) {
  K <- ncol(gamma[[1]])
  J <- nrow(gamma[[1]])
  pi <- rep(0,2)
  pi[1] <- sum(gamma[[1]])/K
  pi[2] <- sum(gamma[[2]])/K
  return(pi)
}

update_a <- function(y,gamma,a.init) {
  a <- sum(gamma[[1]])/(sum(gamma[[1]]*y))
  if (is.na(a)) { a <- a.init }
  return(a)
}

loglikelihood <- function(y,pi,a,mu,sigma,gamma,g.on,g.off,sigma.on,sigma.off,rho) {
  K <- ncol(y)
  J <- nrow(y)
  loglik <- sum(log(pi[1]*dexp(as.matrix(y),a)+pi[2]*dlnorm(as.matrix(y),mu[1]+g.off,sigma[1])+
                      (1-pi[1]-pi[2])*dlnorm(as.matrix(y),mu[2]+g.on,sigma[2])))+
    sum(dmvnorm(cbind(g.off,g.on),mu=c(0,0),
                Sigma=array(c(sigma.off^2,rho*sigma.off*sigma.on,
                              rho*sigma.off*sigma.on,sigma.on^2),dim=c(2,2)),log=TRUE))
  return(loglik)
}