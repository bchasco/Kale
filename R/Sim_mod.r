rm(list=ls())
dist_mat <- matrix(c(0,1,
                     1,0),
                   ncol = 2,
                   byrow = TRUE)

dect <- qlogis(c(0.8,0.5,0))

#Survival matrix
trans_par <- qlogis(1-0.9)

f_gamma <- function(dist_mat,
                    trans_par){
  t_gam <- RTMB::matrix(0,2,2)
  idx <- 1
  for(i in 1:(ncol(dist_mat)-1)){
    for(j in 2:ncol(dist_mat)){
      t_gam[i,j] <- RTMB::plogis(trans_par[idx] * dist_mat[i,j])
      idx <- idx + 1
    }
    t_gam[i,i] <- 1 - t_gam[i,2]
  }
  t_gam[ncol(dist_mat),ncol(dist_mat)] <- 1
  return(t_gam)
}
true_gam <- f_gamma(dist_mat, trans_par)

#Detection matrix
f_omega <- function(par_om){
  t_om <- RTMB::matrix(0,2,3,byrow=TRUE)
  for(k in 1:nrow(t_om)){
    t_om[k,k] <- RTMB::plogis(par_om[k])
    t_om[k,3] <- 1-t_om[k,k]
  }
  om <- (t_om)
  return(t_om)
}
true_om <- f_omega(dect)


n_fish <- 1000
n_weeks <- 8
fish_i <- list()
state_i <- list()
wk_i <- list()
gam <- f_gamma(dist_mat,
               trans_par)
om <- f_omega(dect[1:2])
for(i in 1:n_fish){
  draw <- rmultinom(1,1,rep(1/(ncol(gam)),ncol(gam)))
  init_state <- which(draw==1, arr.ind = TRUE)[1,1]

  state <- init_state
  obs <- init_state

  delta <- rep(0,ncol(gam))
  delta[init_state] <- 1
  for(wk in 2:n_weeks){
    tmp_state <- delta %*% gam
    state[wk] <- which(rmultinom(1,1,as.vector(gam[state[wk-1],]))==1, arr.ind=TRUE)[1,1]
    obs[wk] <- which(rmultinom(1,1,as.vector(om[state[wk],]))==1, arr.ind=TRUE)[1,1]
    if(obs[wk]< 3){
      break;
      obs[wk + 1] <- 3
    }
  }
  fish_i[[i]] <- obs
  state_i[[i]] <- state
}

fwd_algorithm <- function(y, Gamma, Omega, delta){

  T <- length(y)
  if(T < 1) return(0)

  ns <- ncol(Gamma)

  # Track the likelihoods for the states at each observation
  alpha <- matrix(0, nrow=T, ncol=ns)

  # Initialize alpha_1(s)
  obs1 <- y[1]
  for(s in 1:ncol(Gamma)){
    alpha[1,s] <- delta[s] * Omega[s, obs1]
  }

  # 2) Forward recursion for t=2..T
  for(t in 2:T){
    obs_t <- y[t] #This is observed state
    if(obs_t!=3){
      states <- obs_t
    }else{
      states <- 1:ns
    }
    for(s in states){
      tmp_sum <- 0
      for(r in 1:ns){
        #State process
        tmp_sum <- tmp_sum + alpha[t-1, r] * Gamma[r,s]
      }
      #Observation process
      alpha[t,s] <- tmp_sum * Omega[s, obs_t]
    }
  }
  # The total probability of the observed sequence = sum(alpha[T, s])
  return(sum(alpha[T,]))
}

obs_dist <- table(unlist(lapply(fish_i, function(x){paste(x,collapse = "_")})))
obs_i <- names(obs_dist)

data <- list(y = as.vector(obs_dist),
             obs_i = obs_i,
             dist_mat = as.matrix(dist_mat))
parameters <- list(t_vars = rep(-3,1),
                   d_vars = rep(1,2))

# Model of the negative log-likelihood of the data
f <- function(parms){

  RTMB::getAll(data, parms)
  y <- OBS(y)

  f_gamma <- function(dist_mat,
                      trans_par){
    t_gam <- RTMB::matrix(0,2,2)
    idx <- 1
    for(i in 1:(ncol(dist_mat)-1)){
      for(j in 2:ncol(dist_mat)){
        t_gam[i,j] <- exp(trans_par[idx]) * dist_mat[i,j]
        idx <- idx + 1
      }
      t_gam[i,i] <-  -t_gam[i,2]
    }
    t_gam <- Matrix::expm(t_gam)
    return(t_gam)
  }
  Q <- f_gamma(dist_mat, t_vars)

  f_omega <- function(par_om){
    t_om <- RTMB::matrix(0,2,3,byrow=TRUE)
    for(k in 1:nrow(t_om)){
      t_om[k,k] <- RTMB::plogis(par_om[k])
      t_om[k,3] <- 1-t_om[k,k]
    }
    om <- (t_om)
    return(t_om)
  }
  Om <- f_omega(d_vars)

  nll <- 0
  pi_i <- rep(0, length(obs_i))
  for(i in seq_along(obs_i)){
    y_i <- as.integer(strsplit(obs_i[i],"_")[[1]])  # e.g. c(1,3,1,...) in {1,2,3}
    delta <- rep(0,2)
    delta[y_i[1]] <- 1
    pi_i[i] <- fwd_algorithm(y_i, Q, Om, delta)
  }
  REPORT(pi_i)
  REPORT(Q)
  REPORT(Om)
  pi_i <- pi_i/sum(pi_i)
  y %~% dmultinom(prob = pi_i)
}

obj <- MakeADFun(f, parameters)
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he))
rep <- obj$report()
print(round(rep$Q,3))
print(true_gam)
print(round(rep$Om,3))
print(true_om)
plot(round(rep$pi_i/sum(rep$pi_i),2), data$y/sum(data$y))
# chk <- checkConsistency(obj)
