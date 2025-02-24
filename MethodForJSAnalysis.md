---
title: "Methods for JS Lewis analysis"
author: "Brandon"
date: "2025-02-20"
output: 
  html_document:
    keep_md: true
---



# Methods

This is an individual based model describing the estimate of carcass abundance in the Lewis River. This model is different than that of Schwarz. The goal is to estimate the total carcass abundance (i.e., $\lambda$ \ref{GenericLikelihood})


To understand the likelihood, it is best to start with a model that has no spatial component.
\begin{equation} \label{eq:Likelihood}
  L = P(B_{total}|\lambda p)\prod_i^{B_{total}} P(t_i|\boldsymbol{\pi})P(D_{i,t}|t_i,\phi_i,p_i)P(T_{i,t}|D_{i,t},p_i)P(R_{i,t\prime}|T_{i,t},\phi_i,p_i) 
\end{equation}

Breaking down the likelihood into it's components, the probability of the observed total number of carcasses is
\begin{equation}\label{P_PopulationSize}
  P(B_{total}|\lambda p)
\end{equation}
where,  $\lambda$ is total population size, and p is the capture probability of a carcass.
The probability of the arrival process is,
\begin{equation}\label{P_arrival}
  P(t_i|\boldsymbol{\pi})
\end{equation}
where, $t_i$ is a latent effect, and $\boldsymbol{\pi}$ vector of arrival probabilities for each time-step.
The probability of the initial detection for any carcass, 
\begin{equation}\label{P_initial detection}
  P(D_{i,t}|t_i,\phi_i,p_i)
\end{equation}
where, $D_{i,t} = 1$, $t_i$ is the latent estimate of arrival time, $phi_i$ is the survival probability, and $p_i$ is the detection probability.
The probability of carcass being tagged is,
\begin{equation}\label{P_tagging}
  P(T_{i,t}|D_{i,t},\psi)
\end{equation}
where, $T_{i,t}=1$ is a carcass was tagged, $\psi$ is the tagging rate.
Finally, there is probability of recapture,

\begin{equation}\label{P_recapture}
  P(R_{i,t\prime}|T_{i,t},\phi_i,p_i)
\end{equation}

where, $R_{i,t\prime} = 1$ for a carcass recaptured on day $t\prime$.

## Create the data

We begin by reading the data and transforming some of the output to make it easier to work in RTMB.


``` r
#Data list
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

``` r
library(tidyr)

d <- read.csv("data/simpleData2.csv") %>%
  mutate(t_wk = lubridate::week(lubridate::mdy(TagDate)),
         r_wk = lubridate::week(lubridate::mdy(RecapDate)),
         t_yr = lubridate::year(lubridate::mdy(TagDate)),
         r_yr = lubridate::year(lubridate::mdy(RecapDate))) %>%
  filter(t_yr == 2024) %>% 
  filter(t_wk > 10) %>%
  mutate(t_k = t_wk - min(t_wk) + 1,
         r_k = r_wk - min(t_wk) + 1) %>%
  mutate(t_l = TagState,
         r_l = RecapState) %>%
  filter(is.na(r_wk) | r_k>0) %>% 
  mutate(tag = ifelse(Tag1=="",FALSE,TRUE)) %>%
  group_by(t_k,r_k,t_l,r_l,tag) %>%
  summarise(n = n())
```

```
## `summarise()` has grouped output by 't_k', 'r_k', 't_l', 'r_l'. You can
## override using the `.groups` argument.
```


## Just model the CJS $P(R_{i,t\prime}|T_{i,t},\phi_i,p_i)$

It is possible to breaks down the analysis in to separate components. We can start by just looking at the Cormack-Jolly-Seber part of the model. For those tagged individuals that are tagged, we can estimate survival and detection probability. For now, we can simply assume all fish are equal and the detection and survival is constant across time and location. The probability of 

I am going to use a matrix-algebra approach to analyzing the data. And start with the simplest data set possible.


``` r
data <- list(t_l = d$t_l, #tagging location
             r_l = d$r_l, #recapture location
             t_k = d$t_k, #tagging week
             r_k = d$r_k,
             tag = d$tag,
             n = d$n) #recapture week, last week if not recapture

# Initial parameter values
parameters <- list(
  phi_par = 0,
  p_par = 0
  )
```

Now lets create the likelihoods function given the data.



``` r
f <- function(parms){

  RTMB::getAll(data,
               parms)
  
  #Negative
  nll <- rep(0,length(d$t_l))
  
  #Survival
  phi <- matrix(0,2,2)
  phi[1,1] <- -exp(phi_par)
  phi[1,2] <- -phi[1,1]
  phi <- Matrix::expm(phi)
  
  #Detection probability
  p <- matrix(0,2,2)
  p[1,1] <- -exp(p_par)
  p[1,2] <- -p[1,1]
  p <- Matrix::expm(p)
  
  #This is the amount of time to first detection
  for(i in 1:length(t_l)){
    
    #Accumulator
    m <- matrix(0,2,2)
    diag(m) <- 1
    
    #Initial state
    delta <- rep(0,2)
    delta[1] <- 1
    if(tag[i] == TRUE ){
      if(is.na(r_k[i])){
        t_prime <- 16
      }else{
        t_prime <- r_k[i]
      }
      for(j in (t_k[i]+1):t_prime){
        if(j<t_prime){ #not detected
          m <- m %*% phi %*% diag(p[,2])
        }else{#last observation
          if(is.na(r_k[i])){#not detected
            m <- m %*% phi %*% diag(p[,2])
          }else{#Detected
            m <- m %*% phi %*% diag(p[,1])
          }
        }
      }
      tmp_phi <- phi * 0
      tmp_phi[,2] <- phi[,2]
      m <- m %*% tmp_phi

    }
    nll[i] <- log(t(delta) %*% m %*% rep(1,2))
  }

  RTMB::REPORT(nll)  
  RTMB::REPORT(phi)  
  RTMB::REPORT(p)  
  return(-sum(nll*n))
  # return(0)
}
```

Next we can optimize.


### Results

Finally, we can look at the results of $\phi$ and $p$ matrixes,


``` r
#Survival/persistence
print(round(rep$phi,2))
```

```
## 2 x 2 Matrix of class "dgeMatrix"
##      [,1] [,2]
## [1,] 0.35 0.65
## [2,] 0.00 1.00
```

``` r
#Survival/persistence
print(round(rep$p,2))
```

```
## 2 x 2 Matrix of class "dgeMatrix"
##      [,1] [,2]
## [1,] 0.71 0.29
## [2,] 0.00 1.00
```

In this simple model, for the tagged individuals, the "persistence/survival" between time steps is 0.35 and the detection probability is 0.71.


## Just model the CJS $P(R_{i,t\prime}|T_{i,t},\phi_i,p_i)$ plus tagging rate $P(T_{i,t}|D_{i,t},\psi)$ 


``` r
# Initial parameter values
parameters <- list(
  phi_par = 0,
  p_par = 0,
  psi_par = 0
  )
```


``` r
f <- function(parms){

  RTMB::getAll(data,
               parms)
  
  #Negative
  nll <- rep(0,length(d$t_l))
  nll_tag <- rep(0,length(d$t_l))   
  
  #Survival
  phi <- matrix(0,2,2)
  phi[1,1] <- -exp(phi_par)
  phi[1,2] <- -phi[1,1]
  phi <- Matrix::expm(phi)
  
  #Detection probability
  p <- matrix(0,2,2)
  p[1,1] <- -exp(p_par)
  p[1,2] <- -p[1,1]
  p <- Matrix::expm(p)

  
  #This is the amount of time to first detection
  for(i in 1:length(t_l)){
  
    #Probability of tagging
    nll_tag[i] <-  RTMB::dbinom(tag[i]*n[i],
                             n[i],
                             RTMB::plogis(psi_par), 
                             log = TRUE) 
    #Accumulator
    m <- matrix(0,2,2)
    diag(m) <- 1
    
    #Initial state
    delta <- rep(0,2)
    delta[1] <- 1
    if(tag[i] == TRUE ){
      if(is.na(r_k[i])){
        t_prime <- 16
      }else{
        t_prime <- r_k[i]
      }
      for(j in (t_k[i]+1):t_prime){
        if(j<t_prime){ #not detected
          m <- m %*% phi %*% diag(p[,2])
        }else{#last observation
          if(is.na(r_k[i])){#not detected
            m <- m %*% phi %*% diag(p[,2])
          }else{#Detected
            m <- m %*% phi %*% diag(p[,1])
          }
        }
      }
      tmp_phi <- phi * 0
      tmp_phi[,2] <- phi[,2]
      m <- m %*% tmp_phi

    }
    nll[i] <- log(t(delta) %*% m %*% rep(1,2))
  }

  psi <- RTMB::plogis(psi_par)
  
  RTMB::REPORT(nll)  
  RTMB::REPORT(nll_tag)
  RTMB::REPORT(phi)  
  RTMB::REPORT(p)  
  RTMB::REPORT(psi)
  return(-sum(nll*n) - sum(nll_tag))
  # return(0)
}
```

Next we can optimize.



### Results

Not surprisingly, the values for the survival ($\phi$) and detection probability ($p$), 0.35, 0, 0.65, 1 and 0.71, 0, 0.29, 1, respectively, remain unchanged, while the estimate of tagging rate ($\psi$) is equal to 0.27. The estimate for $\psi$ is identical to 3100, the number of carcasses tagged, divided by the 11290, the number of carcasses observed. 


# Tables








