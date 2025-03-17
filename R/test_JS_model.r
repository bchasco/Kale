ns <- 8 #number of strata
ps <- dpois(1:8, lambda = 3)
ps <- ps / sum(ps)

N <- 1000 #population size

si <- unlist(apply(rmultinom(N, 1, ps),
                    2,
                    function(x){ which.max(x) }))

phi <- 0.8
p <- 0.3
ch_i <- list()
for(i in 1:N){
    tmp_ch <- 0
    jcnt <- 1
    for(j in si[i]:8){
        if(runif()<phi){
            tmp_ch[jcnt] <- 1
        }else{
            tmp_ch[jcnt] <- 1
        }
    }

}
## jan messing around
