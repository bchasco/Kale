# f_gamma <- function(gd,trans){
#   gam <- matrix(0,gd,gd)
#
#   idx <- 1
#   for(i in 2:gd){ #rows
#     for(j in 1:(i-1)){ #cols
#       if(i!=j){
#         gam[i,j] <- -exp(trans[idx])
#         idx <- idx + 1
#       }
#     }
#     gam[i,i] <- -sum(gam[i,])
#   }
#   gam_mat <- Matrix::expm(gam)
#   return(gam_mat)
# }

# f_omega <- function(p){
#   od <- length(p)
#   om <- matrix(0,od,od+1)
#
#   for(i in 1:od){
#     om[i,i] <- RTMB::plogis(p[i])
#     om[i,od+1] <- 1 - RTMB::plogis(p[i])
#   }
#   return(om)
# }

