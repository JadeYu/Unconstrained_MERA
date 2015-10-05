source("1sp_growth.R")
R0 <- 100
Dr_seq <- c(1,0.3,0.8,0)
theta <- 1
r <- 0.01

growth_function(r,R0,Dr_seq,theta,nstep=30)

