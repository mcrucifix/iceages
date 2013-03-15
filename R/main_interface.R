# Copyright (c) 2012 Michel Crucifix <michel.crucifix@uclouvain.be>

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject

# the following conditions:

# The above copyright notice and this permission notice shall be
# incluudedin all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND INFRINGEMENT
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR

# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# ------------------------------------------------------------------
# R Code developed for R version 2.15.2 (2012-10-26) -- "Trick or Treat"
# ------------------------------------------------------------------ 


# produce an insolation time series
# -------------------------------------


insol_ts <- function(times,Astro)
{

  N=length(times)

  # times in 10ka units

  pre=rep(0,N)
  forcing=as.numeric(rep(0.,3))
  dforcingdt=as.numeric(rep(0.,3))
  cpre=obl=pre

  for (i in seq(along=times))
  {
    INSOL = .Fortran('astro', Astro$nap, Astro$nao, as.numeric(times[i]), 
                   Astro$amppre, Astro$omepre, Astro$angpre, 
                   Astro$ampobl, Astro$omeobl, Astro$angobl, 
                   forcing=forcing,dforcingdt=dforcingdt, DUP=FALSE)
  #
     pre[i] <- INSOL$forcing[1]
     cpre[i] <- INSOL$forcing[2]
     obl[i] <- INSOL$forcing[3]
   }
   list(pre=pre,cpre=cpre, obl=obl)
}

# ---------------------------------------------------------
# time series from a determinstic model, with astro output
# init : intial conditions ; par : parameters
# -------------------------------------------------------

propagate_d <- function(model,times, init, par, Astro,...)
{  
  ndim <- length(init)
  state <- init
  # copy times and make sure they are numeric 
  times <- as.numeric(times)
  state_out=matrix(0, nrow=length(times), ncol=ndim)
  state_out[1,] = init

  for (i in seq(along=times)[-1])
     {
      state <- propagate_1step_D (model, NULL, state, par, times[i-1], times[i],  Astro=Astro,...)$state
      state_out[i,] <- state
     }
 
   listinsol <- insol_ts(times,Astro)
   return(cbind(state_out, pre=as.matrix(listinsol$pre), obl=as.matrix(listinsol$obl)))
}

# ---------------------------------------------------------
# time series from a stochastic model, with astro output
# -------------------------------------------------------

propagate_s <- function(model,times, init, par, Astro, seed=01425, deltat=.1, isum=1)
{  
  ndim <- length(init)
  state <- init
  # copy times and make sure they are numeric 
  times <- as.numeric(times)
  state_out=matrix(0, nrow=length(times), ncol=ndim)
  state_out[1,] = init
  set.seed(seed)
  for (i in seq(along=times)[-1])
     {
      ix = as.integer(runif(1)*1000)+1
      state <- propagate_1step_S (model, NULL, state, par, times[i-1], 
      times[i], ix, deltat, isum,  Astro=Astro)$state
      state_out[i,] <- state
     }
 
   listinsol <- insol_ts(times,Astro)

   return(cbind(state_out, pre=listinsol$pre, obl=listinsol$obl))
}

# computes pullback attractors (determinstic)

pullback_d <- function (model, par, Astro, t_back=-500, times=seq(0,100,0.2), ...  )
 { 
   t_0 <- times[1]
   K <- basin(model,par, t_back, t_0, Astro,...)

   S <- list()
   for (i in (1:K$nck))  
     { cat(sprintf('attractor %i \n', i))
       S[[i]] <- propagate_d(model, times, K$clusters[i,],  par, Astro,... )
     }
   return(list(K=K,S=S))
 }

# stroboscopic Poincare section


stroboscopic_d <- function (model, par, Astro, t_back=-500, t_I=0, DT=1.89760722271748, n=250, ...) 
 {
   times = t_I + seq(n)*DT
   pullback_d (model, par, Astro, t_back, times=times, ...  )$S
 }

# Rayleigh number

rayleigh <- function (model, par, init, t_0, t_1, DT, n=as.integer(250), ncores=1, ...)

{
  N <- nrow(par)

  state <- matrix(init, N, length(init), byrow = TRUE);

  state <- propagate_1step_D_multicore (model, NULL, state, par, 
            t_0, t_1, ncores=ncores, ...)$state

  if (t_1 <= t_0) stop ('t_1 must be greater than t_0')
  if (DT  <= 0.)  stop ('DT must be strictly positive')
  if (as.integer(n)   <= 0 )  stop ('as.integer(n) must be at least 1')

  norm <- function(state)  sqrt ( apply(state*state, 1, sum) ) 

  sum_states  = state * 0;
  sum_norm = rep(0,N) ; 

  for (i in seq(0,n))
    {
    t <- t_1 + i * DT
    t_2 <- t + DT
    state <- propagate_1step_D_multicore (model, NULL, state, 
              par, t, t_2, ncores=ncores, ...)$state
    sum_states <- sum_states + state
    sum_norm   <- sum_norm   + norm(state)
    }

  return ( norm(sum_states) / sum_norm ) 
 
  }
# Lyapunov plot

lyapunov <- function (model, par, init, t_0=-500, t_1=-400, t_2=0, ncores=1, ...)
{
 N <- nrow(par)
 if (is.vector(init))
  { state <- matrix (init, N, length(init), byrow=TRUE) }
 else if (is.matrix(init))
  { state <- init }
 else stop ('state must be a vector or a matrix')

 ds <- state*0+sqrt(2)
 lyap <- 0. 
 icalclyap <- as.integer(1)

 # 1. find direction 
 OUT <- propagate_1step_D_multicore (model, NULL, state, par, t_0, t_1, 
                                     icalclyap, ds, lyap, ncores=ncores, ...) 
 state <- OUT$state ; ds <- OUT$ds ;
 # 2. compute Lyapunov
 lyap <- propagate_1step_D_multicore (model, NULL, state, par, t_1, t_2, 
                                    icalclyap, ds, lyap, ncores=ncores,...)$lyap
 # scale by deltat
 lyap <- lyap / (t_2 - t_1)

 return(lyap)
} 
