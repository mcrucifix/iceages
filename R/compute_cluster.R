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

DistributeN <- function(N, ncores)
  { 
    mncores = min(N,ncores)
    n <- floor(N / mncores)
    Ni <- rep(n,mncores)

    ## redistribute extra processes 

    Nextra <- N - n*mncores
    if (Nextra>0) { Ni[1:Nextra] =  Ni[1:Nextra] + 1 }
    SNi <- c(0, cumsum(Ni))
    list(SNi = SNi, ncores=mncores)
    }


compute_cluster <- function (model, par, tback, t0, Astro, ncores=14,deltat=0.01, ...)

## Basin Function : function calculating the number of
## clusters (and, optionnaly, the basins of attraction)
## tback and  t0     : initial and final times
## ncores : number of cores for distributing the task, 
##          in multicore environment
## ... : other parameters passed to the propagation routine

## deltat may also be supplied as a vector (as par), which is convenient
## when one of the parameters represents a time stretch factor

{

## defines the parallel process : will compute, for any combination of par, 
## the number of clusters = number of distinct solutions. 
ClusterMatrix <- function(model, par, tback, t0, Astro, deltat, ...) 
 {
   N = nrow(par)
   K = rep(0,nrow(par))
   for (i in (1:N)) K[i]=basin(model, as.numeric(par[i,]),tback, t0, Astro, deltat=deltat[i],...)$nc
   K
 }


## launches the different parallel processes

N <- nrow(par)
if (length(deltat)==1) deltat=rep(deltat, N)

## do we use more than one core (in which case requests multicore) ?
## and if multicore is indeed available

if (ncores > 1 & require(multicore)) 

## devides N into a number of parallel processes
{ 
 
  P <- DistributeN(N,ncores)

  p <- list()

  for (i in seq(1,ncores))
   {
     Nn     <- seq(P$SNi[i]+1, P$SNi[i+1])
     p[[i]] <- parallel( ClusterMatrix(model, par[Nn,], tback, t0, Astro, deltat[Nn], ...) )
   }

   # and collect the information 
   OUT <- collect(p)


   K <- OUT[[1]]

   for (i in seq(2,P$ncores)) K <- c(K, OUT[[i]])
   K
}
else
{
   if (ncores > 1) message ('multicore not available so will run on 1 core')
   K <-   ClusterMatrix(model, par, tback, t0, Astro, deltat, ...) 
}
}


