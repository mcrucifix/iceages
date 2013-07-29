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

# modif history
# 22 march : manages input/output of Wiener process increment (dw)
#            in stochastic propagator

propagate_1step_D <-
function (model="t06_f", 
      N, state, par, told, tnew, icalclyap=as.integer(0), ds=0., lyap=0., deltat=0.003, Astro=NULL) 
  {
     Ns <- size_npar_nrow (N,state,par)
     if (is.vector(par)) par <- matrix(par, Ns$N, length(par), byrow=TRUE)
     func <- getNativeSymbolInfo(model$func)$addr
     .Call("c_propagate_d", func, Ns, state, par, told, tnew, icalclyap, ds, lyap, Astro, deltat)
  }
 
propagate_1step_D_multicore <-
function ( model, N, state, par, told, tnew, 
            icalclyap=as.integer(0), ds=NULL, lyap=NULL, deltat=0.003, Astro=NULL, ncores=1)
{           
  if (!require(multicore))
  {       propagate_1step_D (model, NULL, state, par, told, tnew, icalclyap, 
                  ds[Nn, ], lyap[Nn], deltat, Astro)  }
  else
  {
 
    if (is.null(N)  & !is.vector(state)) N = nrow(state);
    if (is.vector(state)) 
      { state <- as.matrix(state, 1, length(state))
        N = 1
      }
    # here we defintely need to have ds and lyap as appropriately sized matrices
 
    if (is.null(ds)) ds <- state*0
    if (is.null(lyap)) lyap <- rep(0, N)
    if (is.vector(par)) par <- matrix(par, N, length(par), byrow=TRUE)
 
    NN <- DistributeN(N, ncores)
 
    p <- list()
 
    for (i in seq(1,NN$ncores))
    {
      Nn <- seq(NN$SNi[i]+1, NN$SNi[i+1])
      
      p[[i]] <- parallel ( 
        propagate_1step_D (model, NULL, state[Nn, ], par[Nn, ], told, tnew, icalclyap, 
                  ds[Nn, ], lyap[Nn], deltat, Astro) )
    }
 
    OUT <- collect(p)
    for (i in seq(1,NN$ncores))
    {
      Nn <- seq(NN$SNi[i]+1, NN$SNi[i+1])
      state [Nn, ] <- OUT[[i]]$state
      ds [Nn, ] <- OUT[[i]]$ds
      lyap [Nn] <- OUT[[i]]$lyap
    }
 
      return (list(state=state, ds = ds, lyap=lyap));
  }

}

propagate_1step_S <-
function (model="vdp_s", 
      N, state, par, told, tnew, ix,  deltat, isum=1, Astro, rmean=NULL,dw=NULL) 
   {
      # if dw not provided, initialise to zero with the same shape as state
      Ns <- size_npar_nrow (N,state,par)
      if (is.vector(par)) 
       { names_par <- names(par) 
         par <- matrix(par, Ns$N, length(par), byrow=TRUE)
         colnames(par) <- names_par
       }
      if (is.null(dw)) dw = as.numeric(state*0.) 
      if (is.null(rmean)) rmean = rep(0.,Ns$ndim)
      func <- getNativeSymbolInfo(model$func)$addr
      if (is.matrix(par))  scaletime = par[,'omega'] else scaletime = par['omega']
      .Call("c_propagate_s", func, Ns, state, par, as.numeric(scaletime), 
                             told, tnew, deltat, ix, isum, Astro, rmean, dw)
   }
