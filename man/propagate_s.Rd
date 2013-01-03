\name{propagate_s}
\alias{propagate_s}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Propagate stochastic model forward over a series ot timesteps
}
\description{
Essentially a wrapper to \code{\link{propagate_1step_S}}. 
Propagate determinstic model forward over a series ot timesteps. The numerical scheme is supplied
by the model (see original codes). A routine \cite{Npropagate_s} should be provided in a 
forthcoming package version, to integrate simultaneously a bunch of initial conditions}
\usage{
propagate_s(model, times, init, par, Astro, seed=01425, deltat=.1,isum= 1)
}
\arguments{
  \item{model}{
  One of the deterministic models out of the packaged model list (\code{data(models)}) or
  user-supplied.
}
%%     ~~Describe \code{model} here~~
  \item{times}{
  Vector of times. \code{times[1]} at which the model trajectory is being evalueted. Is used as the initial condition
}
  \item{init}{ State vector of initial conditions}
%%     ~~Describe \code{init} here~~
  \item{par}{ Parameter vector}
  \item{Astro}{ Astronomical forcing list, as supplied by routine \code{\link{read_astro}} }
  \item{seed}{ Random number seed (is converted to an integer)}
  \item{deltat}{ Time step used supplied to stochsatic propagator}
  \item{isum}{ Every random number is obtained as the sum of \code{isum} random numbers, then divided
               by \code{sqrt{isum}}. See application  in \code{examples} below}.
%%     ~~Describe \code{\dots} here~~
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Matrix, with columns as follows:
\itemize{
\item{\code{Column 1 .. n:}}{States at specified times (\code{n} is the system dimension)}
\item{\code{Column n+1:}}{precession according to the astronomical solution}
\item{\code{Column n+2:}}{obliquity according to the astronomical solution}}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
M. Crucifix (2012)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{propagate_d}} (deterministic version)
}
\examples{
## check convergence of time series

require(iceages)
data(models)
Astro <- read_astro(34,34)
times=seq(-80,0,0.1)


# test brownian bridge
bb <- function(spar, title='Brownian bridge: vdp mode')
{
  bs = c(1,2,5,10,20)
  for (j in seq(along=bs))
  {
    print(bs[j])
    print(0.1/bs[j])
    STOCH <- propagate_s(models$cr12_s, init=c(0., 0.), 
                       par=spar, times=times, Astro=Astro, 
                       isum=bs[j], deltat=0.01*bs[j], seed=132);
    if (j==1)
    {
       plot(times, STOCH[,1], type='l', main =title)
    } else {   
       lines(times, STOCH[,1], type='l', col=j)
   } 
  }
  legend('topright', legend=bs, lty=1, col=seq(along=bs))
}

# van der pol mode

spar1 = c(alpha = 30., beta0 =  0.7, beta1 = 0,  
         beta2 = 0., delta = 1., 
         gammapre = 0.7, gammaobl = 0.6,
         omega = 3.695, sigmax=0.3, sigmay=0.3)

bb(spar1, title='van der pol mode')

# 1-d mode
spar2 = c(alpha = 30., beta0 =  0.0, beta1 = 0.70,  
         beta2 = 1.4, delta = 0., 
         gammapre = 3.696, gammaobl = 3.20,
         omega = 3.695, sigmax=0.8, sigmay=0.8)


bb(spar2, title='1-d mode')

# intermediate

spar = 0.5*(spar1 + spar2)

bb(spar=0.5*(spar1+spar2), title='intermediate')


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
