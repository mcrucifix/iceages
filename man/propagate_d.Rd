\name{propagate_d}
\alias{propagate_d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Propagate determinstic model forward over a series ot timesteps
}
\description{
Essentially a wrapper to \code{\link{propagate_1step_D}}. 
Propagate determinstic model forward over a series ot timesteps. The numerical schemeis supplied
by the model (see original codes)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
propagate_d(model, times, init, par, Astro, ...)
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
  \item{\dots}{ Other parameters supplied to \code{\link{propagate_1step_D}}
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Matrix, with columns os follows:
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
 \code{\link{propagate_s}} (stochastic version)
}
\examples{
## check convergence of time series
data(models)

times = seq(-80,0,0.2)
deltat = c(0.01, 0.02, 0.1)

# uses 12 components of precession and obliquity
Astro <- read_astro(12,12)

sol = list()
for (i in seq(along=deltat))
 {
  sol[[i]] <- propagate_d (models$vdp_d, models$vdp_d$spar, init=c(0.1, 0.1), 
                          times=times, deltat=deltat[i], Astro=Astro)
 }

plot(times, sol[[1]][,1], type='l')
for (i in seq(along=sol)[-1]) 
 { 
  lines(times, sol[[i]][,1], type='l', col=i)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
