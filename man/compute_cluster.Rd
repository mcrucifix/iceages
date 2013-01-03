\name{compute_cluster}
\alias{compute_cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
automatic computation of the number of pullback attracting trajectories
}
\description{
Wrapper to basin, to count the number of pullback attracting trajectories
for a series of parameters, with multicore support if available
}
\usage{
compute_cluster(model, par, tback, t0, Astro, ncores = 14, deltat = 0.01, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  One of the models out of the packaged model list (\code{data(models)}) or
  user-supplied.
}
  \item{par}{
  Matrix of parameters (possibly with a single row), with row names correspoding
  to parameters to be suppiled to the model.
}
  \item{tback}{
  Initial time, for which a series of initial conditions as defined in \code{model} are set
}
  \item{t0}{
  Time at which states are evaluated and clustered into a number of distinct solutions }
%%     ~~Describe \code{t0} here~~
  \item{Astro}{
  Astronomical forcing list, as supplied by routine \code{read_astro} }
  \item{ncores}{
  Number of cores over which calculations will be distributed. Inactive if \code{package::multicore} 
  is unavailable}
%%     ~~Describe \code{ncores} here~~
  \item{deltat}{
  Optional parameter giving the time step of the numerical scheme for integrating \code{model}.}
%%     ~~Describe \code{deltat} here~~
  \item{\dots}{
  Other models supplied to the \cite{basin} function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of integers of length \code{nrow(par)} with the number of 
distinct solutions as a function of parameters.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
M. Crucifix
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(models)

parlist <- as.list(models$pp04$spar)

parlist$gamma = seq(0.2, 1.2, 0.2)
parlist$omega = seq(0.2, 1.2, 0.2)

pargrid <- as.matrix(expand.grid(parlist))

Astro <- read_astro(1,1)

C <- compute_cluster(models$pp04, pargrid, -100, 0, Astro, ncores=2)

print(C)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
