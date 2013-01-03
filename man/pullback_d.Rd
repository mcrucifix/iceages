\name{pullback_d}
\alias{pullback_d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates determintistic pullback attractor over a series of timesteps
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pullback_d(model, par, Astro, t_back = -500, times = seq(0, 100, 0.2), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{par}{
%%     ~~Describe \code{par} here~~
}
  \item{Astro}{
%%     ~~Describe \code{Astro} here~~
}
  \item{t_back}{
%%     ~~Describe \code{t_back} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
data(models)
times <- seq(0,100,1)
init  <- c(0.2, 0.2, 0.2)
Astro <- read_astro(34,34)
PROP <- propagate_d(models$pp04, init=init, par=c(1., 1.), times=times, Astro=Astro);
plot(times, PROP[,1])
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
