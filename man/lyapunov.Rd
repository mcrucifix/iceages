\name{lyapunov}
\alias{lyapunov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate greatest Lyapunov exponent over time span
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lyapunov(model, par, init, t_0 = -500, t_1 = -400, t_2 = 0, ncores = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{par}{
%%     ~~Describe \code{par} here~~
}
  \item{init}{
%%     ~~Describe \code{init} here~~
}
  \item{t_0}{
%%     ~~Describe \code{t_0} here~~
}
  \item{t_1}{
%%     ~~Describe \code{t_1} here~~
}
  \item{t_2}{
%%     ~~Describe \code{t_2} here~~
}
  \item{ncores}{
%%     ~~Describe \code{ncores} here~~
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
function (model, par, init, t_0 = -500, t_1 = -400, t_2 = 0, 
    ncores = 1, ...) 
{
    N <- nrow(par)
    if (is.vector(init)) {
        state <- matrix(init, N, length(init), byrow = TRUE)
    }
    else if (is.matrix(init)) {
        state <- init
    }
    else stop("state must be a vector or a matrix")
    ds <- state * 0 + sqrt(2)
    lyap <- 0
    icalclyap <- as.integer(1)
    OUT <- propagate_1step_D_multicore(model, NULL, state, par, 
        t_0, t_1, icalclyap, ds, lyap, ncores = ncores, ...)
    state <- OUT$state
    ds <- OUT$ds
    lyap <- propagate_1step_D_multicore(model, NULL, state, par, 
        t_1, t_2, icalclyap, ds, lyap, ncores = ncores, ...)$lyap
    lyap <- lyap/(t_2 - t_1)
    return(lyap)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
