\name{rayleigh}
\alias{rayleigh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates reyleigh number, given a period
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rayleigh(model, par, init, t_0, t_1, DT, n = as.integer(250), ncores = 1, ...)
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
  \item{DT}{
%%     ~~Describe \code{DT} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (model, par, init, t_0, t_1, DT, n = as.integer(250), 
    ncores = 1, ...) 
{
    N <- nrow(par)
    state <- matrix(init, N, length(init), byrow = TRUE)
    state <- propagate_1step_D_multicore(model, NULL, state, 
        par, t_0, t_1, ncores = ncores, ...)$state
    if (t_1 <= t_0) 
        stop("t_1 must be greater than t_0")
    if (DT <= 0) 
        stop("DT must be strictly positive")
    if (as.integer(n) <= 0) 
        stop("as.integer(n) must be at least 1")
    norm <- function(state) sqrt(apply(state * state, 1, sum))
    sum_states = state * 0
    sum_norm = rep(0, N)
    for (i in seq(0, n)) {
        t <- t_1 + i * DT
        t_2 <- t + DT
        state <- propagate_1step_D_multicore(model, NULL, state, 
            par, t, t_2, ncores = ncores, ...)$state
        sum_states <- sum_states + state
        sum_norm <- sum_norm + norm(state)
    }
    return(norm(sum_states)/sum_norm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
