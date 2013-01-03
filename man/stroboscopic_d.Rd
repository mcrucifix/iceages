\name{stroboscopic_d}
\alias{stroboscopic_d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates stroboscopic section of a deterministic model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stroboscopic_d(model, par, Astro, t_back = -500, t_I = 0, DT = 1.89760722271748, n = 250, ...)
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
  \item{t_I}{
%%     ~~Describe \code{t_I} here~~
}
  \item{DT}{
%%     ~~Describe \code{DT} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (model, par, Astro, t_back = -500, t_I = 0, DT = 1.89760722271748, 
    n = 250, ...) 
{
    times = t_I + seq(n) * DT
    pullback_d(model, par, Astro, t_back, times = times, ...)$S
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
