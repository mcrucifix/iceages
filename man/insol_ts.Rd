\name{insol_ts}
\alias{insol_ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes time-series of precession and obliquity
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
insol_ts(times, Astro)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{Astro}{
%%     ~~Describe \code{Astro} here~~
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
function (times, Astro) 
{
    N = length(times)
    pre = rep(0, N)
    forcing = as.numeric(rep(0, 3))
    dforcingdt = as.numeric(rep(0, 3))
    cpre = obl = pre
    for (i in seq(along = times)) {
        INSOL = .Fortran("astro", Astro$nap, Astro$nao, as.numeric(times[i]), 
            Astro$amppre, Astro$omepre, Astro$angpre, Astro$ampobl, 
            Astro$omeobl, Astro$angobl, forcing = forcing, dforcingdt = dforcingdt, 
            DUP = FALSE)
        pre[i] <- INSOL$forcing[1]
        cpre[i] <- INSOL$forcing[2]
        obl[i] <- INSOL$forcing[3]
    }
    list(pre = pre, cpre = cpre, obl = obl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
