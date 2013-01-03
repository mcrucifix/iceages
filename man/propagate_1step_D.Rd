\name{propagate_1step_D}
\alias{propagate_1step_D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Propagates deterministic model foreworkd from told to tnew
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
propagate_1step_D(model = "t06_f", N, state, par, told, tnew, icalclyap = as.integer(0), ds = 0, lyap = 0, deltat = 0.003, Astro = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{state}{
%%     ~~Describe \code{state} here~~
}
  \item{par}{
%%     ~~Describe \code{par} here~~
}
  \item{told}{
%%     ~~Describe \code{told} here~~
}
  \item{tnew}{
%%     ~~Describe \code{tnew} here~~
}
  \item{icalclyap}{
%%     ~~Describe \code{icalclyap} here~~
}
  \item{ds}{
%%     ~~Describe \code{ds} here~~
}
  \item{lyap}{
%%     ~~Describe \code{lyap} here~~
}
  \item{deltat}{
%%     ~~Describe \code{deltat} here~~
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
function (model = "t06_f", N, state, par, told, tnew, icalclyap = as.integer(0), 
    ds = 0, lyap = 0, deltat = 0.003, Astro = NULL) 
{
    Ns <- size_npar_nrow(N, state, par)
    if (is.vector(par)) 
        par <- matrix(par, Ns$N, length(par), byrow = TRUE)
    func <- getNativeSymbolInfo(model$func)$addr
    .Call("c_propagate_d", func, Ns, state, par, told, tnew, 
        icalclyap, ds, lyap, Astro, deltat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
