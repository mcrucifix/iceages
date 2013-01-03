\name{propagate_1step_S}
\alias{propagate_1step_S}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Propagets stochastic model forward from told to tnew
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
propagate_1step_S(model = "vdp_s", N, state, par, told, tnew, ix, deltat, isum = 1, Astro)
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
  \item{ix}{
%%     ~~Describe \code{ix} here~~
}
  \item{deltat}{
%%     ~~Describe \code{deltat} here~~
}
  \item{isum}{
%%     ~~Describe \code{isum} here~~
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
function (model = "vdp_s", N, state, par, told, tnew, ix, deltat, 
    isum = 1, Astro) 
{
    Ns <- size_npar_nrow(N, state, par)
    func <- getNativeSymbolInfo(model$func)$addr
    .Call("c_propagate_s", func, Ns, state, par, told, tnew, 
        deltat, ix, isum, Astro)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
