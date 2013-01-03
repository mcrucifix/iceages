\name{basin}
\alias{basin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
computation of pullback basins of attraction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
basin(model, par, tback, t0, Astro, tol = 0.1, initgrid = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  One of the deterministic models out of the packaged model list (\code{data(models)}) or
  user-supplied.
}
  \item{par}{
  Vector of models parameters
}
  \item{tback}{
  Initial time, for which a series of initial conditions as defined in \code{model} are set
}
  \item{t0}{
  Time at which states are evaluated and clustered into a number of distinct solutions 
}
  \item{Astro}{
  Astronomical forcing list, as supplied by routine \code{\link{read_astro}} }
  \item{tol}{
  Euclidian distance between two solutions, above which they are considered as distinct}
  \item{initgrid}{
  list of vectors defining the initial conditions. Defaults to \code{model$initgrid} in not supplied.
%%     ~~Describe \code{initgrid} here~~
}
  \item{\dots}{
  Other parameters supplied to \cite{propagate_1step_D}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with
\itemize{
\item {\code{state_out}} {Array of final states}
\item {\code{kernel  }}  { Array of states being centers of distinct final solutions (kernels)
                    representing at least 1 \% of the initial conditions}
\item {\code{belong }} {list of :
 \itemize{
 \item {\code{initgrid}} { \code{initgrid} supplied or model-default }
 \item {\code{z}} { Array mappling initial conditions to the kernels }
 }
\item {\code{nc}} {number of distinct kernels representnag at least 1 \% of the initial conditions}
\item {\code{nck}} {number of distinct kernels}
}
}
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
Astro <- read_astro(34,34)
BASIN <- basin(models$pp04, models$pp04$spar, -200, 0, Astro)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
