\name{models}
\alias{models}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List of deterministic and stochastic ice age models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data(models)
}

\details{
List of models, as follows:
\itemize{
\item{\code{i80_d} : }{Imbrie and Imbrie (1980) seminal 1-D model.}
\item{\code{vdp_d} : }{Van der pol oscillator with constant bias.}
\item{\code{cr12_d} : }{Generalisation of the van der Pol, with two stable states on the slow manifold}
\item{\code{vdp_s} : }{Stochastic version of the van der Pol oscillator}
\item{\code{cr12_s} : }{Stochastic version of the generalised van der Pol}
\item{\code{sm90_d} : }{Saltzman and Maasch (1990) 3-D model of ice ages}
\item{\code{sm91_d} : }{Saltsman and Maasch (1991) 3-D model of ice ages}
\item{\code{t06_d} : }{Tziperman et al. (2006) hybrid (1.5-D) model of ice ages}
\item{\code{pp04_d} : }{Paillard and Parrenin (2004) 3-D model of ice ages}
\item{\code{pp12_d} : }{Parrenin and Palliard (2012) 1.5 D model of ice ages}
\item{\code{i11_d} : }{Imbrie et al. (2011) 2-D model of ice ages}
}
}
\value{
A list of models, where each model is a list with
\itemize{
\item {\code{func}} {Reference to \code{fortran} or {C} function. See \code{src} directory for examples of deterministic (ending as \code{_d_f.f90}) and stochastic (ending as \code{_s_f.f90}) models}
\item{\code{name  }}  {String giving model name}
\item {\code{spar  }}  {Named vector of standard model parameters. Should contain \code{omega} as the time-scale factor}
\item {\code{initgrid}} {List of \code{n} items, where \code{n} is the dimension of the state of the system. Each item list is a vector of distinct values that will be used to construct a grid of initial conditions when estimating pullback sections and attractors (see \code{\link{basin}}})
}
}
\references{
\enumerate{
\item{J. Imbrie and J. Z. Imbrie, Modelling the climatic response to orbital variations,  Science, 207, 943-953  1980}
\item{B. Saltzman and K. A. Maasch, A first-order global model of late Cenozoic climate,  Transactions of the Royal Society of Edinburgh Earth Sciences, 81, 315-325  1990}
\item{B. Saltzman and K. A. Maasch, A first-order global model of late Cenozoic climate. II Further analysis based on a simplification of the CO_2 dynamics,  Climate Dynamics, 5, 201-210  1991}
\item{E. Tziperman et al., Consequences of pacing the Pleistocene 100 kyr ice ages by nonlinear phase locking to Milankovitch forcing,  Paleoceanography, 21, PA4206  2006}
\item{D. Paillard and F. Parrenin, The Antarctic ice sheet and the triggering of deglaciations,  Earth Planet. Sc. Lett., 227, 263-271  2004}
\item{John Z. Imbrie, Annabel Imbrie-Moore, and Lorraine E. Lisiecki, A phase-space model for Pleistocene ice volume,  Earth and Planetary Science Letters, 307, 94--102  2011}
\item{M. Crucifix, Why glacial-interglacial cycles could be unpredictable, in prep.}
\item{F. Parrenin and D. Paillard, Terminations VI and VIII (&sim; 530 and &sim; 720 kyr BP) tell us the importance of obliquity and precession in the triggering of deglaciations,  Climate of the Past Discussions, 8, 3143--3157  2012}
}}
\author{
M. Crucifix
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pullback_d}}, \code{\link{propagate_s}}
}
\examples{
data(models)
print(models$vdp_d)
# see pullback_d for an example of use of determinsitic models, 
# and propagate_s for a use of stochastic models.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
