\name{fastANCOM-package}
\alias{fastANCOM-package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A fast method for analysis of composition of microbiomes
}
\description{
\strong{fastANCOM} We implement the original ANCOM and a fast version fastANCOM under the log-linear model.

    The complete list of functions can be displayed with library(help = fastANCOM).

    More information on fastANCOM can be found :
    \url{https://github.com/ZRChao/fastANCOM}
}
%\usage{
%adaANCOM(x)
%}
%- maybe also 'usage' for other objects documented here.
%\arguments{
%  \item{x}{
%%     ~~Describe \code{x} here~~
%}
%}
\details{
      Package: fastANCOM

      Type: R package

      Version: 0.1.0

      Date: 2021-7-1

      License: GPL-2
}

\author{
Chao Zhou

Maintainer : Chao ZHOU <Supdream8@sjtu.edu.cn>
}


\references{
Mandal, Siddhartha, et al. "Analysis of composition of microbiomes: a novel method for studying microbial composition." Microbial ecology in health and disease 26.1 (2015): 27663.
}



%% ~Make other sections like Warning with \section{Warning }{....} ~



\examples{
data <- matrix(rpois(100*60, 10), 60)
b <- rnorm(10, 4)
data[31:60,1:10] <- t(t(data[31:60,1:10])*2^b)
group <- rep(c(0, 1), each=30)
fit <- fastANCOM(Y=data, x=group)
head(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.




