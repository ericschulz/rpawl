\name{getFrequencies}
\alias{getFrequencies}
\title{
    Observed Frequencies in each bin.
}
\description{
    This function provides a convenient way to check whether the target frequencies have been reached.
    Since new bins can be created during the algorithm, this function aggregates them in the right way
    so that the user can know if the initial bins (on which the desired frequencies were specified) 
    have been visited enough.
}
\usage{
    getFrequencies(results, binning)
}
\arguments{
    \item{results}{
        Object of class \code{"list"}: either the output of \code{\link{pawl}} or of \code{\link{adaptiveMH}}.
    }
    \item{binning}{
        Object of class \code{\link{binning}}: the binning on which the frequencies have to be computed.
    }
}
\value{
    The function is supposed to be used for the lines that it prints in the console.
    However it also returns a vector of sampling frequencies corresponding to the initial bins.
}
\author{
  Luke Bornn <l.bornn@stat.ubc.ca>,
  Pierre Jacob <pierre.jacob@ensae.fr>
}
%\note{
%%%  ~~further notes~~
%}
%
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pawl}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.



