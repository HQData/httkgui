\name{llTrapAUC}
\alias{llTrapAUC}
\title{
    Calculate AUC using trapezoid and logarithmic trapezoid methods.
}
\description{
    %%  ~~ A concise (1-5 lines) description of what the function does. ~~
        Provided time and value vectors, will pick linear or logarithmic term for trapezoid estimate of AUC based on whether values go up or down at particular step.
}
\usage{
    llTrapAUC(time, value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{value}{Vector of values over time}
    \item{time}{Vector of times}
}

\details{
    ...
}

\author{
    Witold Wiecek / LA-SER
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Retrieval}
