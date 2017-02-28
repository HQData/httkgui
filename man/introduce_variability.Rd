\name{introduce_variability}
\alias{introduce_variability}
\title{
    Introduce variability in a data frame containing mean estimates.
}
\description{
    %%  ~~ A concise (1-5 lines) description of what the function does. ~~
        Allows user to supply data frame with mean parameter values and a SD or CV information for these parameters. 
        Function processes each column (corresponding e.g. to different species) and draws from a normal or log-normal distribution with same mean and appropriate SD/CV>
}
\usage{
    introduce_variability(input.mean, 
                          input.sd=NULL, 
                          cv=NULL, 
                          log = FALSE,
                          input.name.col = "Parameter") 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{input.mean}{Input data frame with the means.} 
    \item{input.sd}{Input data frame with SDs for each of the elements (same dimension as numerical part of input.mean).}
    \item{cv}{(Named) vector of coefficients of variation. Names should correspond to parameter names specified by input.name.col.}
    \item{log}{Whether draws should be made from log-normal distribution. If not, normal will be used.}
    \item{input.name.col}{Column in input.mean storing names of parameters.}
    
}

\details{
    ...
}

\author{
    Witold Wiecek / LA-SER
}

\examples{
    monte.carlo.cv = c("Total Body Water" = .3,
                       "Plasma Volume" = .3,
                       "Cardiac Output" = .3,
                       "Average BW" = .16, 
                       "Total Plasma Protein" = .14,
                       "Plasma albumin"= .1,
                       "Plasma a-1-AGP"= .3,
                       "Hematocrit"= .3,
                       "Urine"= .3,
                       "Bile"= .3,
                       "GFR"=.3,
                       "Average Body Temperature" = 0)
    ll <- list()
    for(i in 1:10000) 
        ll[[i]] <- introduce_variability(input.mean = physiology.data, 
                                         cv = monte.carlo.cv, 
                                         log = TRUE)[["Human"]]
    
    ll2 <- list()
    for(i in 1:10000) 
        ll2[[i]] <- introduce_variability(input.mean = physiology.data, 
                                          cv = monte.carlo.cv, 
                                          log = FALSE)[["Human"]]
    physiology.data[["Human"]]
    
    apply(do.call(rbind, ll), 1, mean) 
    apply(do.call(rbind, ll), 2, mean) 
    apply(do.call(rbind, ll), 1, sd) 
    apply(do.call(rbind, ll), 2, sd)
    
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Retrieval}
