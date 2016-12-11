#### httk tests ####

# install.packages("abind") #for the step of plotting polygons - we can include in the package later
library(abind)

# parametrization test (we have added Vmax, km) -----------------------------------

inits2 <- parameterize_pbtk(chem.cas=NULL,chem.name = "Chlorpyrifos", species = "Human", default.to.human = F,
                           tissuelist = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                           force.human.clint.fub = F, clint.pvalue.threshold = 0.05)

#single call to introduce_variability function:
introduce_variability(physiology.data, cv=0) #original df
introduce_variability(physiology.data, cv=.1) #same CV on everything
cv.default <- c("Total Body Water" = .3,
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
introduce_variability(physiology.data, cv=cv.default)


# variability - test/demo -----------------------------------------------
pbtk_mclist <- list("inits"=list(), "pbtk_result"=list(), "halflife"=list())
# this is a slow, but for now sufficient implementation of how we can use this function
# I got about 1 minute / 1,000 simulations
for(i in 1:100) {
    inits <- parameterize_pbtk(chem.cas=NULL,chem.name = "Bisphenol-A", species = "Human", default.to.human = F,
                               tissuelist = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                               force.human.clint.fub = F, clint.pvalue.threshold = 0.05, monte.carlo=TRUE)
    # leaving monte.carlo.cv above empty means inits are generated as defaults
    pbtk_mclist[["inits"]][[i]] <- inits
    sol <- solve_pbtk(parameters = inits,
                        chem.name='Bisphenol-A',
                        plots=FALSE,
                        suppress.messages = TRUE,
                        doses.per.day=3,
                        dose=1/3)
    pbtk_mclist[["pbtk_result"]][[i]] <- sol
    # halflife calculation here? or is there no variability in half-life?
    sol_time <- sol[,"time"]
    pbtk_mclist[["halflife"]][[i]] <- apply(sol, 2, function(x) {
            return(0.693/lm(x ~ sol_time)$coefficients[["sol_time"]])
        })
}
# class(pbtk_mclist) <- "httkmc"

# plotting results (not a separate function for now, since this will need to be upgraded):
res <- apply(abind::abind(pbtk_mclist[["pbtk_result"]], along=3), c(1,2), function(x) {
    c("mean"=mean(x, na.rm=T), "lci"=quantile(x,.025, na.rm=T), "uci"=quantile(x,.975, na.rm=T))
})
timevar <- res["mean",,"time"]
# select <- dimnames(res)[3][[1]][-1]
select <- c("Cplasma", "AUC") #etc.
if(!is.null(select)) res <- res[,,select] #you can select some of the columns only

par(mfrow=c(ceiling(dim(res)[3]),3))
for(cd in 1:dim(res)[3]) {
    plot(res["mean",,cd] ~ timevar, type="l", xlab="time", ylab=dimnames(res)[3][[1]][cd])
    polygon(c(timevar, rev(timevar)), c(res["lci.2.5%",,cd], rev(res["uci.97.5%",,cd])), col="gray", border=NA)
    lines(res["mean",,cd] ~ timevar, type="l", lwd=1.2)
}

# plotting distributions of half-lives + summary stats (some of them are meaningless, ie time ~ time):
halflife <- do.call(rbind, pbtk_mclist[["halflife"]])
apply(halflife, 2, function(x) c("lci"=quantile(x, .025), "mean"=mean(x), "uci"=quantile(x, .975)))
# apply(halflife, 2, function(x) hist(x, breaks=100))
par(mfrow=c(3,5))
counter <- 0
apply(halflife, 2, function(x) {
        counter <<- counter + 1
        plot(density(x), main = colnames(halflife)[counter])
    })



# replace GFR with renal flow ---------------------------------------------
inits <- parameterize_pbtk(chem.cas=NULL,chem.name = "Chlorpyrifos", species = "Human", default.to.human = F,
                            tissuelist = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                            force.human.clint.fub = F, clint.pvalue.threshold = 0.05, use.qrenal = TRUE)




# Adding variability on CLh too. ------------------------------------------

# single call with varibility on CLh only
inits <- parameterize_pbtk(chem.cas=NULL,chem.name = "Bisphenol-A", species = "Human", default.to.human = F,
                           tissuelist = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                           force.human.clint.fub = F, clint.pvalue.threshold = 0.05, clh.cv = 0.1)

# variability on both on physiological data and clh.cv (note manual specification of CV values, even though I use the defaults)
inits <- parameterize_pbtk(chem.cas=NULL,chem.name = "Bisphenol-A", species = "Human", default.to.human = F,
                           tissuelist = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                           force.human.clint.fub = F, clint.pvalue.threshold = 0.05, 
                           monte.carlo = TRUE, 
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
                                            "Average Body Temperature" = 0),
                           clh.cv = 0.1)

