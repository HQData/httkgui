#### httk tests ####

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


# variability on GFR - test -----------------------------------------------
pbtk_mclist <- list("inits"=list(), "pbtk_result"=list())
# this is a slow, but for now sufficient implementation of how we can use this function
# I got about 1 minute / 1,000 simulations
for(i in 1:1000) {
    inits <- parameterize_pbtk(chem.cas=NULL,chem.name = "Bisphenol-A", species = "Human", default.to.human = F,
                               tissuelist = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                               force.human.clint.fub = F, clint.pvalue.threshold = 0.05, monte.carlo=TRUE)
    # leaving monte.carlo.cv above empty means inits are generated as defaults
    pbtk_mclist[["inits"]][[i]] <- inits
    pbtk_mclist[["pbtk_result"]][[i]] <- solve_pbtk(parameters = inits,
                                                    chem.name='Bisphenol-A',
                                                    plots=FALSE,
                                                    suppress.messages = TRUE,
                                                    doses.per.day=3,
                                                    dose=1/3)
}




# replace GFR with renal flow ---------------------------------------------
inits <- parameterize_pbtk(chem.cas=NULL,chem.name = "Chlorpyrifos", species = "Human", default.to.human = F,
                            tissuelist = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                            force.human.clint.fub = F, clint.pvalue.threshold = 0.05, use.qrenal = TRUE)


