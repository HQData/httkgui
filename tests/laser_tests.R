#### httk tests ####

# parametrization test (we have added Vmax, km) -----------------------------------

inits2 <- parameterize_pbtk(chem.cas=NULL,chem.name = "Chlorpyrifos", species = "Human", default.to.human = F,
                           tissuelist = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                           force.human.clint.fub = F, clint.pvalue.threshold = 0.05)



# variability on GFR - test -----------------------------------------------
pbtk_mclist <- list("inits"=list(), "pbtk_result"=list())
for(i in 1:10) {
    inits <- parameterize_pbtk(chem.cas=NULL,chem.name = "Bisphenol-A", species = "Human", default.to.human = F,
                               tissuelist = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                               force.human.clint.fub = F, clint.pvalue.threshold = 0.05, monte.carlo=TRUE)
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

