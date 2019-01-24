<!-- README.md is generated from README.Rmd. Please edit that file -->
httkgui
=======

`httkgui` was a proof of concept modification of a popular CRAN toxicokinetics package `httk` (v1.4); `httkgui` modifies some of the TK models, allows rudimentary Monte Carlo simulations of said models and provides GUI ("TKPlate") to run them and review the results

The [original 'httk' (v1.8 as of Jan 2019) package](https://cran.r-project.org/web/packages/httk/index.html) is by John Wambaugh. The code modifications to 'httk' have been developed by Witold Wiecek with concept by Nadia Quignot and Billy Amzal at Analytica Laser (now Certara) as a part of an European Food Safety Authority project *“Integrating toxicokinetics in chemical risk assessment: application to human, animal and environmental risk assessment”* (OC/EFSA/SCER/2014/06)

Installation and use
--------------------

Core functions are same as in `httk`. Changes to models are not documented within code but version controlled.

``` r
devtools::install_github("wwiecek/httkgui")
library(httkgui)
pbtkUI() #Proof-of-concept GUI
```
