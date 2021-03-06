initparms <- function(newParms = NULL){
  parms <- c(
    BW = 70,
    CLmetabolismc = 0,
    hematocrit = 0.44,
    kgutabs = 1,
    Kkidney2plasma = 0,
    Kliver2plasma = 0,
    Krest2plasma = 0,
    Kgut2plasma = 0,
    Klung2plasma = 0,
    Qcardiacc = 4.8,
    Qgfrc = 0,
    Qgutf = 0.205,
    Qkidneyf = 0.221,
    Qliverf = 0.0536,
    Vartc = 0.0487,
    Vgutc = 0.0158,
    Vkidneyc = 0.00119,
    Vliverc = 0.02448,
    Vlungc = 0.00723,
    Vrestc = 0.77654,
    Vvenc = 0.0487,
    Fraction_unbound_plasma = 0.0682,
    Ratioblood2plasma = 0.0,
    CLmetabolism = 0.0,
    Qcardiac = 0.0,
    Qgfr = 0.0,
    Qgut = 0.0,
    Qkidney = 0.0,
    Qliver = 0.0,
    Qrest = 0.0,
    Vart = 0.0,
    Vgut = 0.0,
    Vkidney = 0.0,
    Vliver = 0.0,
    Vlung = 0.0,
    Vrest = 0.0,
    Vven = 0.0,
    Vmax = 0.0,
    km = 0.0
  )
  if (!is.null(newParms)) {
    if (!all(names(newParms) %in% c(names(parms)))) {
      stop("illegal parameter name")
    }
  }
  if (!is.null(newParms)) parms[names(newParms)] <- newParms
  out <- .C("getParms",
   as.double(parms),
  out=double(length(parms)),
  as.integer(length(parms)))$out
  names(out) <- names(parms)
  out
}

#variant with extra Km and Vmax parameters:
initparms_3cl <- function(newParms = NULL){
  parms <- c(
    BW = 70,
    CLmetabolismc = 0,
    hematocrit = 0.44,
    kgutabs = 1,
    Kkidney2plasma = 0,
    Kliver2plasma = 0,
    Krest2plasma = 0,
    Kgut2plasma = 0,
    Klung2plasma = 0,
    Qcardiacc = 4.8,
    Qgfrc = 0,
    Qgutf = 0.205,
    Qkidneyf = 0.221,
    Qliverf = 0.0536,
    Vartc = 0.0487,
    Vgutc = 0.0158,
    Vkidneyc = 0.00119,
    Vliverc = 0.02448,
    Vlungc = 0.00723,
    Vrestc = 0.77654,
    Vvenc = 0.0487,
    Fraction_unbound_plasma = 0.0682,
    Ratioblood2plasma = 0.0,
    CLmetabolism = 0.0,
    CLmetabolism_gut = 0.0,
    CLmetabolism_kidney = 0.0,
    Qcardiac = 0.0,
    Qgfr = 0.0,
    Qgut = 0.0,
    Qkidney = 0.0,
    Qliver = 0.0,
    Qrest = 0.0,
    Vart = 0.0,
    Vgut = 0.0,
    Vkidney = 0.0,
    Vliver = 0.0,
    Vlung = 0.0,
    Vrest = 0.0,
    Vven = 0.0,
    Vmax = 0.0,
    km = 0.0,
    Vmax_liver = 0.0,
    km_liver = 0.0,
    Vmax_gut = 0.0,
    km_liver = 0.0
  )
  if (!is.null(newParms)) {
    if (!all(names(newParms) %in% c(names(parms)))) {
      browser()
      stop("illegal parameter name")
    }
  }
  
  if (!is.null(newParms)) parms[names(newParms)] <- newParms
  out <- .C("getParms_3cl",
    as.double(parms),
    out=double(length(parms)),
    as.integer(length(parms)))$out
  names(out) <- names(parms)
  out
}

#new function to do model scaling
model_scaling_3cl <- function(newParms = NULL) {
  parms <- c(
    BW = 70,
    CLmetabolismc = 0,
    hematocrit = 0.44,
    kgutabs = 1,
    Kkidney2plasma = 0,
    Kliver2plasma = 0,
    Krest2plasma = 0,
    Kgut2plasma = 0,
    Klung2plasma = 0,
    Qcardiacc = 4.8,
    Qgfrc = 0,
    Qgutf = 0.205,
    Qkidneyf = 0.221,
    Qliverf = 0.0536,
    Vartc = 0.0487,
    Vgutc = 0.0158,
    Vkidneyc = 0.00119,
    Vliverc = 0.02448,
    Vlungc = 0.00723,
    Vrestc = 0.77654,
    Vvenc = 0.0487,
    Fraction_unbound_plasma = 0.0682,
    Ratioblood2plasma = 0.0,
    CLmetabolism = 0.0,
    Qcardiac = 0.0,
    Qgfr = 0.0,
    Qgut = 0.0,
    Qkidney = 0.0,
    Qliver = 0.0,
    Qrest = 0.0,
    Vart = 0.0,
    Vgut = 0.0,
    Vkidney = 0.0,
    Vliver = 0.0,
    Vlung = 0.0,
    Vrest = 0.0,
    Vven = 0.0,
    Vmax = 0.0,
    km = 0.0,
    Vmax_gut = 0.0,
    km_gut = 0.0,
    Vmax_kidney = 0.0,
    km_kidney = 0.0,
    CLmetabolism_gut = 0.0,
    CLmetabolism_kidney = 0.0
  )
  if (!is.null(newParms)) {
    if (!all(names(newParms) %in% c(names(parms)))) {
      stop("illegal parameter name")
    }
  }
  if (!is.null(newParms)) parms[names(newParms)] <- newParms
  parms <- within(parms, {
    kgutabs = kgutabs * 24 ;
    CLmetabolism = CLmetabolismc * 24 * BW ;
    CLmetabolism_gut = CLmetabolism_gut * 24 * BW ;
    CLmetabolism_kidney = CLmetabolism_kidney * 24 * BW ;
    Vmax = Vmax * 24 * BW ;
    Vmax_gut = Vmax_gut * 24 * BW ;
    Vmax_kidney = Vmax_kidney * 24 * BW ;
    Qcardiac = Qcardiacc * 24 * BW^0.75 ;
    Qgfr = Qgfrc * BW^0.75 * 24 ;
    Qgut = Qcardiac * Qgutf ;
    Qkidney = Qcardiac * Qkidneyf ;
    Qliver = Qcardiac * Qliverf ;
    Qrest = Qcardiac - ( Qgut + Qkidney + Qliver ) ;
    Vart = Vartc * BW ;
    Vgut = Vgutc * BW ;
    Vkidney = Vkidneyc * BW;
    Vliver = Vliverc * BW ;
    Vlung = Vlungc * BW ;
    Vrest = Vrestc * BW ;
    Vven = Vvenc * BW ;
  })
  parms_numeric <- as.numeric(parms)
  names(parms_numeric) <- names(parms)
  return(parms_numeric)
}

Outputs <- c(
    "Cgut",
    "Cliver",
    "Cven",
    "Clung",
    "Cart",
    "Crest",
    "Ckidney",
    "Cserum",
    "Aserum"
)


initState <- function(parms, newState = NULL) {
  Y <- c(
    Agutlumen = 0.0,
    Agut = 0.0,
    Aliver = 0.0,
    Aven = 0.0,
    Alung = 0.0,
    Aart = 0.0,
    Arest = 0.0,
    Akidney = 0.0,
    Atubules = 0.0,
    Ametabolized = 0.0,
    AUC = 0.0
  )
  Y <- with(as.list(parms), {  Y
  })

  if (!is.null(newState)) {
    if (!all(names(newState) %in% c(names(Y)))) {
      stop("illegal state variable name in newState")
    }
    Y[names(newState)] <- newState
  }
  Y
}
