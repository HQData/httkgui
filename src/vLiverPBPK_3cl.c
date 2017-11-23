/* 
 This is modification of vLiverPBPK .c program with addition of some parameters
 
 New parameters are compartment-specific CL parameters or Vmax, km.
 If CL is to be used, km should be supplied as 0. If km != 0 we automatically use Vmax and km.
 
 Code has been cleaned up and reorganised.
 C1_ prefix is added for convenience to expand the script to mixtures (see separate script).
 
 Author of the modification: WW / LASER
*/

#include <R.h>

/* Model variables: States */
#define C1_ID_Agutlumen    0x0000
#define C1_ID_Agut         0x0001
#define C1_ID_Aliver       0x0002
#define C1_ID_Aven         0x0003
#define C1_ID_Alung        0x0004
#define C1_ID_Aart         0x0005
#define C1_ID_Arest        0x0006
#define C1_ID_Akidney      0x0007
#define C1_ID_Atubules     0x0008
#define C1_ID_Ametabolized 0x0009
#define C1_ID_AUC          0x000a

/* Model variables: Outputs */
#define C1_ID_Cgut         0x0000
#define C1_ID_Cliver       0x0001
#define C1_ID_Cven         0x0002
#define C1_ID_Clung        0x0003
#define C1_ID_Cart         0x0004
#define C1_ID_Crest        0x0005
#define C1_ID_Ckidney      0x0006
#define C1_ID_Cserum       0x0007
#define C1_ID_Aserum       0x0008


static double parms[45];

/* Parameters */
// note that C1_CLmetabolismc is either a stupid mess or I'm just not getting something!
// but it is kept here for compatibility purposes
#define C1_BW parms[0]
#define C1_CLmetabolismc parms[1]
#define C1_hematocrit parms[2]
#define C1_kgutabs parms[3]
#define C1_Kkidney2plasma parms[4]
#define C1_Kliver2plasma parms[5]
#define C1_Krest2plasma parms[6]
#define C1_Kgut2plasma parms[7]
#define C1_Klung2plasma parms[8]
#define C1_Qcardiacc parms[9]
#define C1_Qgfrc parms[10]
#define C1_Qgutf parms[11]
#define C1_Qkidneyf parms[12]
#define C1_Qliverf parms[13]
#define C1_Vartc parms[14]
#define C1_Vgutc parms[15]
#define C1_Vkidneyc parms[16]
#define C1_Vliverc parms[17]
#define C1_Vlungc parms[18]
#define C1_Vrestc parms[19]
#define C1_Vvenc parms[20]
#define C1_Fraction_unbound_plasma parms[21]
#define C1_Ratioblood2plasma parms[22]
#define C1_CLmetabolism parms[23]
#define C1_Qcardiac parms[24]
#define C1_Qgfr parms[25]
#define C1_Qgut parms[26]
#define C1_Qkidney parms[27]
#define C1_Qliver parms[28]
#define C1_Qrest parms[29]
#define C1_Vart parms[30]
#define C1_Vgut parms[31]
#define C1_Vkidney parms[32]
#define C1_Vliver parms[33]
#define C1_Vlung parms[34]
#define C1_Vrest parms[35]
#define C1_Vven parms[36]
#define C1_Vmax_liver parms[37]
#define C1_Km_liver parms[38]
#define C1_Vmax_gut parms[39]
#define C1_Km_gut parms[40]
#define C1_Vmax_kidney parms[41]
#define C1_Km_kidney parms[42]
#define C1_CLmetabolism_gut parms[43]
#define C1_CLmetabolism_kidney parms[44]


/*----- Initializers */ 
void initmod_3cl (void (* odeparms)(int *, double *))
{
  int N=45;
  odeparms(&N, parms);
}

/*----- Dynamics section */

void derivs_3cl (int *neq, double *pdTime, double *y, double *ydot, double *yout, int *ip)
{
  double ametabolised_liver; double ametabolised_kidney; double ametabolised_gut;
  
  //concentrations (amount/volume)
  yout[C1_ID_Cgut] = y[C1_ID_Agut] / C1_Vgut ;
  yout[C1_ID_Cliver] = y[C1_ID_Aliver] / C1_Vliver ;
  yout[C1_ID_Cven] = y[C1_ID_Aven] / C1_Vven ;
  yout[C1_ID_Clung] = y[C1_ID_Alung] / C1_Vlung ;
  yout[C1_ID_Cart] = y[C1_ID_Aart] / C1_Vart ;
  yout[C1_ID_Crest] = y[C1_ID_Arest] / C1_Vrest ;
  yout[C1_ID_Ckidney] = y[C1_ID_Akidney] / C1_Vkidney ;
  yout[C1_ID_Cserum] = y[C1_ID_Aven] / C1_Vven / C1_Ratioblood2plasma ; 
  //(ratioblood2serum the same as for plasma)
  yout[C1_ID_Aserum] = y[C1_ID_Aven] / C1_Ratioblood2plasma * ( 1 - C1_hematocrit ) ;
  ydot[C1_ID_Agutlumen] = - C1_kgutabs * y[C1_ID_Agutlumen] ;
  
  //start flows with calculation of amount metabolised in each of 3 organs: liver, kidney and gut
  //it will depend on whether Km is provided (as != 0)
  
  if(C1_Km_liver != 0) 
    ametabolised_liver  = (C1_Vmax_liver * y[C1_ID_Aliver] / C1_Vliver / C1_Kliver2plasma) / 
      (C1_Km_liver + y[C1_ID_Aliver] / C1_Vliver / C1_Kliver2plasma);
  else 
    ametabolised_liver  = (C1_CLmetabolism        * y[C1_ID_Aliver] / C1_Vliver / C1_Kliver2plasma);
  
  if(C1_Km_kidney!=0)
    ametabolised_kidney = (C1_Vmax_kidney * y[C1_ID_Akidney] / C1_Vkidney / C1_Kkidney2plasma) / 
      (C1_Km_kidney + y[C1_ID_Akidney] / C1_Vkidney / C1_Kkidney2plasma);
  else
    ametabolised_kidney = (C1_CLmetabolism_kidney * y[C1_ID_Akidney] / C1_Vkidney / C1_Kkidney2plasma);
  
  if(C1_Km_gut != 0)
    ametabolised_gut = (C1_Vmax_gut * y[C1_ID_Agut] / C1_Vgut / C1_Kgut2plasma) / 
      (C1_Km_gut + y[C1_ID_Agut] / C1_Vgut / C1_Kgut2plasma);
  else
    ametabolised_gut = (C1_CLmetabolism_gut    * y[C1_ID_Agut] / C1_Vgut / C1_Kgut2plasma);
  
  /* Agut flow: modified by subtracting term gut-specific metabolism */
  ydot[C1_ID_Agut] = C1_kgutabs * y[C1_ID_Agutlumen] + C1_Qgut * ( y[C1_ID_Aart] / C1_Vart - y[C1_ID_Agut] / C1_Vgut * C1_Ratioblood2plasma / C1_Kgut2plasma / C1_Fraction_unbound_plasma );
  ydot[C1_ID_Agut] -= ametabolised_gut;
  
  /* Aliver flow: modified by subtracting liver-specific metabolism */
  ydot[C1_ID_Aliver] = C1_Qliver * y[C1_ID_Aart] / C1_Vart + C1_Qgut * y[C1_ID_Agut] / C1_Vgut * C1_Ratioblood2plasma / C1_Kgut2plasma / C1_Fraction_unbound_plasma - 
  ( C1_Qliver + C1_Qgut ) * y[C1_ID_Aliver] / C1_Vliver / C1_Kliver2plasma / C1_Fraction_unbound_plasma * C1_Ratioblood2plasma;
  ydot[C1_ID_Aliver] -= ametabolised_liver;
  
  
  ydot[C1_ID_Aven]  = ( ( C1_Qliver + C1_Qgut ) * y[C1_ID_Aliver] / C1_Vliver / C1_Kliver2plasma + C1_Qkidney * y[C1_ID_Akidney] / C1_Vkidney / C1_Kkidney2plasma + 
    C1_Qrest * y[C1_ID_Arest] / C1_Vrest / C1_Krest2plasma ) * C1_Ratioblood2plasma / C1_Fraction_unbound_plasma - C1_Qcardiac * y[C1_ID_Aven] / C1_Vven ;
  ydot[C1_ID_Alung] = C1_Qcardiac * ( y[C1_ID_Aven]  / C1_Vven - y[C1_ID_Alung] / C1_Vlung * C1_Ratioblood2plasma / C1_Klung2plasma / C1_Fraction_unbound_plasma ) ;
  ydot[C1_ID_Aart]  = C1_Qcardiac * ( y[C1_ID_Alung] / C1_Vlung * C1_Ratioblood2plasma / C1_Klung2plasma / C1_Fraction_unbound_plasma - y[C1_ID_Aart] / C1_Vart ) ;
  ydot[C1_ID_Arest] = C1_Qrest    * ( y[C1_ID_Aart]  / C1_Vart - y[C1_ID_Arest] / C1_Vrest * C1_Ratioblood2plasma / C1_Krest2plasma / C1_Fraction_unbound_plasma ) ;
  
  /* Akidney flow: modified by subtracting Qgfr and liver-specific metabolism */
  ydot[C1_ID_Akidney]  = C1_Qkidney * y[C1_ID_Aart] / C1_Vart - C1_Qkidney * y[C1_ID_Akidney] / C1_Vkidney / C1_Kkidney2plasma * C1_Ratioblood2plasma / C1_Fraction_unbound_plasma;
  ydot[C1_ID_Akidney] -= C1_Qgfr * y[C1_ID_Akidney] / C1_Vkidney / C1_Kkidney2plasma; //see Atubules
  ydot[C1_ID_Akidney] -= ametabolised_kidney;
  /* End Akidney flow. */
  
  ydot[C1_ID_Atubules] = C1_Qgfr * y[C1_ID_Akidney] / C1_Vkidney / C1_Kkidney2plasma ;
  
  /* Modification: contribution to amount metabolised partitioned into three organs */
  ydot[C1_ID_Ametabolized] = ametabolised_liver + ametabolised_gut + ametabolised_kidney;
  
  ydot[C1_ID_AUC] = y[C1_ID_Aven] / C1_Vven / C1_Ratioblood2plasma ; 
  
} /* derivs */
