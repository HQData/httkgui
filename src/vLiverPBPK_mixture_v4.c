/* 
 This is modification of vLiverPBPK .c program with 
 addition of new metabolism parameters
 as well as option to run both compartments at once
 
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
#define C2_ID_Agutlumen    0x000b
#define C2_ID_Agut         0x000c
#define C2_ID_Aliver       0x000d
#define C2_ID_Aven         0x000e
#define C2_ID_Alung        0x000f
#define C2_ID_Aart         0x0010
#define C2_ID_Arest        0x0011
#define C2_ID_Akidney      0x0012
#define C2_ID_Atubules     0x0013
#define C2_ID_Ametabolized 0x0014
#define C2_ID_AUC          0x0015

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
#define C2_ID_Cgut         0x0009
#define C2_ID_Cliver       0x000a
#define C2_ID_Cven         0x000b
#define C2_ID_Clung        0x000c
#define C2_ID_Cart         0x000d
#define C2_ID_Crest        0x000e
#define C2_ID_Ckidney      0x000f
#define C2_ID_Cserum       0x0010
#define C2_ID_Aserum       0x0011

static double parms[102];

/* Parameters: first compound */

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

/* Parameters: second compound  */
#define C2_BW parms[45]
#define C2_CLmetabolismc parms[46]
#define C2_hematocrit parms[47]
#define C2_kgutabs parms[48]
#define C2_Kkidney2plasma parms[49]
#define C2_Kliver2plasma parms[50]
#define C2_Krest2plasma parms[51]
#define C2_Kgut2plasma parms[52]
#define C2_Klung2plasma parms[53]
#define C2_Qcardiacc parms[54]
#define C2_Qgfrc parms[55]
#define C2_Qgutf parms[56]
#define C2_Qkidneyf parms[57]
#define C2_Qliverf parms[58]
#define C2_Vartc parms[59]
#define C2_Vgutc parms[60]
#define C2_Vkidneyc parms[61]
#define C2_Vliverc parms[62]
#define C2_Vlungc parms[63]
#define C2_Vrestc parms[64]
#define C2_Vvenc parms[65]
#define C2_Fraction_unbound_plasma parms[66]
#define C2_Ratioblood2plasma parms[67]
#define C2_CLmetabolism parms[68]
#define C2_Qcardiac parms[69]
#define C2_Qgfr parms[70]
#define C2_Qgut parms[71]
#define C2_Qkidney parms[72]
#define C2_Qliver parms[73]
#define C2_Qrest parms[74]
#define C2_Vart parms[75]
#define C2_Vgut parms[76]
#define C2_Vkidney parms[77]
#define C2_Vliver parms[78]
#define C2_Vlung parms[79]
#define C2_Vrest parms[80]
#define C2_Vven parms[81]
#define C2_Vmax_liver parms[82]
#define C2_Km_liver parms[83]
#define C2_Vmax_gut parms[84]
#define C2_Km_gut parms[85]
#define C2_Vmax_kidney parms[86]
#define C2_Km_kidney parms[87]
#define C2_CLmetabolism_gut parms[88]
#define C2_CLmetabolism_kidney parms[89]

/* Parameters of interaction */
#define C1_alpha_liver  parms[90]
#define C1_ki_liver     parms[91]
#define C1_alpha_gut    parms[92]
#define C1_ki_gut       parms[93]
#define C1_alpha_kidney parms[94]
#define C1_ki_kidney    parms[95]
#define C2_alpha_liver  parms[96]
#define C2_ki_liver     parms[97]
#define C2_alpha_gut    parms[98]
#define C2_ki_gut       parms[99]
#define C2_alpha_kidney parms[100]
#define C2_ki_kidney    parms[101]


/*----- Initializers */ 
void initmod_mixture (void (* odeparms)(int *, double *))
{
  int N=102;
  odeparms(&N, parms);
}

/*----- Dynamics section */

void derivs_mixture (int *neq, double *pdTime, double *y, double *ydot, double *yout, int *ip)
{
  double ametabolised_liver; double ametabolised_kidney; double ametabolised_gut;
  double prop_liver, prop_kidney, prop_gut; //for mixture part
  double C1_C_liver, C1_C_kidney, C1_C_gut; //available concentrations
  double C2_C_liver, C2_C_kidney, C2_C_gut; //available concentrations
  
  /* Available concentrations */
  
  //concentrations (amount/volume)
  yout[C1_ID_Cgut] = y[C1_ID_Agut] / C1_Vgut ;
  yout[C1_ID_Cliver] = y[C1_ID_Aliver] / C1_Vliver ;
  yout[C1_ID_Cven] = y[C1_ID_Aven] / C1_Vven ;
  yout[C1_ID_Clung] = y[C1_ID_Alung] / C1_Vlung ;
  yout[C1_ID_Cart] = y[C1_ID_Aart] / C1_Vart ;
  yout[C1_ID_Crest] = y[C1_ID_Arest] / C1_Vrest ;
  yout[C1_ID_Ckidney] = y[C1_ID_Akidney] / C1_Vkidney ;
  yout[C1_ID_Cserum] = y[C1_ID_Aven] / C1_Vven / C1_Ratioblood2plasma ; 
  yout[C1_ID_Aserum] = y[C1_ID_Aven] / C1_Ratioblood2plasma * ( 1 - C1_hematocrit ) ;
  ydot[C1_ID_Agutlumen] = - C1_kgutabs * y[C1_ID_Agutlumen] ;
  yout[C2_ID_Cgut] = y[C2_ID_Agut] / C2_Vgut ;
  yout[C2_ID_Cliver] = y[C2_ID_Aliver] / C2_Vliver ;
  yout[C2_ID_Cven] = y[C2_ID_Aven] / C2_Vven ;
  yout[C2_ID_Clung] = y[C2_ID_Alung] / C2_Vlung ;
  yout[C2_ID_Cart] = y[C2_ID_Aart] / C2_Vart ;
  yout[C2_ID_Crest] = y[C2_ID_Arest] / C2_Vrest ;
  yout[C2_ID_Ckidney] = y[C2_ID_Akidney] / C2_Vkidney ;
  yout[C2_ID_Cserum] = y[C2_ID_Aven] / C2_Vven / C2_Ratioblood2plasma ; 
  yout[C2_ID_Aserum] = y[C2_ID_Aven] / C2_Ratioblood2plasma * ( 1 - C2_hematocrit ) ;
  ydot[C2_ID_Agutlumen] = - C2_kgutabs * y[C2_ID_Agutlumen] ;
  
  //available concentrations
  C1_C_liver  = (yout[C1_ID_Cliver]  / C1_Kliver2plasma );
  C1_C_kidney = (yout[C1_ID_Ckidney] / C1_Kkidney2plasma);
  C1_C_gut    = (yout[C1_ID_Cgut]    / C1_Kgut2plasma   );
  C2_C_liver  = (yout[C2_ID_Cliver]  / C2_Kliver2plasma );
  C2_C_kidney = (yout[C2_ID_Ckidney] / C2_Kkidney2plasma);
  C2_C_gut    = (yout[C2_ID_Cgut]    / C2_Kgut2plasma   );
  
  
  
  /* Compund 1 */
  
  //start flows with calculation of amount metabolised in each of 3 organs: liver, kidney and gut
  //it will depend on whether Km is provided (as != 0)
  if(C1_Km_liver != 0) {
    if((C2_alpha_liver == 0) | (C2_ki_liver == 0)) {
      ametabolised_liver  = (C1_Vmax_liver * y[C1_ID_Aliver] / C1_Vliver / C1_Kliver2plasma) / 
        (C1_Km_liver + y[C1_ID_Aliver] / C1_Vliver / C1_Kliver2plasma);
    } else {
      prop_liver = 1 + (C2_C_liver / C2_alpha_liver / C2_ki_liver);
      ametabolised_liver = 
        (C1_Vmax_liver / prop_liver) * C1_C_liver /                               // Vmax scaled by some factor * C available
        ((((1 + ((yout[C2_ID_Cliver]/C2_Kliver2plasma)/C2_ki_liver))/prop_liver)  // Km scaled ...
            * C1_Km_liver) + C1_C_liver);                                         // ...times Km + C available
    }
  } else 
    ametabolised_liver  = (C1_CLmetabolism        * y[C1_ID_Aliver] / C1_Vliver / C1_Kliver2plasma);
  
  if(C1_Km_kidney != 0) {
    if((C2_alpha_kidney == 0) | (C2_ki_kidney == 0)) {
      ametabolised_kidney  = (C1_Vmax_kidney * y[C1_ID_Akidney] / C1_Vkidney / C1_Kkidney2plasma) / 
        (C1_Km_kidney + y[C1_ID_Akidney] / C1_Vkidney / C1_Kkidney2plasma);
    } else {
      prop_kidney = 1 + (C2_C_kidney / C2_alpha_kidney / C2_ki_kidney);
      ametabolised_kidney = 
        (C1_Vmax_kidney / prop_kidney) * C1_C_kidney /                               // Vmax scaled by some factor * C available
        ((((1 + ((yout[C2_ID_Ckidney]/C2_Kkidney2plasma)/C2_ki_kidney))/prop_kidney)  // Km scaled ...
            * C1_Km_kidney) + C1_C_kidney);                                         // ...times Km + C available
    }
  } else
    ametabolised_kidney = (C1_CLmetabolism_kidney * y[C1_ID_Akidney] / C1_Vkidney / C1_Kkidney2plasma);
  
  if(C1_Km_gut != 0) {
    if((C2_alpha_gut == 0) | (C2_ki_gut == 0)) {
      ametabolised_gut  = (C1_Vmax_gut * y[C1_ID_Agut] / C1_Vgut / C1_Kgut2plasma) / 
        (C1_Km_gut + y[C1_ID_Agut] / C1_Vgut / C1_Kgut2plasma);
    } else {
      prop_gut = 1 + (C2_C_gut / C2_alpha_gut / C2_ki_gut);
      ametabolised_gut = 
        (C1_Vmax_gut / prop_gut) * C1_C_gut /                               // Vmax scaled by some factor * C available
        ((((1 + ((yout[C2_ID_Cgut]/C2_Kgut2plasma)/C2_ki_gut))/prop_gut)  // Km scaled ...
            * C1_Km_gut) + C1_C_gut);                                         // ...times Km + C available
    }
  } else
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
  
  
  
  /* Compound 2 */
  

  //start flows with calculation of amount metabolised in each of 3 organs: liver, kidney and gut
  //it will depend on whether Km is provided (as != 0)
  
  if(C2_Km_liver != 0) {
    if((C2_alpha_liver == 0) | (C2_ki_liver == 0)) {
      ametabolised_liver  = (C2_Vmax_liver * y[C2_ID_Aliver] / C2_Vliver / C2_Kliver2plasma) / 
        (C2_Km_liver + y[C2_ID_Aliver] / C2_Vliver / C2_Kliver2plasma);
    } else {
      prop_liver = 1 + (C2_C_liver / C2_alpha_liver / C2_ki_liver);
      ametabolised_liver = 
        (C2_Vmax_liver / prop_liver) * C2_C_liver /                               // Vmax scaled by some factor * C available
        ((((1 + ((yout[C2_ID_Cliver]/C2_Kliver2plasma)/C2_ki_liver))/prop_liver)  // Km scaled ...
            * C2_Km_liver) + C2_C_liver);                                         // ...times Km + C available
    }
  } else 
    ametabolised_liver  = (C2_CLmetabolism        * y[C2_ID_Aliver] / C2_Vliver / C2_Kliver2plasma);
  
  if(C2_Km_kidney != 0) {
    if((C2_alpha_kidney == 0) | (C2_ki_kidney == 0)) {
      ametabolised_kidney  = (C2_Vmax_kidney * y[C2_ID_Akidney] / C2_Vkidney / C2_Kkidney2plasma) / 
        (C2_Km_kidney + y[C2_ID_Akidney] / C2_Vkidney / C2_Kkidney2plasma);
    } else {
      prop_kidney = 1 + (C2_C_kidney / C2_alpha_kidney / C2_ki_kidney);
      ametabolised_kidney = 
        (C2_Vmax_kidney / prop_kidney) * C2_C_kidney /                               // Vmax scaled by some factor * C available
        ((((1 + ((yout[C2_ID_Ckidney]/C2_Kkidney2plasma)/C2_ki_kidney))/prop_kidney)  // Km scaled ...
            * C2_Km_kidney) + C2_C_kidney);                                         // ...times Km + C available
    }
  } else
    ametabolised_kidney = (C2_CLmetabolism_kidney * y[C2_ID_Akidney] / C2_Vkidney / C2_Kkidney2plasma);
  
  if(C2_Km_gut != 0) {
    if((C2_alpha_gut == 0) | (C2_ki_gut == 0)) {
      ametabolised_gut  = (C2_Vmax_gut * y[C2_ID_Agut] / C2_Vgut / C2_Kgut2plasma) / 
        (C2_Km_gut + y[C2_ID_Agut] / C2_Vgut / C2_Kgut2plasma);
    } else {
      prop_gut = 1 + (C2_C_gut / C2_alpha_gut / C2_ki_gut);
      ametabolised_gut = 
        (C2_Vmax_gut / prop_gut) * C2_C_gut /                               // Vmax scaled by some factor * C available
        ((((1 + ((yout[C2_ID_Cgut]/C2_Kgut2plasma)/C2_ki_gut))/prop_gut)  // Km scaled ...
            * C2_Km_gut) + C2_C_gut);                                         // ...times Km + C available
    }
  } else
    ametabolised_gut = (C2_CLmetabolism_gut    * y[C2_ID_Agut] / C2_Vgut / C2_Kgut2plasma);
  
  
  /* Agut flow: modified by subtracting term gut-specific metabolism */
  ydot[C2_ID_Agut] = C2_kgutabs * y[C2_ID_Agutlumen] + C2_Qgut * ( y[C2_ID_Aart] / C2_Vart - y[C2_ID_Agut] / C2_Vgut * C2_Ratioblood2plasma / C2_Kgut2plasma / C2_Fraction_unbound_plasma );
  ydot[C2_ID_Agut] -= ametabolised_gut;
  
  /* Aliver flow: modified by subtracting liver-specific metabolism */
  ydot[C2_ID_Aliver] = C2_Qliver * y[C2_ID_Aart] / C2_Vart + C2_Qgut * y[C2_ID_Agut] / C2_Vgut * C2_Ratioblood2plasma / C2_Kgut2plasma / C2_Fraction_unbound_plasma - 
  ( C2_Qliver + C2_Qgut ) * y[C2_ID_Aliver] / C2_Vliver / C2_Kliver2plasma / C2_Fraction_unbound_plasma * C2_Ratioblood2plasma;
  ydot[C2_ID_Aliver] -= ametabolised_liver;
  
  
  ydot[C2_ID_Aven]  = ( ( C2_Qliver + C2_Qgut ) * y[C2_ID_Aliver] / C2_Vliver / C2_Kliver2plasma + C2_Qkidney * y[C2_ID_Akidney] / C2_Vkidney / C2_Kkidney2plasma + 
    C2_Qrest * y[C2_ID_Arest] / C2_Vrest / C2_Krest2plasma ) * C2_Ratioblood2plasma / C2_Fraction_unbound_plasma - C2_Qcardiac * y[C2_ID_Aven] / C2_Vven ;
  ydot[C2_ID_Alung] = C2_Qcardiac * ( y[C2_ID_Aven]  / C2_Vven - y[C2_ID_Alung] / C2_Vlung * C2_Ratioblood2plasma / C2_Klung2plasma / C2_Fraction_unbound_plasma ) ;
  ydot[C2_ID_Aart]  = C2_Qcardiac * ( y[C2_ID_Alung] / C2_Vlung * C2_Ratioblood2plasma / C2_Klung2plasma / C2_Fraction_unbound_plasma - y[C2_ID_Aart] / C2_Vart ) ;
  ydot[C2_ID_Arest] = C2_Qrest    * ( y[C2_ID_Aart]  / C2_Vart - y[C2_ID_Arest] / C2_Vrest * C2_Ratioblood2plasma / C2_Krest2plasma / C2_Fraction_unbound_plasma ) ;
  
  /* Akidney flow: modified by subtracting Qgfr and liver-specific metabolism */
  ydot[C2_ID_Akidney]  = C2_Qkidney * y[C2_ID_Aart] / C2_Vart - C2_Qkidney * y[C2_ID_Akidney] / C2_Vkidney / C2_Kkidney2plasma * C2_Ratioblood2plasma / C2_Fraction_unbound_plasma;
  ydot[C2_ID_Akidney] -= C2_Qgfr * y[C2_ID_Akidney] / C2_Vkidney / C2_Kkidney2plasma; //see Atubules
  ydot[C2_ID_Akidney] -= ametabolised_kidney;
  /* End Akidney flow. */
  
  ydot[C2_ID_Atubules] = C2_Qgfr * y[C2_ID_Akidney] / C2_Vkidney / C2_Kkidney2plasma ;
  
  /* Modification: contribution to amount metabolised partitioned into three organs */
  ydot[C2_ID_Ametabolized] = ametabolised_liver + ametabolised_gut + ametabolised_kidney;
  
  ydot[C2_ID_AUC] = y[C2_ID_Aven] / C2_Vven / C2_Ratioblood2plasma ; 
  
} /* derivs */
  