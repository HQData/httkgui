/* 
 This is modification of vLiverPBPK .c program with addition of new metabolism parameters
 
 New parameters are compartment-specific CL parameter, Vmax, km.
 * If Vmax & km are available (by passing km != 0) we use them for calculation
 * Otherwise CL parameter used
 
 Author of the modification: WW / LASER
*/

#include <R.h>

/* Model variables: States */
#define ID_Agutlumen 0x0000
#define ID_Agut 0x0001
#define ID_Aliver 0x0002
#define ID_Aven 0x0003
#define ID_Alung 0x0004
#define ID_Aart 0x0005
#define ID_Arest 0x0006
#define ID_Akidney 0x0007
#define ID_Atubules 0x0008
#define ID_Ametabolized 0x0009
#define ID_AUC 0x000a

/* Model variables: Outputs */
#define ID_Cgut 0x0000
#define ID_Cliver 0x0001
#define ID_Cven 0x0002
#define ID_Clung 0x0003
#define ID_Cart 0x0004
#define ID_Crest 0x0005
#define ID_Ckidney 0x0006
#define ID_Cserum 0x0007
#define ID_Aserum 0x0008

/* Parameters  */
static double parms[45];

#define BW parms[0]
#define CLmetabolismc parms[1]
#define hematocrit parms[2]
#define kgutabs parms[3]
#define Kkidney2plasma parms[4]
#define Kliver2plasma parms[5]
#define Krest2plasma parms[6]
#define Kgut2plasma parms[7]
#define Klung2plasma parms[8]
#define Qcardiacc parms[9]
#define Qgfrc parms[10]
#define Qgutf parms[11]
#define Qkidneyf parms[12]
#define Qliverf parms[13]
#define Vartc parms[14]
#define Vgutc parms[15]
#define Vkidneyc parms[16]
#define Vliverc parms[17]
#define Vlungc parms[18]
#define Vrestc parms[19]
#define Vvenc parms[20]
#define Fraction_unbound_plasma parms[21]
#define Ratioblood2plasma parms[22]
#define CLmetabolism parms[23]
#define Qcardiac parms[24]
#define Qgfr parms[25]
#define Qgut parms[26]
#define Qkidney parms[27]
#define Qliver parms[28]
#define Qrest parms[29]
#define Vart parms[30]
#define Vgut parms[31]
#define Vkidney parms[32]
#define Vliver parms[33]
#define Vlung parms[34]
#define Vrest parms[35]
#define Vven parms[36]

//new additions (LASER):
#define Vmax_liver parms[37]
#define Km_liver parms[38]
#define Vmax_gut parms[39]
#define Km_gut parms[40]
#define Vmax_kidney parms[41]
#define Km_kidney parms[42]
#define CLmetabolism_gut parms[43]
#define CLmetabolism_kidney parms[44]

double lastterm; //for calculating this 'last term to subtract' 

/*----- Initializers */ 
void initmod_3cl (void (* odeparms)(int *, double *))
{
  int N=45;
  odeparms(&N, parms);
}

/*----- Dynamics section */

void derivs_3cl (int *neq, double *pdTime, double *y, double *ydot, double *yout, int *ip)
{
  yout[ID_Cgut] = y[ID_Agut] / Vgut ;

  yout[ID_Cliver] = y[ID_Aliver] / Vliver ;

  yout[ID_Cven] = y[ID_Aven] / Vven ;

  yout[ID_Clung] = y[ID_Alung] / Vlung ;

  yout[ID_Cart] = y[ID_Aart] / Vart ;

  yout[ID_Crest] = y[ID_Arest] / Vrest ;

  yout[ID_Ckidney] = y[ID_Akidney] / Vkidney ;

  yout[ID_Cserum] = y[ID_Aven] / Vven / Ratioblood2plasma ; 
  //(ratioblood2serum the same as for plasma)

  yout[ID_Aserum] = y[ID_Aven] / Ratioblood2plasma * ( 1 - hematocrit ) ;

  ydot[ID_Agutlumen] = - kgutabs * y[ID_Agutlumen] ;

  /* Agut flow: modified by subtracting term (lastterm) proportional to gut-specific metabolism */
  if(Km_gut==0) {
      lastterm = CLmetabolism_gut * y[ID_Agut] / Vgut / Kgut2plasma ;
  } else {
      lastterm = (Vmax_gut * y[ID_Agut] / Vgut / Kgut2plasma) / (Km_gut + y[ID_Agut] / Vgut / Kgut2plasma);
  }
  ydot[ID_Agut] = kgutabs * y[ID_Agutlumen] + Qgut * ( y[ID_Aart] / Vart - y[ID_Agut] / Vgut * Ratioblood2plasma / Kgut2plasma / Fraction_unbound_plasma ) - lastterm;
  /* End Agut flow. */
  
  
  /* Aliver flow: modified by subtracting term (lastterm) proportional to liver-specific metabolism */
  if(Km_liver==0) {
      lastterm = CLmetabolism * y[ID_Aliver] / Vliver / Kliver2plasma ;
  } else {
      lastterm = (Vmax_liver * y[ID_Aliver] / Vliver / Kliver2plasma) / (Km_liver + y[ID_Aliver] / Vliver / Kliver2plasma);
  }
  ydot[ID_Aliver] = Qliver * y[ID_Aart] / Vart + Qgut * y[ID_Agut] / Vgut * Ratioblood2plasma / Kgut2plasma / Fraction_unbound_plasma - 
    ( Qliver + Qgut ) * y[ID_Aliver] / Vliver / Kliver2plasma / Fraction_unbound_plasma * Ratioblood2plasma - lastterm;
  /* End Aliver flow. */
  
  ydot[ID_Aven] = ( ( Qliver + Qgut ) * y[ID_Aliver] / Vliver / Kliver2plasma + Qkidney * y[ID_Akidney] / Vkidney / Kkidney2plasma + 
  Qrest * y[ID_Arest] / Vrest / Krest2plasma ) * Ratioblood2plasma / Fraction_unbound_plasma - Qcardiac * y[ID_Aven] / Vven ;

  ydot[ID_Alung] = Qcardiac * ( y[ID_Aven] / Vven - y[ID_Alung] / Vlung * Ratioblood2plasma / Klung2plasma / Fraction_unbound_plasma ) ;

  ydot[ID_Aart] = Qcardiac * ( y[ID_Alung] / Vlung * Ratioblood2plasma / Klung2plasma / Fraction_unbound_plasma - y[ID_Aart] / Vart ) ;

  ydot[ID_Arest] = Qrest * ( y[ID_Aart] / Vart - y[ID_Arest] / Vrest * Ratioblood2plasma / Krest2plasma / Fraction_unbound_plasma ) ;

  /* Akidney flow: modified by subtracting term (lastterm) proportional to liver-specific metabolism */
  if(Km_kidney==0) {
    lastterm = - CLmetabolism_kidney * y[ID_Akidney] / Vkidney / Kkidney2plasma;
    ydot[ID_Akidney] = Qkidney * y[ID_Aart] / Vart - Qkidney * y[ID_Akidney] / Vkidney / Kkidney2plasma * Ratioblood2plasma / Fraction_unbound_plasma - 
      Qgfr * y[ID_Akidney] / Vkidney / Kkidney2plasma - lastterm;
  } else {
    //here not just subtracting an extra term but instead of 'Qgfr term' in the above (and sans extra subtraction):
    ydot[ID_Akidney] = Qkidney * y[ID_Aart] / Vart - Qkidney * y[ID_Akidney] / Vkidney / Kkidney2plasma * Ratioblood2plasma / Fraction_unbound_plasma - 
      ((Vmax_kidney * y[ID_Akidney])/(Km_kidney + y[ID_Akidney]));
  }
  /* End Akidney flow. */
  

  ydot[ID_Atubules] = Qgfr * y[ID_Akidney] / Vkidney / Kkidney2plasma ;
  
  /* Modification: CL partitioned into three terms */
  if(Km_kidney!=0 && Km_gut != 0 && Km_liver != 0) {
    ydot[ID_Ametabolized] = 
      (Vmax_liver * y[ID_Aliver] / (Km_liver + y[ID_Aliver])) + 
      (Vmax_kidney * y[ID_Akidney] / (Km_kidney + y[ID_Akidney])) + 
      (Vmax_gut * y[ID_Agut] / (Km_gut + y[ID_Agut])); 
      
  } else {
    ydot[ID_Ametabolized] = 
      (CLmetabolism * y[ID_Aliver] / Vliver / Kliver2plasma) + 
      (CLmetabolism_kidney * y[ID_Akidney] / Vkidney / Kkidney2plasma) + 
      (CLmetabolism_gut * y[ID_Agut] / Vgut / Kgut2plasma);
  }
  ydot[ID_AUC] = y[ID_Aven] / Vven / Ratioblood2plasma ; 
  //(not very helpful here or used to calculate later?)

} /* derivs */
