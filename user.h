/* Physical parameters */
#ifdef USER
int iqd, iqu0, iqe, iqk, iqeps, iqal0, iqie, iqpot;
float g, pfloor, dfloor, avisp, avise, tcond, sdiff, vis, ce1,
 ce2, cs, ct, cd, ce, ceps, kmin, epsmin, mtl, gconst, rat, heat, g1, 
  g2, g3, g4, g5, g6, g7, g8;
bool mon, fxp, inten;
FILE *monfile;
#else
extern int iqd, iqu0, iqe, iqk, iqeps, iqal0, iqie, iqpot;
extern float g, pfloor, dfloor, tfloor, avisp, avise, tcond, sdiff, vis,
 ce1, ce2, cs, ct, cd, ce, ceps, kmin, epsmin, mtl, gconst, rat, heat,
  g1, g2, g3, g4, g5, g6, g7, g8;
extern bool mon, fxp, inten;
extern FILE *monfile;
#endif
