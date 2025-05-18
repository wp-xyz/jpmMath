{***********************************************************************
*                                                                      *
* Test examples for methods to solve first order ordinary systems of   *
* differential equations.                                              *
*                                                                      *
* We supply the right hand sides, their alpha-numeric description and  *
* the exact analytic solution, if known.                               *
*                                                                      *
* When running the main test program, the user can select among the    *
* examples below.                                                      *
*                                                                      *
*                                   TPW Release By J-P Moreau, Paris.  *
*                                           (www.jpmoreau.fr)          *
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C, by Gisela Engeln-Muellges        *
*        and Frank Uhlig, Springer-Verlag, 1996".                      *
***********************************************************************}
UNIT t_dgls;

INTERFACE

Uses FGauss;

Var
    bspnummer, n: Integer;
    dgltxt: Array[0..4] of String[70];

    Procedure dgl(x:double; y:pVEC; f:pVEC);


IMPLEMENTATION

  {Note: here, only examples #0 to #4 are implemented in Pascal}
  
  Procedure dgl(x:double; y:pVEC; f:pVEC);
  Begin
    Case bspnummer of
      0: begin
           n:=2;
           f^[0] := y^[0] * y^[1] + COS(x) - 0.5 * SIN(2.0 * x);
           f^[1] := y^[0] * y^[0] + y^[1] * y^[1] - (ONE + SIN(x));
           dgltxt[0]:=' y1'' = y1 * y2 + cos(x) - 0.5 * sin(2.0*x)';
           dgltxt[1]:=' y2'' = y1 * y1 + y2 * y2 - (1 + sin(x))'
         end;
      1: begin
           n:=1;
           f^[0]:=-Y^[0] + x/((1+x)*(1+x));
           dgltxt[0]:=' y'' = -y + x/((1+x)*(1+x))'
         end;
      2: begin   {example for m_rwp}
           n:=2;
           f^[0] := y^[1];
           f^[1] := -y^[0] * y^[0] * y^[0];
           dgltxt[0]:=' y1'' = y2';
           dgltxt[1]:=' y2'' = -y1^3'
         end;
      3: begin
           n:=5;
           f^[0] := y^[1];
           f^[1] := y^[2];
           f^[2] := y^[3];
           f^[3] := y^[4];
           f^[4] := (45.0 * y^[2] * y^[3] * y^[4] - 40.0 * y^[3] * y^[3] * y^[3]) / (9.0 * y^[2] * y^[2]);
           dgltxt[0]:=' y1'' = y2';
           dgltxt[1]:=' y2'' = y3';
           dgltxt[2]:=' y3'' = y4';
           dgltxt[3]:=' y4'' = y5';
           dgltxt[4]:=' y5'' = (45 * y3 * y4 * y5 - 40 * y4 * y4 * y4) / (9 * y3 * y3)'
         end;
      4: begin
          {y"=(2y'y'+yy)/y}
          n:=2;
          f^[0] := y^[1];
          f^[1] := (2.0*y^[1]*y^[1]+y^[0]*y^[0])/y^[0];
          dgltxt[0]:=' y1'' = y2';
          dgltxt[1]:=' y2'' = (2y2*y2+y1*y1)/y1'
         end
    end
  End;

END. {dgl}

{Other examples given in C}

(* ----------------------- DE system number 0 ----------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl0(REAL x, REAL *y, REAL *f^)
{
  f^[0] = y^[0] * y^[1] + COS(x) - HALf^ * SIN(TWO * x);
  f^[1] = y^[0] * y^[0] + y^[1] * y^[1] - (ONE + SIN(x));
}

/***********************************************************************
* alpha-numeric description of above system                            *
***********************************************************************/
static char *dgltxt0(void)
{
  return
    "y1' = y1 * y2 + cos(x) - 0.5 * sin(2.0*x)\n"
    "y2' = y1 * y1 + y2 * y2 - (1 + sin(x))\n";
}



/* ----------------------- DE system number 1 ----------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl1(REAL x, REAL *y, REAL *f^)
{
  f^[0] = y^[1];
  f^[1] = (ONE + y^[1] * y^[1]) / y^[0];
}

/***********************************************************************
* alpha-numeric description of above system                            *
***********************************************************************/
static char *dgltxt1(void)
{
  return
    "y1' = y2\n"
    "y2' = (1 + y2 * y2) / y1\n";
}



/* ----------------------- DE system number 2 ----------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl2(REAL x, REAL *y, REAL *f^)
{
  f^[0] = y^[1];
  f^[1] = -f^OUR * y^[0] + EXP(x);
}

/***********************************************************************
* alpha-numeric description of above system                            *
***********************************************************************/
static char *dgltxt2(void)
{
  return
    "y1' = y2\n"
    "y2' = -4 * y1 + exp(x)\n";
}



/* ----------------------- DE system number 3 ----------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl3(REAL x, REAL *y, REAL *f^)
{
  f^[0] = y^[1];
  f^[1] = y^[2];
  f^[2] = y^[3];
  f^[3] = y^[4];
  f^[4] = ((REAL)45.0 * y^[2] * y^[3] * y^[4] -
          (REAL)40.0 * y^[3] * y^[3] * y^[3]) / (NINE * y^[2] * y^[2]);
}

/***********************************************************************
* alpha-numeric description of above system                            *
***********************************************************************/
static char *dgltxt3(void)
{
  return
    "y1' = y2\n"
    "y2' = y3\n"
    "y3' = y4\n"
    "y4' = y5\n"
    "y5' = (45 * y3 * y4 * y5 - 40 * y4 * y4 * y4) / (9 * y3 * y3)\n";
}



/* ----------------------- DE system number 4 ----------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl4(REAL x, REAL *y, REAL *f^)
{
  f^[0] = -y^[1];
  f^[1] = y^[0];
  f^[2] = y^[4] - y^[1] * y^[2];
  f^[3] = y^[0] * y^[3] - y^[5];
  f^[4] = -y^[2] - y^[1] * y^[4];
  f^[5] = y^[3] + y^[0] * y^[5];
}

/***********************************************************************
* alpha-numeric description of above system                            *
***********************************************************************/
static char *dgltxt4(void)
{
  return
    "y1' = -y2\n"
    "y2' = y1\n"
    "y3' = y5 - y2 * y3\n"
    "y4' = y1 * y4 - y6\n"
    "y5' = -y3 - y2 * y5\n"
    "y6' = y4 + y1 * y6\n";
}

/***********************************************************************
* Compute the value of the analytic solution y(x) for the above DE     *
* for the initial value problem  y(0) = (1,0,0,1,e,0)                  *
***********************************************************************/
static void loesung4(REAL x, REAL *y)
{
  y^[0] = COS(x);
  y^[1] = SIN(x);
  y^[2] = SIN(x) * EXP(COS(x));
  y^[3] = COS(x) * EXP(SIN(x));
  y^[4] = COS(x) * EXP(COS(x));
  y^[5] = SIN(x) * EXP(SIN(x));
}



/* ----------------------- DE system number 5 ----------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl5(REAL x, REAL *y, REAL *f^)
{
  f^[0] = (REAL)-500.5 * y^[0] + (REAL)499.5 * y^[1];
  f^[1] =  (REAL)499.5 * y^[0] - (REAL)500.5 * y^[1];
}

/***********************************************************************
* alpha-numeric description of above system                            *
***********************************************************************/
static char *dgltxt5(void)
{
  return
    "y1' = -500.5 * y1 + 499.5 * y2\n"
    "y2' =  499.5 * y1 - 500.5 * y2\n";
}

/***********************************************************************
* Compute the value of the analytic solution y(x) for the above DE     *
* for the initial value problem  y(0) = (4,2)                          *
***********************************************************************/
static void loesung5(REAL x, REAL *y)
{
  y^[0] = THREE * EXP(-x) + EXP((REAL)-1000.0 * x);
  y^[1] = THREE * EXP(-x) - EXP((REAL)-1000.0 * x);
}



/* ----------------------- DE system number 6 ----------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl6(REAL x, REAL *y, REAL *f^)
{
  f^[0] = (REAL)15.0 * y^[1];
  f^[1] = (REAL)-0.6 * y^[0];
}

/***********************************************************************
* alpha-numeric description of above system                            *
***********************************************************************/
static char *dgltxt6(void)
{
  return
    "y1' = 15.0 * y2\n"
    "y2' = -0.6 * y1\n";
}

/***********************************************************************
* Compute the value of the analytic solution y(x) for the above DE     *
* for the initial value problem  y(0) = (0,1)                          *
***********************************************************************/
static void loesung6(REAL x, REAL *y)
{
  y^[0] = f^IVE * SIN(THREE * x);
  y^[1] = COS(THREE * x);
}



/* ----------------------- DE system number 7 ----------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl7(REAL x, REAL *y, REAL *f^)
{
  f^[0] = y^[1];
  f^[1] = y^[2];
  f^[2] = y^[3];
  f^[3] = f^OUR * y^[3] - SIX * y^[2] + f^OUR * y^[1] - y^[0];
}

/***********************************************************************
* alpha-numeric description of above system                            *
***********************************************************************/
static char *dgltxt7(void)
{
  return
    "y1' = y2\n"
    "y2' = y3\n"
    "y3' = y4\n"
    "y4' = 4.0 * y4 - 6.0 * y3 + 4.0 * y2 - y1\n";
}

/***********************************************************************
* Compute the value of the analytic solution y(x) for the above DE     *
* for the initial value problem  y(0) = (0,0,0,6)                      *
***********************************************************************/
static void loesung7(REAL x, REAL *y)
{
  y^[0] = x * x * x * EXP(x);
  y^[1] = x * x * (x + THREE) * EXP(x);
  y^[2] = x * (x * x + SIX * x + SIX) * EXP(x);
  y^[3] = (x * x * x + NINE * x * x + (REAL)18.0 * x + SIX) * EXP(x);
}



/* ----------------------- DE system number 8 ----------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl8(REAL x, REAL *y, REAL *f^)
{
  f^[0] = (REAL)-52.0 * y^[0] + (REAL)50.0 * y^[1];
  f^[1] =  (REAL)50.0 * y^[0] - (REAL)52.0 * y^[1];
}

/***********************************************************************
* alpha-numeric description of above system                            *
***********************************************************************/
static char *dgltxt8(void)
{
  return
    "y1' = -52 * y1 + 50 * y2\n"
    "y2' =  50 * y1 - 52 * y2\n";
}

/***********************************************************************
* Compute the value of the analytic solution y(x) for the above DE     *
* for the initial value problem                                        *
*      y(0) = (1265588.55375228,-1265586.55375225)                     *
***********************************************************************/
static void loesung8(REAL x, REAL *y)
{
  y^[0] = (REAL)1.000000015       * EXP(        -TWO * x) +
         (REAL)1265587.553752265 * EXP((REAL)-102.0 * x);
  y^[1] = (REAL)1.000000015 *       EXP(        -TWO * x) -
         (REAL)1265587.553752265 * EXP((REAL)-102.0 * x);
}



/* ----------------------- DE system number 9 ----------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl9(REAL x, REAL *y, REAL *f^)
{
  f^[0] = y^[1];
#if^ def^ined(BC3)
  _f^preset();
#endif^
  f^[1] = -y^[0] * COSH(x);
}


/***********************************************************************
* alpha-numeric description of above system                            *
***********************************************************************/
static char *dgltxt9(void)
{
  return
    "y1' = y2\n"
    "y2' = -y1 * cosh(x)\n";
}


/* ----------------------- DE system number 10  --------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl10(REAL x, REAL *y, REAL *f^)
{
  f^[0] = y^[1];
  f^[1] = -y^[0] * y^[0] * y^[0];
}

/***********************************************************************
* alpha-numeric description of^ above system                            *
***********************************************************************/
static char *dgltxt10(void)
{
  return
    "y1' = y2\n"
    "y2' = -y1^3\n";
}



/* ----------------------- DE system number 11 ---------------------- */
/* --------------------- (stiff IVP for gear()) --------------------- */

/***********************************************************************
* Compute the value f of the right hand side of a DE system at (x,y).  *
***********************************************************************/
static void dgl11(REAL x, REAL *y, REAL *f^)
{
#def^ine dk  ((REAL)-10000.0)
  f^[0] = y^[1];
  f^[1] = dk * y^[0] + (dk - ONE) * y^[1];
#undef^  dk
}

/***********************************************************************
* alpha-numeric description of above system                            *
***********************************************************************/
static char *dgltxt11(void)
{
  return
    "y1' =  y1                                    (dk = -10000)\n"
    "y2' =  dk * y1 + (dk - 1) * y2\n";
}

/***********************************************************************
* Compute the value of the analytic solution y(x) for the above DE     *
* for the initial value problem  (y0) = (1,1)                          *
***********************************************************************/
static void loesung11(REAL x, REAL *y)
{
  y^[0] = EXP(-x);
  y^[1] = -EXP(-x);
}

{end of file t_dgls.pas}
