{***********************************************************************
* This program computes the value of the integral of func(x) over the  *
* interval (a,b) by using the summed Clenshaw-Curtis formula           *
* -------------------------------------------------------------------- *
* SAMPLE RUN:                                                          *
* (Integrate function func(x) = cos(x) - x*sin(x) from x=0 to x=PI)    *
*                                                                      *
* Here parameters m=4 (number of subintervals in [0, PI]),             *
*                 n=6 (number of nodes for Clenshaw-Curtis quadrature) *
*                                                                      *
* Integral = -3.14159265357197E+0000                                   *
* Error code = 0                                                       *
*                                                                      *
* (Exact value is -PI).                                                *
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C,  By Gisela Engeln-Muellges       *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
*                                                                      *
*                                  TPW Release By J-P Moreau, Paris.   *
*                                         (www.jpmoreau.fr)            *
***********************************************************************}
Program Test_Clencurt;

Uses WinCrt;

Const
      ONE  = 1.0;
      ZERO = 0.0;
      SIZE = 10;

Type
      VECT = Array[0..SIZE] of Double;

Var
      error, m, n: Integer;
      t, Tk, Ak  : VECT;
      Result     : Double;

{Function to integrate}
Function func(x:Double): Double;
Begin
  func := cos(x)-x*sin(x)
End;

Function ClenCurtStGew(n:integer; Var StStelle, Gewicht:VECT): Integer;
{***********************************************************************
* Computes the nodes and weights of a Clenshaw-Curtis quadrature       *
* formula of local error order  n+3 for the reference interval [-1,1]. *
* -------------------------------------------------------------------- *
* Parameters:                                                          *
*   integer n            n+1 is the number of nodes and weights:       *
*                        n > 1, n odd                                  *
*   double StStelle[]    computed nodes                                *
*   double Gewicht[]     computed weights                              *
*                                                                      *
* Return value :                                                       *
*   0:                   all ok                                        *
*   1:                   n too small or odd                            *
*                                                                      *
* Author:                Uli Eggermann,  9.24.1990                     *
***********************************************************************}
Var
  k, j, m: integer;
  p, f, g, h, i, d: Double;
Begin
  i:=ONE;
  if (n < 2) or ((n mod 2) <> 0) then             { n - Test          }
  begin
    ClenCurtStGew := 1;
    exit
  end;
  m := n Div 2;
  d := ONE*(n * n - 1);
  h := 2.0 / (n * d);
  f := 4.0 / (ONE*n);
  p := PI/(ONE*n);

  Gewicht[0] := ONE/d; Gewicht[n] := ONE/d;         {  left end point }
  StStelle[0] := -ONE; StStelle[n] := ONE;          { right end point }

  for k := 1 to m do                                { Interval [-1,0] }
  begin
    g := ZERO;
    for j := 1 to m-1 do
      g := g + cos(2.0 * j * k * p) / (4.0 * j * j - 1);
    Gewicht[k] := h * (d + i) - f * g;
    StStelle[k] := - cos(k * p);
    i := -i
  end;

  StStelle[m] := ZERO;                          { center of interval }

  for k := m+1 to n-1 do                           { Interval [0, 1] }
  begin
    Gewicht[k] :=  Gewicht[n - k];
    StStelle[k] := -StStelle[n - k]
  end;
  ClenCurtStGew := 0
End; {ClenCurtStGew}

Function ClenCurt (
              m: integer;
              t: VECT;
              n: integer;
              Tk, Ak: VECT;
          Var Resultat: double
              ): Integer;
{***********************************************************************
* Computes the value of the integral of func (x) over the interval     *
* (a,b) with the partition                                             *
*              t: a = t[0] < t[1] < .. < t[m] = b                      *
* by using the summed  Clenshaw-Curtis formula.                        *
* This program uses precomputed [0..n] weight vectors and Chebyshev    *
* nodes.                                                               *
*                                                                      *
* Parameter:                                                           *
*   int    m              number of subintervals                       *
*   REAL   t []           partition                                    *
*   int    n              n + 1 = number of nodes, n > 1, n even       *
*                                   (n + 2 = global error order)       *
*   REAL   Tk []          Chebyshev nodes                              *
*   REAL   Ak []          Weights                                      *
*                           Tk and Ak must be made available before    *
*                           calling this function by the procedure     *
*                           ClenCurtStGew for example                  *
*   REAL   *Resultat      Compute integral value                       *
*                                                                      *
* Return value :                                                       *
*   0:                    o.k.                                         *
*   1:                    improper number of nodes                     *
*   2:                    improper number of sub intervals             *
*                                                                      *
* Author:                 Uli Eggermann, 10.3.1991                     *
***********************************************************************}
Var
  j, k: integer;
  v, h, sum: Double;
Begin
  if (n < 2) or ((n Mod 2) <> 0) then        { n positive ? n even  ? }
  begin
    ClenCurt := 1;
    exit
  end;
  if m < 1 then                              { partition              }
  begin
    ClenCurt := 2;
    exit
  end;

  Resultat := ZERO;
  for j := 0 to m-1 do                       { loop over intervals    }
  begin
    v := 0.5 * (t[j+1] - t[j]);              { half the interval size }
    h := 0.5 * (t[j+1] + t[j]);              { Interval center        }

    sum := ZERO;
    for k := 0 to n do                       { Chebyshev loop         }
      sum := sum + Ak[k]*func(v*Tk[k] + h);

    Resultat := Resultat + v*sum;
  end;
  ClenCurt :=  0
End;

{main program}
BEGIN

  m:=4; n:=6;

  t[0]:=ZERO; t[1]:=PI/4.0; t[2]:=PI/2.0; t[3]:=3.0*PI/4.0; t[4]:=PI;

  ClenCurtStGew(n,Tk,Ak);

  error := ClenCurt(m,t,n,Tk,Ak,Result);

  writeln;
  writeln(' Integral = ', Result);
  writeln(' Error code = ', error);

  ReadKey;
  DoneWinCrt

END.

{ ------------------------- END clencurt.pas ------------------------ }