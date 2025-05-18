{ ------------------------ UNIT kubnec.pas -------------------------- *
*     "Numerical Algorithms with C, By Gisela Engeln-Muellges         *
*      and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].            *
*                                                                     *
*                                  TPW Release By J-P Moreau, Paris.  *
*                                         (www.jpmoreau.fr)           *
----------------------------------------------------------------------}
UNIT KUBNEC;

INTERFACE

  Uses WinCrt;

  Const NMAX = 25;

  Type VEC = array[0..NMAX] of Double;


  Function Kub4NeCn(
              a, b: Double; Nx: Integer;
              c, d: Double; Ny: Integer;
              Verfahren: Integer;
          Var Wert: Double;
              Schaetzen: Integer;
          Var FehlerSch: Double
             ):Integer;

  Function Kub3RoRi ( Px, Py, Qx, Qy, Rx, Ry: Double;
                      n: integer;
                      Var Wert, FehlerSch: Double
                    ): Integer;


IMPLEMENTATION

{**********************************************************************
* Global  constants and variables                                     *
**********************************************************************}
Var
  K_1: Array[0..1] of Double;
  K_2: Array[0..2] of Double;
  K_3: Array[0..3] of Double;
  K_4: Array[0..4] of Double;
  K_5: Array[0..5] of Double;
  K_6: Array[0..6] of Double;
  K_7: Array[0..7] of Double;

  KubVer    : Integer;           {global variables:  method number}
  KubX, KubY: Integer;                  {number of x-, y-intervals}

  {user defined function(x,y) }
  Function func(x,y:Double):Double;
  begin
    func := x*y
  end;


  Function K4KnotGew(i, j: Integer): Double; Forward;

  Function Kub4NeCn(
              a, b: Double; Nx: Integer;
              c, d: Double; Ny: Integer;
              Verfahren: Integer;
          Var Wert: Double;
              Schaetzen: Integer;
          Var FehlerSch: Double
             ):Integer;
{***********************************************************************
* Cubature over rectangles using Newton-Cotes formulas.                *
*                                                                      *
* Integrate the function f(x,y) over a rectangle (a,b) x (c,d) using   *
* the summed Newton-Cotes cubature formulas for sub-rectangles.        *
*                                                                      *
* Parameters:                                                          *
*   Double a,b        left, right x-end points                         *
*   int    Nx         number of x-intervals                            *
*   Double c,d        ditto for y                                      *
*   int    Ny                                                          *
*   int    Verfahren  Number of method :                               *
*                       1: trapezoidal rule                            *
*                       2: Simpson's rule                              *
*                       3: 3/8                                         *
*                       4: 4/90                                        *
*                       5: 5/288                                       *
*                       6: 6/840                                       *
*                       7: 7/17280 formula                             *
*   double  f (x,y)   function                                         *
*   Double  Wert      value for integral                               *
*   int     Schaetzen 0: no estimation                                 *
*                     1: estimate error by repeating cubature for half *
*                        the step size                                 *
*   Double  FehlerSch error estimate for Wert                          *
*                                                                      *
* Return value :                                                       *
*   0:              o.k.                                               *
*   1:              Nx improper                                        *
*   2:              Ny improper                                        *
*   3:              method number incorrect                            *
*   4:              Integration interval of length zero                *
*                                                                      *
* Author:           Uli Eggermann, 3.31.1996                           *
************************************************************************

{**********************************************************************
* Cubature over rectangles via  Newton-Cotes                          *
**********************************************************************}
Label Return;
Var
  i,j,k,kend,Ordnung: Integer;
  Hab, Hcd, Wert1: Double;
Begin

  Wert1 := 0.0;
  if Nx < 1 then
  begin
    Kub4NeCn := 1;
    goto Return
  end;
  if Ny < 1 then
  begin
    Kub4NeCn := 2;
    goto Return
  end;
  if (Verfahren < 1) or (Verfahren > 7) then
  begin
    Kub4NeCn := 3;
    goto Return
  end;
  if (a = b) or (c = d) then
  begin
    Kub4NeCn := 4;
    goto Return
  end;
                                           {method:     1,2,3,4,5,6,7}
  Ordnung := (Verfahren Div 2) * 2 + 2;    {order:      2,4,4,6,6,8,8}

  KubVer := Verfahren;                       {copy to global variable}

  if Schaetzen <> 0 then kend := 2 else kend := 1;
  for k := 1 to kend do
  begin
    KubX := k * Nx * KubVer;
    KubY := k * Ny * KubVer;
    Wert := 0.0;                                        { initialize }
    Hab  := (b - a) / KubX;                             { step sizes }
    Hcd  := (d - c) / KubY;

    for i := 0 to KubX do                                 { Cubature }
      for j := 0 to KubY do
        Wert := Wert + K4KnotGew(i,j) * func(a + i * Hab, c + j * Hcd);
      
    Wert := Wert * Hab * Hcd * KubVer * KubVer;      { final multipl.}
    if (Schaetzen<>0) and (k = 1) then Wert1 := Wert   { store value }
  end;

  if Schaetzen <> 0 then                              {estimate}
    FehlerSch := (Wert1 - Wert) / ((1 SHL Ordnung) - 1);

  Kub4NeCn := 0;

Return: End;


Function K4KnotGew(i, j: Integer): Double;
{***********************************************************************
* Local function for finding node weights                              *
*                                                                      *
* Weight the functional values at the nodes depending on their location*
* (edge, interior, ...) in the summed cubature formulas.               *
***********************************************************************}
Var
  f: Double;
  k: integer;

  Function Faktor(a,b:Integer):Double;
  Begin
    if (k=0) and (a>0) and (a<b) then
      Faktor := 2.0
    else
      Faktor := 1.0
  End;

Begin
  {***********************************************************
  * for x-direction :                                        *
  *   1) node at interval end  (k == 0)                      *
  *   2) node not at left end (a > 0) nor at right end (a<b) *
  * for y-direction :                                        *
  *   1) node at interval end  (k == 0)                      *
  *   2) node not at bottom (a > 0) and not at top (a < b)   *
  ***********************************************************}
  k := i Mod KubVer;
  Case KubVer of
   1: f := K_1[k] * Faktor (i, KubX);
   2: f := K_2[k] * Faktor (i, KubX);
   3: f := K_3[k] * Faktor (i, KubX);
   4: f := K_4[k] * Faktor (i, KubX);
   5: f := K_5[k] * Faktor (i, KubX);
   6: f := K_6[k] * Faktor (i, KubX);
   7: f := K_7[k] * Faktor (i, KubX)
  end;
  k := j Mod KubVer;
  Case KubVer of
   1: f := f * K_1[k] * Faktor (j, KubY);
   2: f := f * K_2[k] * Faktor (j, KubY);
   3: f := f * K_3[k] * Faktor (j, KubY);
   4: f := f * K_4[k] * Faktor (j, KubY);
   5: f := f * K_5[k] * Faktor (j, KubY);
   6: f := f * K_6[k] * Faktor (j, KubY);
   7: f := f * K_7[k] * Faktor (j, KubY)
  end;

  K4KnotGew := f

End;


Function Kub3NeC3 (
                    Px, Py, Qx, Qy, Rx, Ry: Double;
                    n: integer;
                    Var Wert: Double
                  ): Integer;
{***********************************************************************
* Cubature over triangles using the  Newton-Cotes formulas             *
*                                                                      *
* This function integrates f (x,y) over the triangle PQR using the     *
* summed 3 point cubature formulas of Newton-Cotes on sub-triangles.   *
* The sub-triangles are obtained by a regular partition of the edges   *
* of PQR into n equal parts.                                           *
*                                                                      *
* Input parameters:                                                    *
*   Double   Px,Py   coordinates of P                                  *
*   Double   Qx,Ry   coordinates of Q                                  *
*   Double   Rx,Ry   coordinates of R                                  *
*   integer  n       partion number along edges                        *
*   Double   Wert    value of integral                                 *
*                                                                      *
* Return value :                                                       *
*   0:               o.k.                                              *
*   1:               n improper                                        *
*   2:               corners  P, Q and R are collinear                 *
*                                                                      *
* Author:            Uli Eggermann, 8.1.1990 (C version)               *
***********************************************************************}
Label return;
Const
       Kub3NeC3Epsilon = 0.000001;               {for collinearity test}

{***********************************************************************
* Cubature over triangle via summed 3-point Newton-Cotes formula       *
***********************************************************************}
Var
   hPQx, hPQy, hPRx, hPRy, Area: Double;
   X, Y: Double;
   i,j, wieoft: integer;
Begin

   i:=0;
   if n < 1 then
   begin
     Kub3NeC3 := 1;
     goto return                                            {n ok ? }
   end;

   Area :=  Px * (Qy - Ry)
          + Qx * (Ry - Py)
          + Rx * (Py - Qy);                             {P, Q and R }
   if ABS(Area) < Kub3NeC3Epsilon then                  {collinear? }
   begin
     Kub3NeC3 := 2;
     goto return
   end;          

   n := n * 2;                               {number of halved edges}
   hPQx := (Qx - Px) / n;
   hPQy := (Qy - Py) / n;                       {halve the vector PQ}
   hPRx := (Rx - Px) / n;
   hPRy := (Ry - Py) / n;                       {halve the vector PR}

   Wert := 0.0;                                         {integral = 0 }

   for j := 0 to n-1 do                              {j moves along PR}
     for i := 0 to  n-j do                           {i moves along PQ}
       if ((i mod 2) <> 0) or ((j mod 2) <> 0) then     {sum if needed}
       begin
         X := Px + hPQx*i + hPRx*j;
         Y := Py + hPQy*i + hPRy*j;
         if (i=0) or (j=0) or (i=n-j) then wieoft := 1
         else wieoft := 2;
         Wert := Wert + wieoft * func(X,Y)
       end;

   X := Px + hPQx*i + hPRx*j;
   Y := Py + hPQy*i + hPRy*j;
   if (i=0) or (j=0) or (i=n-j) then wieoft := 1
   else wieoft := 2;
   Wert := Wert + wieoft * func(X,Y);                           { sum }
   Wert := Wert * Area / (1.5 * n * n);                   {last factor}

   Kub3NeC3 := 0;

return: End;


Function  Kub3Nec3n (
                      Px, Py, Qx, Qy, Rx, Ry: Double;
                      n: integer;
                      Var W: VEC
                    ): Integer;
{***********************************************************************
* computing and storing the n cubature values in given array W.        *
***********************************************************************}
Label 10, return;
Var
  i,ret: integer;
Begin
  if n < 1 then
  begin
    Kub3Nec3n := 1;
    goto return
  end;

  for i := 0 to n-1 do
  begin
    ret:=Kub3NeC3(Px,Py,Qx,Qy,Rx,Ry, 1 SHL i, W[i]);
    if ret <> 0 then goto 10;  {break}
  end;
10: Kub3Nec3n := ret;
return: end;


Function RoRiExtr (
                    RoFo: VEC;
                    nRoFo:Integer;
                    Ordnung:Integer;
                    Var Wert: Double;
                    var FehlerSch: Double
                  ): Integer;
{***********************************************************************
* Richardson extrapolation for a given Romberg sequence                *
*                                                                      *
* Parameter:                                                           *
*   Double  RoFo []   given  Romberg sequence                          *
*   integer nRoFo     length of  RoFo                                  *
*   integer Ordnung   Order of method                                  *
*   Double  Wert      final value of extrapolation                     *
*   Double  FehlerSch error estimate of Wert                           *
*                                                                      *
* Return value :                                                       *
*   0:              o.k.                                               *
*   1:              lack of available memory (not used here)           *
*   2:              nRoFo too small                                    *
*                                                                      *
* Author          Uli Eggermann, 3.31.1996                             *
***********************************************************************}
Label return;
Var
  j, k, p, s: Integer;
  RiEx: VEC;
Begin

  if nRoFo < 2 then
  begin
    RoRiExtr := 2;
    goto return
  end;

  for k := 0 to nRoFo - 1 do                                    { copy }
    RiEx[k] := RoFo[k];

  s := 1 SHL Ordnung;                                      { 2 ^ order }
  p := s;

  for k := nRoFo - 1 Downto 1 do
  begin
    p := p * s;
    for j := 0 to k-1 do 
      RiEx[j] := (p * RiEx[j+1] - RiEx[j]) / (p - 1)
  end;

  Wert      := RiEx[0];
  FehlerSch := RiEx[0] - RiEx[1];

  RoRiExtr := 0;

return: End;


Function Kub3RoRi ( Px, Py, Qx, Qy, Rx, Ry: Double;
                    n: integer;
                    Var Wert, FehlerSch: Double
                  ): Integer;
{***********************************************************************
* cubature over triangular regions using the summed 3 point formula    *
* and Romberg-Richardson extrapolation:                                *
*                                                                      *
* Function f (x,y) is integrated approximately over the triangle PQR   *
* using the summed 3 point cubature formula on similar subtriangles.   *
* The number of cubatures with halved subtriangle edge length each is  *
* defined by n.                                                        *
*                                                                      *
* Input parameters:                                                    *
*   Double   Px,Py   coordinaten of P                                  *
*   Double   Qx,Ry   coordinaten of Q                                  *
*   Double   Rx,Ry   coordinaten of R                                  *
*   integer  n       number of cubatures                               *
*                                                                      *
* Output parameters:                                                   *
*   Double   Wert    approximation of the double integral              *
*   Double   FehlSch error estimate for `Wert'                         *
*                                                                      *
* Return value:                                                        *
*   0:               o.k.                                              *
*   1:               n improper                                        *
*   2:               corners P, Q and R are collinear                  *
*   3:               no more memory for auxiliary vector               *
*                    (not used here).                                  *
*                                                                      *
* Subroutines used:                                                    *
*   RoRiExtr      Romberg-Richardson extrapolation                     *
*   Kub3NeC3      computation of cubature value                        *
*   Kub3Nec3n     computation and storing of the n cubature values     *
*                                                                      *
* author         Uli Eggermann, 07.05.1990 (C version)                 *
***********************************************************************}
Label return;
Var
  i,ii: integer;
  W: VEC;
Begin
  if n < 1 then
  begin
    Kub3RoRi:=1;
    goto return
  end;

  {storing Newton-Cotes values in the first column}
  i := Kub3Nec3n (Px,Py, Qx,Qy, Rx,Ry, n, W);

                         {computing the resulting Richardson columns, }
  if i = 0 then                  {but only if the first column is ok. }
    if n = 1 then
    begin
      Wert := W[0];
      FehlerSch := 0.0;
      i := 0
    end
    else
      i := RoRiExtr (W, n, 2, Wert, FehlerSch);

{     write(' W='); For ii:=0 to n-1 do write(' ',W[ii]); writeln;
      writeln(' n=',n,' Wert=',Wert,' Fehler=', FehlerSch);
      ReadKey; } 

  Kub3RoRi:=i;
Return: end;


BEGIN

  {Initialize constants}
  K_1[0] := 1.0/2.0;  K_1[1] := 1.0/2.0;
  K_2[0] := 1.0/6.0;  K_2[1] := 4.0/6.0;   K_2[2] := 1.0/6.0;
  K_3[0] := 1.0/8.0;  K_3[1] := 3.0/8.0;   K_3[2] := 3.0/8.0;   K_3[3] := 1.0/8.0;
  K_4[0] := 7.0/90.0; K_4[1] := 32.0/90.0; K_4[2] := 12.0/90.0; K_4[3] := 32.0/90.0; K_4[4] := 7.0/90.0;
  K_5[0] := 19.0/288.0; K_5[1] := 75.0/288.0; K_5[2] := 50.0/288.0;
  K_5[3] := 50.0/288.0; K_5[4] := 75.0/288.0; K_5[5] := 19.0/288.0;
  K_6[0] := 41.0/840.0; K_6[1] := 216.0/840.0; K_6[2] := 27.0/840.0; K_6[3] := 272.0/840.0;
  K_6[4] := 27.0/840.0; K_6[5] := 216.0/840.0; K_6[6] := 41.0/840.0;
  K_7[0] :=  751.0/17280.0; K_7[1] := 3577.0/17280.0; K_7[2] := 1323.0/17280.0;
  K_7[3] := 2989.0/17280.0; K_7[4] := 2989.0/17280.0; K_7[5] := 1323.0/17280.0;
  K_7[6] := 3577.0/17280.0; K_7[7] :=  751.0/17280.0;

  KubVer := -1

END.

{ ------------------------ END kubnec.pas --------------------------- }