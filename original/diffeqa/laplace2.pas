{*************************************************************************
* Solve Laplace Equation by relaxation Method: d2T/dx2 + d2T/dy2 = 0 (3) *
* ---------------------------------------------------------------------- *
* Example #3: Determine the temperatures T(x,y) in a domain [x,y] made   *
*             of one square of length C with the following conditions:   *
*             along OC: T=0, along AB: T=1, along OA and CB: dT/dy = 0   *
*             (see sketch below).                                        *
*             Main parameters:                                           *
*             C: length of square domain.                                *
*             Number of subdivisions 8 in ox and 8 in Oy.                *
*             integration steps dx=C/8.                                  *
*             weight coefficient w = 1.3                                 *
*             Maximum residual error ER = 0.001                          *
*             Maximum number of iterations IT = 100                      *
* ---------------------------------------------------------------------- *
* SAMPLE RUN:                                                            *
*                                                                        *
* Niter  Max. Residual     W                                             *
*  43      8.6E-0004      1.30                                           *
* Temperature                                                            *
*  0.000  0.125  0.249  0.374  0.499  0.624  0.749  0.875  1.000         *
*  0.000  0.125  0.249  0.374  0.499  0.624  0.749  0.875  1.000         *
*  0.000  0.124  0.249  0.374  0.499  0.624  0.749  0.875  1.000         *
*  0.000  0.124  0.249  0.374  0.499  0.624  0.749  0.875  1.000         *
*  0.000  0.124  0.249  0.374  0.499  0.624  0.749  0.875  1.000         *
*  0.000  0.124  0.249  0.374  0.499  0.624  0.749  0.875  1.000         *
*  0.000  0.124  0.249  0.374  0.499  0.624  0.749  0.875  1.000         *
*  0.000  0.124  0.249  0.373  0.498  0.624  0.749  0.875  1.000         *
*  0.000  0.124  0.249  0.373  0.498  0.624  0.749  0.874  1.000         *
*                                                                        *
* ---------------------------------------------------------------------- *
* REFERENCE:  "Méthode de calcul numérique- Tome 2 - Programmes en Basic *
*              et en Pascal By Claude Nowakowski, Edition du P.S.I.,     *
*              1984" [BIBLI 04].                                         *
*                                                                        *
*                                   TPW Release By J-P Moreau, Paris.    *
*                                           (www.jpmoreau.fr)            *
**************************************************************************
{ Note: Here, what is new is the limit conditions involving partial derivatives dT/dy.
        Inside domain, the classical scheme is applied to calculate the residual value at current
        point (i,j): 4*T(i,j)-T(i-1,j)-T(i+1,j)-Ti,j-1)-T(i,j+1). Must be near zero (< ER).
        For limit lines OA and CB, the partial derivative is approximated by:
                        dT/dy = (T(i,j+1)-T(i,j-1))/(2*dy)
        But the point T(i,j+1) for segment CB is out of the domain and becomes a fictitious point,
        so we apply the following scheme:  r=4*T(i,j)-2*T(i,j-1)-T(i-1,j)-T(i+1,j).
        For the other segment OA, it becomes: r=4*T(i,j)-2*T(i,j+1)-T(i-1,j)-T(i+1,j).
        The exact solution for this problem is T=K*x (0<x<C), here k=0.25.
        This allows verifying the results accuracy.
--------------------------------------------------------------------------------------------------------}          
Program laplace2;

Uses WinCrt;

Const                 {           dT/dy=0           }
      IT = 100;       {       C-------------B       }
      ER = 0.001;     {       |             |       }
      OM = 1.3;       {       |             |       }
      C  = 4;         {   T=0 |             |T=1    }
      NM = 8;         {       |             |       }
                      {       |             |       }
                      {       O-------------A       }
                      {           dT/dy=0           }

Label 10, 20;
Var
      T:Array[0..NM,0..NM] of Real;
      i,j,l,n1: Integer;
      r,rm: Real;


      {calculate residual for inside square domain}
      Procedure Schema;
      Begin
        r:=T[i+1,j]+T[i-1,j]+T[i,j+1]+T[i,j-1]-4*T[i,j];
        T[i,j]:=T[i,j]+0.25*OM*r;
        If Abs(r) > rm then rm:=Abs(r)
      End;

      {calculate residual for lower limit}
      Procedure Schema1;
      Begin
        r:=T[i+1,j]+T[i-1,j]+2*T[i,j+1]-4*T[i,j];
        T[i,j]:=T[i,j]+0.25*OM*r;
        If Abs(r) > rm then rm:=Abs(r)
      End;

      {calculate residual for upper limit}
      Procedure Schema2;
      Begin
        r:=T[i+1,j]+T[i-1,j]+2*T[i,j-1]-4*T[i,j];
        T[i,j]:=T[i,j]+0.25*OM*r;
        If Abs(r) > rm then rm:=Abs(r)
      End;


BEGIN

  {Fix limit conditions}
  For i:=0 to NM do
    For j:=0 to NM do
      T[i,j]:=0.0;
  For j:=0 to NM do
    T[NM,j]:=1.0;

  n1:=NM-1;
  {main calculation loop}
  For l:=1 to IT do
  begin
    rm:=0.0;
    j:=0;
    For i:=1 to n1 do Schema1;
    For j:=1 to n1 do
      For i:=1 to n1 do
        Schema;
    j:=NM;
    For i:=1 to n1 do Schema2;
    if rm < ER then goto 10;
  end;
  writeln(' No convergence after ',IT,' iterations.');
  writeln(' Residual rm = ', rm);
  goto 20;
10:writeln;
  writeln(' Niter  Max. Residual     W');
  writeln(l:4, '     ', rm:-1, '    ',OM:6:2);
  writeln(CHR(29),' Temperature');
  for j:=NM downto 0 do
  begin
    for i:=0 to NM do write(T[i,j]:7:3);
    writeln
  end;

20:Readkey;
  DoneWinCrt 

END.

{end of file laplace2.pas}