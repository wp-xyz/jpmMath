{***********************************************************************
*  Solve Laplace Equation by relaxation Method: d2T/dx2 + d2T/dy2 = 0  *
* -------------------------------------------------------------------- *
* Example: Determine the temperatures T(x,y) in a square plate [x,y]   *
*          of side length C with following limit conditions:           *
*          1) T = 0 for sides x = 0 and y = 0.                         *
*          2) T = x^3 for side y = C                                   *
*          3) T = 16*x for side x = C  (see sketch below).             *
*          Other parameters:                                           *
*          square length C = 4 (any unit)                              *
*          Number of subdivisions NM = 8 (9 knots per line from 0 to 8)*
*          integration steps dx or dy = C/NM                           *
*          weight coefficient w = 1.45                                 *
*          Maximum residual error ER = 0.01                            *
*          Maximum number of iterations IT = 40                        *
* -------------------------------------------------------------------- *
* SAMPLE RUN:                                                          *
*                                                                      *
* Niter  Max. Residual     W                                           *
*  18      3.1E-0003      1.45                                         *
* Temperature                                                          *
*  0.000  0.125  1.000  3.375   8.000  15.625  27.000  42.875  64.000  *
*  0.000  1.401  3.401  6.586  11.515  18.693  28.509  41.106  56.000  *
*  0.000  2.078  4.616  8.052  12.783  19.122  27.237  37.040  48.000  *
*  0.000  2.296  4.932  8.224  12.442  17.774  24.278  31.818  40.000  *
*  0.000  2.175  4.592  7.469  10.987  15.254  20.283  25.953  32.000  *
*  0.000  1.814  3.790  6.075   8.782  11.974  15.646  19.711  24.000  *
*  0.000  1.291  2.681  4.259   6.092   8.213  10.616  13.244  16.000  *
*  0.000  0.668  1.384  2.189   3.113   4.171   5.359   6.651   8.000  *
*  0.000  0.000  0.000  0.000   0.000   0.000   0.000   0.000   0.000  *
*                                                                      *
* -------------------------------------------------------------------- *
* REFERENCE:  "Méthode de calcul numérique- Tome 2 - Programmes en     *
*              Basic et en Pascal By Claude Nowakowski, Edition du     *
*              P.S.I., 1984" [BIBLI 04].                               *
*                                                                      *
*                                 TPW Release By J-P Moreau, Paris.    *
*                                         (www.jpmoreau.fr)            *
***********************************************************************}
Program Laplace;
Uses WinCrt;

Const

                   {    y   T=x^3            }
      IT = 40;     {     ------------        }           
      ER = 0.01;   {    |     C      |       }
      OM = 1.45;   {    |            |       }
      NM = 8;      { T=0|           C|T=16x  }
      C = 4;       {    |            |       }
                   {    |            |       }
                   {     ------------ x      }
                   {         T=0             }

Var
      dx,r,rm : REAL;
      i,j,l,n1: INTEGER;

      T: Array[0..8,0..8] of REAL;

{calculate Residual at current point (i,j) }
Procedure Schema;
Begin
  r:=T[i+1,j]+T[i-1,j]+T[i,j+1]+T[i,j-1]-4*T[i,j]
End;


{main program}
BEGIN
  dx:=C/NM;
  {Begin limit conditions}
  for i:=0 to NM do
    for j:=0 to NM do
      T[i,j]:=0.0;
  for i:=0 to NM do
  begin
    T[i,NM]:=dx*i*dx*i*dx*i;
    T[NM,i]:=16*dx*i
  end;
  {End limit conditions}
  n1:=NM-1; l:=0;
  Repeat
    rm:=0.0; l:=l+1;
    for i:=1 to n1 do
      for j:=1 to n1 do
      begin
        Schema;
        T[i,j]:=T[i,j]+0.25*OM*r;
        if ABS(r) > rm then rm:=ABS(r)
      end;
  Until (rm < ER) or (l > IT);
  if l<IT then
  begin
    writeln;
    writeln(' Niter  Max. Residual     W');
    writeln(l:4, '     ', rm:-1, '    ',OM:6:2);
    writeln(CHR(29),' Temperature');
    for j:=NM downto 0 do
    begin
      for i:=0 to NM do write(T[i,j]:7:3);
      writeln
    end
  end
  else
    write('No Convergence, RMAXI=', rm:-1);
  Writeln;
  ReadKey;
  DoneWinCrt
END.

{end of file Laplace.pas}