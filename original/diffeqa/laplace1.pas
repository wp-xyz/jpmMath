{********************************************************************************************************
*                Solve Laplace Equation by relaxation Method: d2T/dx2 + d2T/dy2 = 0  (2)                *
* ----------------------------------------------------------------------------------------------------- *
* Example #2: Determine the temperatures T(x,y) in a domain [x,y] made of two concentrical rectangles   *
*             (see sketch below) with following limit conditions:                                       *
*             1) T = 0 for sides of inner rectangle.                                                    *
*             2) T = 1 for sides of outer rectangle.                                                    *
*             Main parameters:                                                                          *
*             IM,JM: length, width of outer rectangle.                                                  *
*             IB-IA,JB-JA: length, width of inner rectangle.                                            *
*             Number of subdivisions 16 in  ox and 8 in Oy.                                             *
*             integration steps dx=IM/16, dy = JM/8.                                                    *
*             weight coefficient w = 1.45                                                               *
*             Maximum residual error ER = 0.05                                                          *
*             Maximum number of iterations IT = 40                                                      *
*             (See file laplace.pas for example #1).                                                    *
* ----------------------------------------------------------------------------------------------------- *
* SAMPLE RUN:                                                                                           *
*                                                                                                       *
* Output file laplace.lst contains:                                                                     *
*                                                                                                       *
* Niter  Max. Residual     W                                                                            *
*  10      8.4E-0003      1.45                                                                          *
* Temperature                                                                                           *
* 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 *
*                                                                                                       *
* 1.000 0.965 0.925 0.871 0.796 0.688 0.551 0.514 0.508 0.515 0.551 0.688 0.796 0.871 0.924 0.965 1.000 *
*                                                                                                       *
* 1.000 0.936 0.862 0.764 0.624 0.404 0.000 0.000 0.000 0.000 0.000 0.405 0.625 0.765 0.862 0.936 1.000 *
*                                                                                                       *
* 1.000 0.917 0.822 0.700 0.533 0.304 0.000                   0.000 0.305 0.535 0.700 0.822 0.917 1.000 *
*                                                                                                       *
* 1.000 0.910 0.809 0.679 0.506 0.278 0.000                   0.000 0.280 0.508 0.679 0.808 0.910 1.000 *
*                                                                                                       *
* 1.000 0.919 0.823 0.700 0.534 0.304 0.000                   0.000 0.303 0.535 0.701 0.822 0.917 1.000 *
*                                                                                                       *
* 1.000 0.940 0.865 0.766 0.625 0.404 0.000 0.000 0.000 0.000 0.000 0.412 0.629 0.767 0.863 0.937 1.000 *
*                                                                                                       *
* 1.000 0.965 0.928 0.872 0.797 0.688 0.551 0.514 0.507 0.514 0.550 0.688 0.803 0.874 0.926 0.966 1.000 *
*                                                                                                       *
* 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 *
*                                                                                                       *
* ----------------------------------------------------------------------------------------------------- *
* REFERENCE:  "Méthode de calcul numérique- Tome 2 - Programmes en Basic et en Pascal By Claude         *
*              Nowakowski, Edition du P.S.I., 1984" [BIBLI 04].                                         *
*                                                                                                       *
*                                                                  TPW Release By J-P Moreau, Paris.    *
*                                                                          (www.jpmoreau.fr)            *
********************************************************************************************************}
PROGRAM Laplace1;

Uses WinCrt;

Const               {   Y                                 }
      IT=40;        { JM ---------------------------      }
      ER=0.05;      {   |                           |     }      
      OM=1.45;      { JB|........--------           |     }
      IM=16;        {   |        | T=0  |           |     }
      JM=8;         {   |        |      |           |T=1  }
      IA=6;         { JA|........-------            |     }            
      IB=10;        {   |        .     .            |     }
      JA=2;         {   |        .     .            |     }
      JB=6;         {    --------------------------- X    }
                    {   0       IA    IB           IM     }

Var
      r,rm: REAL;
      i,j,k,l: Integer;
      T:Array[0..16,0..8] of REAL;

      fp: TEXT;

      {calculate residual error at (i,j) }
      Procedure Schema;
      Begin
        r:=T[i+1,j]+T[i-1,j]+T[i,j+1]+T[i,j-1]-4*T[i,j];
        T[i,j]:=T[i,j]+0.25*OM*r;
        if ABS(r) > rm then rm:=r
      End;

      {print current value of temperature array T}
      Procedure Print;
      Begin
        write(fp,T[i,j]:6:3)
      End;

  {main program}
  BEGIN

    Assign(fp,'laplace.lst'); Rewrite(fp);

    {Fix limit conditions}
    For i:=0 to IM do
      For j:=0 to JM do
        T[i,j]:=0.0;
    For i:=0 to IM do
    begin
      T[i,0]:=1.0;
      T[i,JM]:=1.0
    end;
    For j:=0 to JM do
    begin
      T[0,j]:=1.0;
      T[IM,j]:=1.0
    end;
    l:=0;
    {main calculation loop}
    Repeat
      rm:=0.0; Inc(l);
      For i:=1 to IM-1 do
        For j:=1 to JA-1 do
          Schema;
      For j:=JA to JB do
      begin
        For i:=1 to IA-1 do Schema;
        For i:=IB+1 to IM-1 do Schema
      end;
      For i:=1 to IM-1 do
        For j:=JB+1 to JM-1 do
          Schema;
    Until (rm < ER) or (l > IT);
    {print results}
    If l < IT then
    begin
      writeln(fp);
      writeln(fp,' Niter  Max. Residual     W');
      writeln(fp,l:4, '     ', rm:-1, '    ',OM:6:2);
      writeln(fp,' Temperature');
      For J:=JM Downto JB do
      begin
        For i:=0 to IM do Print;
        writeln(fp);
        writeln(fp)
      end;
      For j:=JB-1 Downto JA+1 do
      begin
        For i:=0 to IA do Print;
        For i:=IA+1 to IB-1 do write(fp,' ':6);
        For i:=IB to IM do Print;
        writeln(fp);
        writeln(fp)
      end;
      For J:=JA Downto 0 do
      begin
        For i:=0 to IM do Print;
        writeln(fp);
        writeln(fp)
      end;
    end
    else
      write(fp,' No convergence, R max. =', rm);

    close(fp);
    writeln;
    writeln(' See results in file laplace.lst...');
    ReadKey;
    DoneWinCrt

  END.

{end of file laplace1.pas}