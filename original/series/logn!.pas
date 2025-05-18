{***************************************************
*     Program to demonstrate LN(X!) procedure      *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*           Pascal Version by J.-P. Moreau, Paris. *
*                      (www.jpmoreau.fr)           *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
*    X     LN(X!)        EXP(LN(X!))               *
*  ----------------------------------              *
*    1    -0.000307                1               *
*    2     0.693146                2               *
*    3     1.791759                6               *
*    4     3.178054               24               *
*    5     4.787492              120               *
*    6     6.579251              720               *
*    7     8.525161             5040               *
*    8    10.604603            40320               *
*    9    12.801827           362880               *
*   10    15.104413          3628800               *
*   11    17.502308         39916800               *
*   12    19.987214        479001600               *
*   13    22.552164       6227020800               *
*   14    25.191221      87178291200               *
*   15    27.899271    1307674368000               *
*                                                  *  
****************************************************
See file chi-sq.txt, equation (2.3.9)              }

PROGRAM LNFACTX;
Uses WinCrtMy, Typedef, Utilit;

Var  i   : integer;
     x,y : ar_reel;


{**************************************************
* Series approximation subroutine for LN(X!)      *
* Accuracy better then 6 places for x>=3          *
* Accuracy better than 12 places for x>10         *
* Advantage is that very large values of the      *
* argument can be used without fear of over flow. *
* Reference: CRC Math Tables.                     *
* x is the input, y is the output.                *
**************************************************}
PROCEDURE LN_FACTX;
Var x1 : ar_reel;
Begin
  x1 := 1.0 / (x * x);
  y := (x + 0.5) * LN(x) - x * (1 - x1 / 12 + x1 * x1 / 360 - x1 * x1 * x1 / 1260 + x1 * x1 * x1 * x1 / 1680);
  y := y + 0.918938533205
End;


{main program}
BEGIN
  CLRSCR;
  Writeln;
  Writeln('   X       LN(X!)      EXP(LN(X!)) ');
  Writeln(' ----------------------------------');
  FOR i := 1 TO 15 DO
  begin
    x := i;

    LN_FACTX;

    Write('  '); Write(x:2:0); Write('    '); Aff_reel(y); Write('    '); Writeln(EXP(y):13:0)
  end;
  Writeln(' ----------------------------------');
  Writeln;
  ReadKey; DoneWinCrt
END.

{End of file logn!.pas}