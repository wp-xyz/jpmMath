{**************************************************************
*   Calculate the median value of an array without sorting    *
* ----------------------------------------------------------- *
* REFERENCE:                                                  *
*      "NUMERICAL RECIPES by W.H. Press, B.P. Flannery,       *
*       S.A. Teukolsky and W.T. Vetterling, Cambridge         *
*       University Press, 1986".                              *
* ----------------------------------------------------------- *
* SAMPLE RUN:  (see input file tmdian1.dat).                  *
*                                                             *
* N=80                                                        *
* Given table:                                                *
*                                                             *
* 407.8 192.8 851.1 604.4 932.3 799.4 914.5 965.8 453.7 295.1 *
* 154.5 977.4 410.2 916.2 934.7 504.8 823.2 225.2 456.6  49.0 *
* 933.5 663.0 335.3 346.6 568.7 956.1 654.7 300.7 379.6 591.9 *
* 992.9 689.6 644.7 305.4 148.2 257.2 664.6 612.1 713.0  99.7 *
*  46.5 167.6 984.6 847.2  55.4  82.7 999.0  10.7 877.7 929.4 *
* 398.1 972.8 874.1 755.1 472.1 122.8 671.4  35.5 128.8  76.8 *
* 454.2 959.2 510.1 791.3 122.8 176.6 237.9 995.8 548.3 309.8 *
* 162.6 996.5 750.0 250.6 577.7 761.1 101.9 797.1 539.0 723.5 *
*                                                             *
* Median value: 558.5000                                      *
*                                                             *
*                        Pascal Release By J-P Moreau, Paris. *
*                                (www.jpmoreau.fr)            *
**************************************************************}
PROGRAM TMDIAN1;
Uses WinCrt;

Const NMAX = 1024;          {max. item number}

Type  pVEC = ^VEC;
       VEC = Array[1..NMAX] of Double;

Var   A:pVEC;              {pointer to input random Table}
      I,N:Integer;
      AMED:Double;
      fp:TEXT;


      Function MAX(a,b:Double):Double;
      Begin
        if a>=b then MAX := a
        else MAX:=b
      End;

      Function MIN(a,b:Double):Double;
      Begin
        if a<=b then MIN := a
        else MIN:=b
      End;


{******************************************************
* Given an array X of N numbers, returns their median *
* value XMED. The array X is not modified, and is     *
* accessed sequentially in each consecutive pass.     *
******************************************************}
Procedure MDIAN1(X:pVEC; N:Integer; Var XMED:Double);
Label 1;
Const BIG:Double = 1.2e16; AFAC:Double = 1.5; AMP:Double = 1.5;
{ Here, AMP is an overconvergence factor: on each iteration,
  we move the guess by this factor. AFAC is a factor used to
  optimize the size of the "smoothing constant" EPS at each
  iteration.            }
Var A,AA,AM,AP,DUM,EPS,SUM,SUMX,XM,XP,XX:Double;
    J,NM,NP:Integer;
Begin
  A:=0.5*(X^[1]+X^[N]);
  EPS:=ABS(X^[N]-X^[1]);
  AP:=BIG;
  AM:=-BIG;
1:SUM:=0.0;
  SUMX:=0.0;
  NP:=0;
  NM:=0;
  XP:=BIG;
  XM:=-BIG;
  For J:=1 to N do
  begin
    XX:=X^[J];
    if XX <> A then
    begin
      if XX > A then
      begin
        NP:=NP+1;
	if XX < XP then XP:=XX
      end
      else if XX < A then
      begin
        NM:=NM+1;
	if XX > XM then XM:=XX
      end;
      DUM:=1.0/(EPS+ABS(XX-A));
      SUM:=SUM+DUM;
      SUMX:=SUMX+XX*DUM
    end;
  end;
  if NP-NM >= 2 then       {guess is too low, make another pass}
  begin
    AM:=A;
    AA:=XP+MAX(0.0,SUMX/SUM-A)*AMP;
    if AA > AP then AA:=0.5*(A+AP);
    EPS:=AFAC*ABS(AA-A);
    A:=AA;
    goto 1
  end
  else if NM-NP >= 2 then  {guess is too high}
  begin
    AM:=A;
    AA:=XM+MIN(0.0,SUMX/SUM-A)*AMP;
    if AA < AM then AA:=0.5*(A+AM);
    EPS:=AFAC*ABS(AA-A);
    A:=AA;
    goto 1
  end
  else                     {guess is ok}
  begin
    if (N MOD 2) = 0 then  {for even N median is an average}
      if NP = NM then
	XMED:=0.5*(XP+XM)
      else if NP > NM then
	XMED:=0.5*(A+XP)
      else
	XMED:=0.5*(XM+A)
    else
      if NP = NM then
	XMED:=A
      else if NP > NM then
	XMED:=XP
      else
	XMED:=XM
  end
End;

{write table of size N to standard output}
Procedure TWRIT(N:Integer; ARR:pVEC);
Var i: Integer;
Begin
  writeln;
  For i:=1 to N do
  begin
    write(' ',ARR^[i]:6:1);
    if (i MOD 10) = 0 then writeln
  end;
End;


{main program}
BEGIN

  New(A);
  Assign(fp,'tmdian1.dat'); Reset(fp);
  readln(fp,N);
  writeln(' N=', N);
  For i:=1 to N do
  begin
    read(fp,A^[i]);
    if (i mod 10) = 0 then readln(fp)
  end;
  close(fp);

  writeln(' Given table:');
  TWRIT(N,A);
 
  {call MDIAN1 subroutine}
  MDIAN1(A,N,AMED);

  writeln;
  writeln(' Median value: ', AMED:8:4);
  writeln;
  ReadKey;

  Dispose(A);
  DoneWinCrt

END.

{end of file tmdian1.pas}