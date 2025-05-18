{*********************************************************
* Calculate the mean and first third moments of a set of *
* data Y(i)                                              *
* ------------------------------------------------------ *
* Description:                                           *
*                                                        *
* Input data:                                            *
*		ndata:	Integer size of data set                 *
*		Y:	Real Vector of ndata ordinates               *
* Output data:                                           *
*		S1:	ndata-Mean of ndata values                   *
*       S2:	Sum (Y(i)-S1)^2 for i=1..ndata               *
*		S3:	Sum (Y(i)-S1)^3 for i=1..ndata               *
*		S4:	Sum (Y(i)-S1)^4 for i=1..ndata               *
* ------------------------------------------------------ *
* SAMPLE RUN:                                            *
*  Number of data: 5                                     *
*  1: 12                                                 *
*  2: 9                                                  *
*  3: 7                                                  *
*  4: 15                                                 *
*  5: 6                                                  *
*                                                        *
*  S1=  9.80000000000291E+0000                           *
*  S2=  5.48000000000466E+0001                           *
*  S3=  7.39200000003912E+0001                           *
*  S4=  1.02497600000165E+0003                           *
*                                                        *
*  Error code: 0                                         *
*                                                        *
* ------------------------------------------------------ *
* Ref.: "Journal of Applied Statistics (1972) vol.21,    *
*        page 226".                                      *
*                                                        *
*                   Pascal Release By J-P Moreau, Paris. *
*                            (www.jpmoreau.fr)           *
*********************************************************} 
PROGRAM TMOMNTS;
Uses WinCrt;

Const NMAX = 25;

Var
    error,i,ndata: Integer;
    Y: Array[1..NMAX] of REAL;
    S1,S2,S3,S4: REAL;



    Procedure MOMNTS(X:REAL; K, N:Integer; Var S1, S2, S3, S4:REAL; Var IFAULT:Integer);
{ ------------------------------------------------------------------
         ALGORITHM AS 52  APPL. STATIST. (1972) VOL.21, P.226

         ADDS A NEW VALUE, X, WHEN CALCULATING A MEAN, AND SUMS
         OF POWERS OF DEVIATIONS. N IS THE CURRENT NUMBER OF
         OBSERVATIONS, AND MUST BE SET TO ZERO BEFORE FIRST ENTRY
  ------------------------------------------------------------------ }
    Label 10,20,30,40,50,60,Return;
    Var
         AN, AN1, DX, DX2, ZERO, ONE, TWO, THREE, FOUR, SIX: REAL;
    Begin

      ZERO:=0.0; ONE:=1.0; TWO:=2.0; THREE:=3.0; FOUR:=4.0; SIX:=6.0;

      IF (K > 0) AND (K < 5) AND (N >= 0) THEN GOTO 10;
      IFAULT := 1;
      goto RETURN;
10:   IFAULT := 0;
      N := N + 1;
      IF N > 1 THEN GOTO 20;

{     FIRST ENTRY, SO INITIALISE }

      S1 := X;
      S2 := ZERO;
      S3 := ZERO;
      S4 := ZERO;
      GOTO RETURN;

{     SUBSEQUENT ENTRY, SO UPDATE }

20:   AN := N;
      AN1 := AN - ONE;
      DX := (X - S1) / AN;
      DX2 := DX * DX;
      Case K of
        1:goto 60;   {calculate only mean S1}
        2:goto 50;   {calculate S1 and S2}
        3:goto 40;   {calculate S1, S2 and S3}
        4:goto 30    {calculate S1..S4}
      End;
30:   S4 := S4 - DX * (FOUR * S3 - DX * (SIX * S2 + AN1 *
           (ONE + AN1*AN1*AN1) * DX2));
40:   S3 := S3 - DX * (THREE * S2 - AN * AN1 * (AN - TWO) * DX2);
50:   S2 := S2 + AN * AN1 * DX2;
60:   S1 := S1 + DX;

RETURN: End;


{main program}
BEGIN

  writeln;
  write('  Number of data: '); readln(ndata);
  For i:=1 to ndata do
  begin
    write(' ',i:2,': '); readln(Y[i])
  end;

  For i:=1 to ndata do
    MOMNTS(Y[i],4,i-1,S1,S2,S3,S4,error);

  writeln;
  writeln('  S1= ', S1);
  writeln('  S2= ', S2);
  writeln('  S3= ', S3);
  writeln('  S4= ', S4);
  writeln;
  writeln('  Error code: ', error);

  ReadKey;
  DoneWinCrt

END.

{end of file momnts.pas}