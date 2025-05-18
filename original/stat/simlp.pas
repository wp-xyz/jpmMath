{**********************************************************
* Simple Linear Regression Minimizing the Sum of Absolute *
* Deviation (SAD)                                         *
* ------------------------------------------------------- *
* SAMPLE RUN:                                             *
* (Find Line Alpha + Beta X, minimizing SAD for next set  *
*  of data):                                              *
*                                                         *
* N = 10                                                  *
* X(1)=1  Y(1)=1.15                                       *
* X(2)=2  Y(2)=2.22                                       *
* X(3)=3  Y(3)=2.94                                       *
* X(4)=4  Y(4)=3.85                                       *
* X(5)=5  Y(5)=5.10                                       *
* X(6)=6  Y(6)=5.99                                       *
* X(7)=7  Y(7)=7.25                                       *
* X(8)=8  Y(8)=8.04                                       *
* X(9)=9  Y(9)=9.003                                      *
* X(10)=10 Y(10)=9.999                                    *
*                                                         *
*  Alpha =  0.1667778                                     *
*  Beta  =  0.9832222                                     *
*                                                         *
*  SAD =  0.8270000                                       *
*  ITER = 3                                               *
*  Error code: 0                                          *
*                                                         *
* ------------------------------------------------------- *
* Ref.: Journal of Applied Statistics (1978) vol.27 no.3. *
*                                                         *
*                   Pascal Release By J-P Moreau, Paris.  *
*                            (www.jpmoreau.fr)            *
**********************************************************} 
PROGRAM TSIMLP;
Uses WinCrt;

Const
      NMAX = 25;
      ACU: Double = 1E-6;
      BIG: Double = 1E19;
      HALF: Double = 0.5;
      ZERO: Double = 0.0;
      ONE: Double = 1.0;
      TWO: Double = 2.0;

Type
      pVEC = ^VEC;
       VEC = Array[1..NMAX] of Double;
     pIVEC = ^IVEC;
      IVEC = Array[1..NMAX] of Integer;

Var
      N:Integer;         {Number of points}
      X, Y: pVEC;        {Input Data Set}
      D: pVEC;           {Utility vector}
      INEXT: pIVEC;      {    idem      }
      SAD:Double;        {Sum of Absolute Deviation}
      Alpha,Beta:Double; {Coefficients of linear regression}
      Iter:Integer;      {Number of iterations}
      Ifault:Integer;    {Error code (must be zero) }
      I:Integer;         {loop index}

    Function ISign(a,b : Integer) : Integer;
    Begin
      if (b <0.0) then ISign := - Abs(a)
                  else ISign :=   Abs(a)
    End;

    Function Sign(a,b : Double) : Double;
    Begin
      if (b <0.0) then Sign := - Abs(a)
                  else Sign :=   Abs(a)
    End;


    Procedure SIMLP(N:Integer; X, Y:pVEC; Var SAD, ALPHA, BETA:Double; D: pVEC;
                    Var ITER:Integer; INEXT:pIVEC; Var IFAULT:Integer);

{     ALGORITHM AS 132  APPL. STATIST. (1978) VOL.27, NO.3

      SIMPL:   Fit  Y := ALPHA + BETA.X + error             }

    Label 10,20,50,60,70,80,90,100,110,120,130,140,150,160, Return;
    Var A1,A2,AAAA,BBBB,AHALF,AONE,DDD,DET,TOT1,TOT2,Y1,Y2: Double;
        AAA,BBB,RATIO,RHO,SUBT,SUM,T,TEST,ZZZ: Double;
        I,IBAS1,IBAS2,IFLAG,IIN,IOUT,ISAVE,J: Integer;
    Begin
{     Initial settings }

      IFAULT := 0;
      ITER := 0;
      AHALF := HALF + ACU;
      AONE := AHALF + AHALF;

{     Determine initial basis }

      D^[1] := ZERO;
      Y1 := Y^[1];
      IBAS1 := 1;
      A1 := X^[1];
      For I := 2 to N do
      begin
	IF ABS(A1 - X^[I]) < ACU Then GOTO 10;
	A2 := X^[I];
	IBAS2 := I;
	Y2 := Y^[I];
	GOTO 20;
10:   end;
      IFAULT := 1;
      GOTO RETURN;

{     Calculate initial beta value }

20:   DET := ONE / (A2 - A1);
      AAAA := (A2 * Y1 - A1 * Y2) * DET;
      BBBB := (Y2 - Y1) * DET;

{     Calculate initial D-vector }

      For I := 2 to N do
      begin
	DDD := Y^[I] - (AAAA + BBBB * X^[I]);
	D^[I] := Sign(ONE, DDD)
      end;
      TOT1 := ONE;
      TOT2 := X^[IBAS2];
      D^[IBAS2] := - ONE;
      For I := 2 to N do
      begin
	TOT1 := TOT1 + D^[I];
	TOT2 := TOT2 + D^[I] * X^[I]
      end;
      T := (A2 * TOT1 - TOT2) * DET;
      IF ABS(T) < AONE Then GOTO 50;
      DET := - DET;
      GOTO 70;

{     Main iterative loop begins }

50:   T := (TOT2 - A1 * TOT1) * DET;
      IF ABS(T) < AONE Then GOTO 160;
      IFLAG := 2;
      IOUT := IBAS2;
      X^[IOUT] := A1;
      AAA := A1;
      BBB := A2;
      GOTO 80;
60:   T := (TOT2 - A2 * TOT1) * DET;
      IF ABS(T) < AONE Then GOTO 160;
70:   IFLAG := 1;
      BBB := A1;
      AAA := A2;
      IOUT := IBAS1;
80:   RHO := Sign(ONE, T);
      T := HALF * ABS(T);
      DET := DET * RHO;

{     Perform partial sort of ratios }

      INEXT^[IBAS1] := IBAS2;
      RATIO := BIG;
      SUM := AHALF;
      For I := 1 to N do
      begin
	DDD := (X^[I] - AAA) * DET;
	IF DDD * D^[I] <= ACU Then  GOTO 120;
	TEST := (Y^[I] - AAAA - BBBB * X^[I]) / DDD;
	IF TEST >= RATIO Then GOTO 120;
	J := IBAS1;
	SUM := SUM + ABS(DDD);
90:     ISAVE := ABS(INEXT^[J]);
	IF TEST >= D^[ISAVE] Then GOTO 110;
	IF SUM < T Then GOTO 100;
	SUBT := ABS((X^[ISAVE] - AAA) * DET);
	IF SUM - SUBT <  T Then GOTO 100;
	SUM := SUM - SUBT;
	D^[ISAVE] := Sign(1, INEXT^[J]);
	INEXT^[J] := INEXT^[ISAVE];
	GOTO 90;
100:    J := ISAVE;
	ISAVE := ABS(INEXT^[J]);
	IF TEST < D^[ISAVE] Then GOTO 100;
110:    INEXT^[I] := INEXT^[J];
	INEXT^[J] := ISign(I, Round(D^[I]));
	D^[I] := TEST;
	IF SUM <  T Then GOTO 120;
	IIN := ABS(INEXT^[IBAS1]);
	RATIO := D^[IIN];
120:  end;

{     Update basic indicators }

      IIN := ABS(INEXT^[IBAS1]);
      J := IIN;
130:  ISAVE := ABS(INEXT^[J]);
      IF ISAVE = IBAS2 Then GOTO 140;
      ZZZ := ISign(1, INEXT^[J]);
      TOT1 := TOT1 - ZZZ - ZZZ;
      TOT2 := TOT2 - TWO * ZZZ * X^[ISAVE];
      D^[ISAVE] := - ZZZ;
      J := ISAVE;
      GOTO 130;
140:  ZZZ := ISign(1, INEXT^[IBAS1]);
      TOT1 := TOT1 - RHO - ZZZ;
      TOT2 := TOT2 - RHO * BBB - ZZZ * X^[IIN];
      D^[IOUT] := - RHO;
      ITER := ITER + 1;
      IF IFLAG = 1 Then GOTO 150;
      X^[IBAS2] := A2;
      IBAS2 := IIN;
      D^[IBAS2] := - ONE;
      A2 := X^[IIN];
      Y2 := Y^[IIN];
      DET := ONE / (A1 - A2);
      AAAA := (A1 * Y2 - A2 * Y1) * DET;
      BBBB := (Y1 - Y2) * DET;
      GOTO 60;
150:  IBAS1 := IIN;
      A1 := X^[IIN];
      D^[IBAS1] := ZERO;
      Y1 := Y^[IIN];
      DET := ONE / (A2 - A1);
      AAAA := (A2 * Y1 - A1 * Y2) * DET;
      BBBB := (Y2 - Y1) * DET;
      GOTO 50;

{     Calculate optimal sum of absolute deviations }

160:  SAD := ZERO;
      For I := 1 to N do
      begin
	D^[I] := Y^[I] - (AAAA + BBBB * X^[I]);
	SAD := SAD + ABS(D^[I])
      end;
      ALPHA := AAAA;
      BETA := BBBB;

Return: End;


{main program}
BEGIN

  New(X); New(Y); New(D); New(INEXT);

  N:=10;

  For I:=1 to N do X^[I]:=1.0*I;
  Y^[1]:=1.15; Y^[2]:=2.22; Y^[3]:=2.94; Y^[4]:=3.85; Y^[5]:=5.10;
  Y^[6]:=5.99; Y^[7]:=7.25; Y^[8]:=8.04; Y^[9]:=9.003; Y^[10]:=9.999;

  SIMLP(N,X,Y,SAD,ALPHA,BETA,D,ITER,INEXT,IFAULT);

  writeln;
  writeln(' Alpha = ', ALPHA:10:7);
  writeln(' Beta  = ', BETA:10:7);
  writeln;
  writeln(' SAD = ', SAD:10:7);
  writeln(' ITER = ', ITER);
  writeln(' Error code: ', IFAULT);
  writeln;

  ReadKey;
  Dispose(X); Dispose(Y); Dispose(D); Dispose(INEXT);
  DoneWinCrt  

END.

{end of file simlp.pas}