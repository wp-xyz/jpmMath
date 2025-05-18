{*******************************************************
*      Calculate Incomplete Beta Function Ix(a,b)      *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
*                                                      *
*  Incomplete Beta Function Ix(a,b)                    *
*                                                      *
*  A = 5.00000000000000E-0001                          *
*  B = 5.00000000000000E+0000                          *
*  x = 2.00000000000000E-0001                          *
*                                                      *
*  Y = 8.55072399728097E-0001                          *
*                                                      *
* ---------------------------------------------------- *
* Reference:                                           * 
* "Numerical Recipes,  By W.H. Press, B.P. Flannery,   *
*  S.A. Teukolsky and T. Vetterling, Cambridge         *
*  University Press, 1986" [BIBLI 08].                 *
*                                                      *
*                 Pascal Release By J-P Moreau, Paris. *
*                          (www.jpmoreau.fr)           *
*******************************************************}
PROGRAM TEST_BETAI;

Uses WinCrt;

Var
    a,b,x,y: Double;

    Function GAMMLN(XX:Double): Double; Forward;
    Function BETACF(A,B,X:Double): Double; Forward;


Function BETAI(A,B,X:Double): Double;
{ returns the incompleteBeta function Ix(a,b) }
Label Return;
Var  BT: Double;
Begin
  IF (X < 0.0) OR (X > 1.0) Then
  begin
    writeln(' BETAI: Bad argument X (must be 0<=X<=1).');
    goto Return
  end;

  IF (X = 0.0) OR (X = 1.0) Then
    BT:=0.0
  ELSE
    BT:=EXP(GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)+A*Ln(X)+B*Ln(1.0-X));

  IF X < (A+1.0)/(A+B+2.0) THEN
  begin
    BETAI:=BT*BETACF(A,B,X)/A;
    goto RETURN
  end
  ELSE
  begin  
    BETAI:=1.0-BT*BETACF(B,A,1.0-X)/B;
    goto RETURN
  end;
Return: End;

Function BETACF(A,B,X:Double): Double;
{ continued fraction for incomplete Beta function, used by BETAI }
Label 1, Return;
Const ITMAX:Integer = 100;
      EPS:Double = 3e-7;
Var   AM,BM,AZ,QAB,QAP,QAM,BZ,EM,TEM,D,AP,BP,APP,BPP,AOLD: Double;
      M: Integer;
Begin
  AM:=1.0;
  BM:=1.0;
  AZ:=1.0;
  QAB:=A+B;
  QAP:=A+1.0;
  QAM:=A-1.0;
  BZ:=1.0-QAB*X/QAP;
  For M:=1 to ITMAX do
  begin
    EM:=M;
    TEM:=EM+EM;
    D:=EM*(B-M)*X/((QAM+TEM)*(A+TEM));
    AP:=AZ+D*AM;
    BP:=BZ+D*BM;
    D:=-(A+EM)*(QAB+EM)*X/((A+TEM)*(QAP+TEM));
    APP:=AP+D*AZ;
    BPP:=BP+D*BZ;
    AOLD:=AZ;
    AM:=AP/BPP;
    BM:=BP/BPP;
    AZ:=APP/BPP;
    BZ:=1.0;
    IF ABS(AZ-AOLD) < EPS*ABS(AZ) Then GOTO 1
  end;
  writeln(' BETACF: A or B too big, or ITMAX too small.');
  goto RETURN;
1:BETACF:=AZ;
Return: End;
  
Function GAMMLN(XX:Double): Double;
{ returns the value Ln(Gamma(XX) for XX>0 }
Var STP,HALF,ONE,FPF,X,TMP,SER: Double;
    COF: Array[1..6] of Double;
    J: Integer;
Begin
  COF[1]:= 76.18009173; COF[2]:=-86.50532033; COF[3]:=24.01409822;
  COF[4]:=-1.231739516; COF[5]:=0.120858003E-2; COF[6]:=-0.536382E-5;
  STP:=2.50662827465;
  HALF:=0.5; ONE:=1.0; FPF:=5.5;
  X:=XX-ONE;
  TMP:=X+FPF;
  TMP:=(X+HALF)*Ln(TMP)-TMP;
  SER:=ONE;
  For J:=1 to 6 do
  begin
    X:=X+ONE;
    SER:=SER+COF[J]/X
  end;
  GAMMLN:=TMP+Ln(STP*SER)
End;


{main program}
BEGIN

  a:=0.50; b:=5.0; x:=2e-01;

  y := BETAI(a,b,x);

  writeln;
  writeln(' Incomplete Beta Function Ix(A,B)');
  writeln;
  writeln(' A =', a);
  writeln(' B =', b);
  writeln(' x =', x);
  writeln;
  writeln(' Y =', y);
  writeln;

  ReadKey;
  DoneWinCrt

END.

{ end of file ibeta.pas }