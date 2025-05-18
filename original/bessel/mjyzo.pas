{********************************************************************
*   Compute the zeros of Bessel functions Jn(x), Yn(x), and their   *
*   derivatives using subroutine JYZO                               *
* ----------------------------------------------------------------- *
* SAMPLE RUN:                                                       *
* (Compute 10 zeroes for n=1).                                      *
*                                                                   *
*  Please enter order n and number of zeroes: 1 10                  *
*                                                                   *
*  Zeros of Bessel functions Jn(x), Yn(x) and their derivatives     *
*                       ( n = 1 )                                   *
*   m       jnm           j'nm          ynm           y'nm          *
*  -----------------------------------------------------------      *
*   1     3.8317060     1.8411838     2.1971413     3.6830229       *
*   2     7.0155867     5.3314428     5.4296810     6.9415000       *
*   3    10.1734681     8.5363164     8.5960059    10.1234047       *
*   4    13.3236919    11.7060049    11.7491548    13.2857582       *
*   5    16.4706301    14.8635886    14.8974421    16.4400580       *
*   6    19.6158585    18.0155279    18.0434023    19.5902418       *
*   7    22.7600844    21.1643699    21.1880689    22.7380347       *
*   8    25.9036721    24.3113269    24.3319426    25.8843146       *
*   9    29.0468285    27.4570506    27.4752950    29.0295758       *
*  10    32.1896799    30.6019230    30.6182865    32.1741182       *
*  -----------------------------------------------------------      *
*                                                                   *
* ----------------------------------------------------------------- *
* Ref.: www.esrf.fr/computing/expg/libraries/smf/PROGRAMS/MJYZO.FOR *
*                                                                   *
*                             TPW Release 1.0 By J-P Moreau, Paris. *
*                                       (www.jpmoreau.fr)           *
********************************************************************}
      PROGRAM MJYZO;

      USES Wincrt;

      Type
           VEC = Array[0..30] of Double;
      Var
           RJ0, RJ1, RY0, RY1: VEC;
           M, N, NT: Integer;

      FUNCTION Log10(x:Double): Double;
      begin
        if x>0 then
          log10:= ln(x)/ln(10)
        else
          log10:= -1e15
      end;

      {y^x}
      Function Power(y,x:Double):Double;
      Begin
        IF x<0 THEN EXIT;
        Power:=Exp(x*Ln(y))
      End;

      Procedure JYNDD(N:Integer; X: Double; Var BJN,DJN,FJN,BYN,DYN,FYN:Double);
      Forward;

      Procedure JYZO(N,NT: integer; Var RJ0,RJ1,RY0,RY1:VEC);
{     ========================================================
!       Purpose: Compute the zeros of Bessel functions Jn(x),
!                Yn(x), and their derivatives
!       Input :  n  --- Order of Bessel functions (0 to 100)
!                NT --- Number of zeros (roots)
!       Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT
!                RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT
!                RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT
!                RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT
!       Routine called: JYNDD for computing Jn(x), Yn(x), and
!                       their first and second derivatives
!     ======================================================== }
      Label 10,15,20,25;
      Var
          BJN,DJN,FJN,BYN,DYN,FYN: Double;
          X,X0:Double;
          L:Integer;
      Begin
        IF  N <= 20 THEN
          X:=2.82141+1.15859*N
        ELSE
          X:=N+Power(1.85576*N,0.33333)+Power(1.03315/N,0.33333);
        L:=0;
10:     X0:=X;
        JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN);
        X:=X-BJN/DJN;
        IF ABS(X-X0) > 1e-9 then GOTO 10;
        L:=L+1;
        RJ0[L]:=X;
        X:=X+3.1416+(0.0972+0.0679*N-0.000354*N*N)/L;
        IF L < NT then GOTO 10;
        IF N <= 20 THEN
           X:=0.961587+1.07703*N
        ELSE
           X:=N+Power(0.80861*N,0.33333)+Power(0.07249/N,0.33333);
        IF N = 0 then X:=3.8317;
        L:=0;
15:     X0:=X;
        JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN);
        X:=X-DJN/FJN;
        IF ABS(X-X0) > 1e-9 then GOTO 15;
        L:=L+1;
        RJ1[L]:=X;
        X:=X+3.1416+(0.4955+0.0915*N-0.000435*N*N)/L;
        IF L < NT then GOTO 15;
        IF N <= 20 THEN
           X:=1.19477+1.08933*N
        ELSE
           X:=N+Power(0.93158*N,0.33333)+Power(0.26035/N,0.33333);
        L:=0;
20:     X0:=X;
        JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN);
        X:=X-BYN/DYN;
        IF ABS(X-X0) > 1e-9 then GOTO 20;
        L:=L+1;
        RY0[L]:=X;
        X:=X+3.1416+(0.312+0.0852*N-0.000403*N*N)/L;
        IF L < NT then  GOTO 20;
        IF N <= 20 THEN
           X:=2.67257+1.16099*N
        ELSE
           X:=N+Power(1.8211*N,0.33333)+Power(0.94001/N,0.33333);
        L:=0;
25:     X0:=X;
        JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN);
        X:=X-DYN/FYN;
        IF ABS(X-X0) > 1e-9 then  GOTO 25;
        L:=L+1;
        RY1[L]:=X;
        X:=X+3.1416+(0.197+0.0643*N-0.000286*N*N)/L;
        IF L < NT then GOTO 25
      End;

      Procedure JYNDD(N:Integer; X: Double; Var BJN,DJN,FJN,BYN,DYN,FYN:Double);
{     =========================================================
!       Purpose: Compute Bessel functions Jn(x) and Yn(x), and
!                their first and second derivatives 
!       Input:   x   ---  Argument of Jn(x) and Yn(x) ( x > 0 )
!                n   ---  Order of Jn(x) and Yn(x)
!       Output:  BJN ---  Jn(x)
!                DJN ---  Jn'(x)
!                FJN ---  Jn"(x)
!                BYN ---  Yn(x)
!                DYN ---  Yn'(x)
!                FYN ---  Yn"(x)
!     ========================================================= }
      Label 15;
      Var  
          BJ, BY: VEC;
          K,M, MT, NT: Integer;
          F,F0,F1,BS,SU: Double;
          E0,EC,S1: Double;
      Begin
        For NT:=1 to 900 do
        begin
          MT:=Round(0.5*Log10(6.28*NT)-NT*Log10(1.36*ABS(X)/NT));
          IF MT > 20 then GOTO 15;
        end;
15:     M:=NT;
        BS:=0.0;
        F0:=0.0;
        F1:=1e-35;
        SU:=0.0;
        For K:=M Downto 0 do
        begin
           F:=2.0*(K+1.0)*F1/X-F0;
           IF K <= N+1 then BJ[K+1]:=F;
           IF K MOD 2 = 0 THEN
           begin
             BS:=BS+2.0*F;
             IF K <> 0 then
               if (K Div 2) MOD 2 = 0 then
                 SU:=SU + F/K
               else
                 SU:=SU - F/K;
           end;
           F0:=F1;
           F1:=F
        end;
        For K:=0 to N+1 do BJ[K+1]:=BJ[K+1]/(BS-F);
        BJN:=BJ[N+1];
        EC:=0.5772156649015329;
        E0:=0.3183098861837907;
        S1:=2.0*E0*(LN(X/2.0)+EC)*BJ[1];
        F0:=S1-8.0*E0*SU/(BS-F);
        F1:=(BJ[2]*F0-2.0*E0/X)/BJ[1];
        BY[1]:=F0;
        BY[2]:=F1;
        for K:=2 to N+1 do
        begin
          F:=2.0*(K-1.0)*F1/X-F0;
          BY[K+1]:=F;
          F0:=F1;
          F1:=F
        end;
        BYN:=BY[N+1];
        DJN:=-BJ[N+2]+N*BJ[N+1]/X;
        DYN:=-BY[N+2]+N*BY[N+1]/X;
        FJN:=(N*N/(X*X)-1.0)*BJN-DJN/X;
        FYN:=(N*N/(X*X)-1.0)*BYN-DYN/X
      End;

    {main program}
    BEGIN

      Writeln;
      Write('  Please enter order n and number of zeroes: ');
      Readln(N, NT);
      Writeln;

      JYZO(N,NT,RJ0,RJ1,RY0,RY1);

      WRITELN('  Zeros of Bessel functions Jn(x), Yn(x), and their derivatives');
      WRITELN('                       ( n = ',n,' )');                 
      WRITELN('   m       jnm           j''nm          ynm          y''nm');
      WRITELN('  -----------------------------------------------------------');
      for M:=1 to NT do
        WRITELN(' ',M:3,RJ0[M]:14:7,RJ1[M]:14:7,RY0[M]:14:7,RY1[M]:14:7);
      WRITELN('  -----------------------------------------------------------');
      ReadKey;
      DoneWinCrt

    END.

{ End of file mjyzo.pas }