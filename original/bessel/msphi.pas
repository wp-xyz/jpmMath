{**************************************************************
*     Purpose: This program computes the modified spherical   *
*              Bessel functions of the first kind in(x) and   *
*              in'(x) using subroutine SPHI                   *
* ----------------------------------------------------------- *
*     Input :  x --- Argument of in(x)                        *
*              n --- Order of in(x) ( 0 to 250 )              *
*     Output:  SI(n) --- in(x)                                *
*              DI(n) --- in'(x)                               *
*     Example: x = 10.0                                       *
*                n          in(x)               in'(x)        *
*              --------------------------------------------   *
*                0     .1101323287D+04     .9911909633D+03    *
*                1     .9911909633D+03     .9030850948D+03    *
*                2     .8039659985D+03     .7500011637D+03    *
*                3     .5892079640D+03     .5682828129D+03    *
*                4     .3915204237D+03     .3934477522D+03    *
*                5     .2368395827D+03     .2494166741D+03    *
* ----------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special     *
*             Functions,                                      *
*             jin.ece.uiuc.edu/routines/routines.html".       *
*                                                             *
*                        Pascal Release By J-P Moreau, Paris. *
*                               (www.jpmoreau.fr)             *
**************************************************************}
Program msphi;

Uses WinCrt;

Type VEC = Array[0..250] of double;

Const MM = 2.30258509;

Var
     SI, DI: VEC;
     k, n, nm, ns: Integer;
     X: double;


Function ENVJ(N:integer; X:double): double;
Begin
  ENVJ := (0.5*Ln(6.28*N)-N*Ln(1.36*X/N))/MM
End;

Function MSTA1(X: double; MP: integer): Integer;

{      ===================================================
!      Purpose: Determine the starting point for backward  
!               recurrence such that the magnitude of    
!               Jn(x) at that point is about 10^(-MP)
!      Input :  x     --- Argument of Jn(x)
!               MP    --- Value of magnitude
!      Output:  MSTA1 --- Starting point   
!      =================================================== }
Label 20;
Var
        A0,F,F0,F1:Double; IT,NN,N0,N1:Integer;
Begin
        A0:=abs(X);
        N0:=Round(1.1*A0)+1;
        F0:=ENVJ(N0,A0)-MP;
        N1:=N0+5;
        F1:=ENVJ(N1,A0)-MP;
        for IT:=1 to 20 do
        begin             
           NN:=Round(N1-(N1-N0)/(1.0-F0/F1));                  
           F:=ENVJ(NN,A0)-MP;
           if abs(NN-N1) < 1 then goto 20;
           N0:=N1;
           F0:=F1;
           N1:=NN;
           F1:=F
        end;
20:	MSTA1 := NN
End;

Function MSTA2(X:double; N, MP: integer): Integer;

{      ===================================================
!      Purpose: Determine the starting point for backward
!               recurrence such that all Jn(x) has MP
!               significant digits
!      Input :  x  --- Argument of Jn(x)
!               n  --- Order of Jn(x)
!               MP --- Significant digit
!      Output:  MSTA2 --- Starting point
!      =================================================== }
Label 20;
Var
        A0,EJN,F,F0,F1,HMP,OBJ:Double; IT,N0,N1,NN: Integer;
Begin
        A0:=abs(X);
        HMP:=0.5*MP;
        EJN:=ENVJ(N,A0);
        if EJN <= HMP then
        begin
           OBJ:=MP;
           N0:=Round(1.1*A0)
        end
        else
        begin
           OBJ:=HMP+EJN;
           N0:=N
        end;
        F0:=ENVJ(N0,A0)-OBJ;
        N1:=N0+5;
        F1:=ENVJ(N1,A0)-OBJ;
        for IT:=1 to 20 do
        begin
           NN:=Round(N1-(N1-N0)/(1.0-F0/F1));
           F:=ENVJ(NN,A0)-OBJ;
           if abs(NN-N1) < 1 then goto 20;
           N0:=N1;
           F0:=F1;
           N1:=NN;
           F1:=F
        end;
20:     MSTA2 := NN+10
End;

Function cosh(x:double): double;
begin
  cosh:=(exp(x)+exp(-x))/2.0
end;

Function sinh(x:double): double;
begin
  sinh:=(exp(x)-exp(-x))/2.0
end;


Procedure SPHI(N:integer; X:double; Var NM:integer; Var SI, DI: VEC);

{      ========================================================
!      Purpose: Compute modified spherical Bessel functions
!               of the first kind, in(x) and in'(x)
!      Input :  x --- Argument of in(x)
!               n --- Order of in(x) ( n = 0,1,2,... )
!      Output:  SI(n) --- in(x)
!               DI(n) --- in'(x)
!               NM --- Highest order computed
!      Routines called:
!               MSTA1 and MSTA2 for computing the starting
!               point for backward recurrence
!      ======================================================== }
Label   10;
Var
        K, M: Integer; CS, F, F0, F1, SI0: Double;
Begin		
        NM:=N;
        if abs(X) < 1E-100 then
        begin
	   for K:=0 to N do
           begin
              SI[K]:=0.0;
              DI[K]:=0.0
           end; 
           SI[0]:=1.0;
           DI[1]:=0.333333333333333;
           goto 10
        end;
        SI[0]:=sinh(X)/X;
        SI[1]:=-(sinh(X)/X-cosh(X))/X;
        SI0:=SI[0];
        if N >= 2 then
        begin
           M:=MSTA1(X,200);
           if M < N then
              NM:=M
           else
              M:=MSTA2(X,N,15);
           F0:=0.0;
           F1:=1.0-100;
           for K:=M Downto 0 do
           begin
              F:=(2.0*K+3.0)*F1/X+F0;
              if K <= NM then SI[K]:=F;
              F0:=F1;
              F1:=F
           end;
           CS:=SI0/F;
           for K:=0 to NM do  SI[K] := SI[K]*CS
        end;
        DI[0]:=SI[1];
        for K:=1 to NM do
          DI[K]:=SI[K-1]-(K+1.0)/X*SI[K];
        
10: End;


{main program}
BEGIN
        InitWinCrt;
        writeln;
        write(' Please enter n and x: ');
        readln(n, X);
        
	if n <= 10 then
           ns:=1
        else
        begin
           write(' Please enter order step Ns: ');
           readln(ns)
        end;
        
	SPHI(n,X,nm,SI,DI);
        
	writeln;
        writeln('  n          in(x)             in''(x)');
        writeln('--------------------------------------------');
        k:=0;
        Repeat 
          writeln(k:3, SI[k]:18:6, DI[k]:18:6);
          k:=k+ns
        until k = nm+1;
        ReadKey;
        DoneWinCrt
END.

{end of file msphi.pas}
