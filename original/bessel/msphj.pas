{***************************************************************
*     Purpose: This program computes the spherical Bessel      *
*              functions jn(x) and jn'(x) using subroutine     *
*              SPHJ                                            *
*     Input :  x --- Argument of jn(x)                         *
*              n --- Order of jn(x)  ( n = 0 to 250 )          * 
*     Output:  SJ(n) --- jn(x)                                 *
*              DJ(n) --- jn'(x)                                * 
*     Example:   x =10.0                                       *
*                n          jn(x)              jn'(x)          *
*              --------------------------------------------    *
*                0    -.5440211109D-01    -.7846694180D-01     * 
*                1     .7846694180D-01    -.7009549945D-01     *
*                2     .7794219363D-01     .5508428371D-01     *
*                3    -.3949584498D-01     .9374053162D-01     *
*                4    -.1055892851D+00     .1329879757D-01     *
*                5    -.5553451162D-01    -.7226857814D-01     *
* ------------------------------------------------------------ *
* REFERENCE: "Fortran Routines for Computation of Special      *
*             Functions,                                       *
*             jin.ece.uiuc.edu/routines/routines.html".        *
*                                                              *
*                         Pascal Release By J-P Moreau, Paris. *
*                                 (www.jpmoreau.fr)            *
***************************************************************}
Program MSPHJ;

Uses WinCrt;

Type VEC = Array[0..250] of double;

Const MM = 2.30258509;

Var
     SJ, DJ: VEC;
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

Procedure SPHJ(N:integer; X:double; Var NM:integer; Var SJ, DJ: VEC);
{      =======================================================
!      Purpose: Compute spherical Bessel functions jn(x) and
!               their derivatives
!      Input :  x --- Argument of jn(x)
!               n --- Order of jn(x)  ( n = 0,1,úúú )
!      Output:  SJ(n) --- jn(x)
!               DJ(n) --- jn'(x)
!               NM --- Highest order computed
!      Routines called:
!               MSTA1 and MSTA2 for computing the starting
!               point for backward recurrence
!      ======================================================= }
Label   10;
Var
        K, M: Integer; CS, F, F0, F1, SA, SB: Double;
Begin
        NM:=N;
        if abs(X) < 1e-100 then
        begin
	   for K:=0 to N do
           begin
             SJ[K]:=0.0;
             DJ[K]:=0.0
           end; 
           SJ[0]:=1.0;
           DJ[1]:=0.333333333333333;
           goto 10;
        end;
        SJ[0]:=sin(X)/X;
        SJ[1]:=(SJ[0]-cos(X))/X;
        if N >= 2 then
        begin
	   SA:=SJ[0];
           SB:=SJ[1];
           M:=MSTA1(X,200);
           if M < N then
              NM:=M
           else
              M:=MSTA2(X,N,15);
           F0:=0.0;
           F1:=1.0-100;
           for K:=M Downto 0 do
           begin
              F:=(2.0*K+3.0)*F1/X-F0;
              if K <= NM then SJ[K]:=F;
              F0:=F1;
              F1:=F
           end;
           if abs(SA) > abs(SB) then CS:=SA/F;
           if abs(SA) <= abs(SB) then CS:=SB/F0;
	   for K:=0 to NM do  SJ[K] := SJ[K]*CS
        end;      
        DJ[0]:=(cos(X)-sin(X)/X)/X;
	for K:=1 to NM do
          DJ[K]:=SJ[K-1]-(K+1.0)*SJ[K]/X;
10:End;


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
        
	SPHJ(n,X,nm,SJ,DJ);
        
	writeln;
        writeln('  n          jn(x)             jn''(x)');
        writeln('--------------------------------------------');
        k:=0;
        Repeat 
          writeln(k:3, SJ[k]:18:6, DJ[k]:18:6);
          k:=k+ns
        until k = nm+1;
        ReadKey;
        DoneWinCrt
END.

{end of file msphj.pas}