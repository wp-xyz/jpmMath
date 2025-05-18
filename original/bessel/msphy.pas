{**************************************************************
*      Purpose: This program computes the spherical Bessel    *
*               functions yn(x) and yn'(x) using subroutine   *
*               SPHY                                          *
*      Input :  x --- Argument of yn(x) ( x > 0 )             *
*               n --- Order of yn(x) ( n = 0 to 250 )         *
*      Output:  SY(n) --- yn(x)                               *
*               DY(n) --- yn'(x)                              *
*      Example:   x = 10.0                                    *
*                 n          yn(x)               yn'(x)       *
*               --------------------------------------------  *
*                 0     .8390715291D-01    -.6279282638D-01   *
*                 1     .6279282638D-01     .7134858763D-01   *
*                 2    -.6506930499D-01     .8231361788D-01   *
*                 3    -.9532747888D-01    -.2693831344D-01   *
*                 4    -.1659930220D-02    -.9449751377D-01   *
*                 5     .9383354168D-01    -.5796005523D-01   *
* ----------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special     *
*             Functions,                                      *
*             jin.ece.uiuc.edu/routines/routines.html".       *
*                                                             *
*                        Pascal Release By J-P Moreau, Paris. *
*                                (www.jpmoreau.fr)            *
**************************************************************}
Program MSPHY;

Uses WinCrt;

Type VEC = Array[0..250] of double;

Var
     SY, DY: VEC;
     k, n, nm, ns: Integer;
     X: double;


Procedure SPHY(N:integer; X:double; Var NM:integer; Var SY, DY: VEC);
{       ======================================================
!       Purpose: Compute spherical Bessel functions yn(x) and
!                their derivatives
!       Input :  x --- Argument of yn(x) ( x ò 0 )
!                n --- Order of yn(x) ( n = 0,1,úúú )
!       Output:  SY(n) --- yn(x)
!                DY(n) --- yn'(x)
!                NM --- Highest order computed
!       ====================================================== }
Label 10, 20;
Var
        K:Integer; F, F0, F1:Double;
Begin
	NM:=N;
        if X < 1E-60 then
        begin
           for K:=0 to N do
           begin
             SY[K]:=-1E+60;
             DY[K]:=1E+60
           end;
           goto 10
        end;
        SY[0]:=-cos(X)/X;
        SY[1]:=(SY[0]-sin(X))/X;
        F0:=SY[0];
        F1:=SY[1];
        for K:=2 to N do
        begin
           F:=(2.0*K-1.0)*F1/X-F0;
           SY[K]:=F;
           if abs(F) >= 1E+60 then goto 20;              
           F0:=F1;
           F1:=F
        end;
20:     NM:=K-1;
        DY[0]:=(sin(X)+cos(X)/X)/X;
        for K:=1 to NM do
           DY[K]:=SY[K-1]-(K+1.0)*SY[K]/X;
10: End;


{main program}
BEGIN
        InitWinCrt;
        writeln;
        write(' Please enter n and x: ');
        readln(n, X);

        Inc(n);

	if n <= 10 then
           ns:=1
        else
        begin
           write(' Please enter order step Ns: ');
           readln(ns)
        end;
        
	SPHY(n,X,nm,SY,DY);
        
	writeln;
        writeln('  n          yn(x)             yn''(x)');
        writeln('--------------------------------------------');
        k:=0;
        Repeat 
          writeln(k:3, SY[k]:18:6, DY[k]:18:6);
          k:=k+ns
        until k = nm+1;
        ReadKey;
        DoneWinCrt
END.

{end of file msphy.pas}