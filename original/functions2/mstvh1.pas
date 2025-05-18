{**************************************************************
!*       Purpose: This program computes Struve function       *
!*                H1(x) using subroutine STVH1                *
!*       Input :  x   --- Argument of H1(x) ( x ò 0 )         *
!*       Output:  SH1 --- H1(x)                               *
!*       Example:                                             *
!*                   x          H1(x)                         *
!*                -----------------------                     *
!*                  0.0       .00000000                       *
!*                  5.0       .80781195                       *
!*                 10.0       .89183249                       *
!*                 15.0       .66048730                       *
!*                 20.0       .47268818                       *
!*                 25.0       .53880362                       *
!* ---------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special    *
!*             Functions, jin.ece.uiuc.edu/routines/routines  *
!*             .html".                                        *
!*                                                            *
!*                       Pascal Release By J-P Moreau, Paris. *
!*                                (www.jpmoreau.fr)           *
!*************************************************************} 
Program MSTVH1;

Uses WinCrt;

Var  X, SH1: double;

Procedure STVH1(X: double; Var SH1: double);
{       =============================================
!       Purpose: Compute Struve function H1(x)
!       Input :  x   --- Argument of H1(x) ( x ò 0 )
!       Output:  SH1 --- H1(x)
!       ============================================= }
Label   15, 25;
Var     A0,BY1,P1,Q1,R,S,T,T2,TA1: double;
	K, KM: integer;
Begin
        R:=1.0;
        if (X <= 20.0) then
        begin
           S:=0.0;
           A0:=-2.0/PI;
           for K:=1 to 60 do
           begin
              R:=-R*X*X/(4.0*K*K-1.0);
              S:=S+R;
              if (abs(R) < abs(S)*1.0E-12) then goto 15
           end;
15:        SH1:=A0*S
        end
        else
        begin
           S:=1.0;
           KM:=Round(0.5*X);
           if (X > 50.0)  then KM:=25;
           for K:=1 to KM do
           begin
              R:=-R*(4.0*K*K-1.0)/(X*X);
              S:=S+R;
              if (abs(R) < abs(S)*1.0E-12) then goto 25
           end;
25:        T:=4.0/X;
           T2:=T*T;
           P1:=((((0.42414e-5*T2-0.20092e-4)*T2+0.580759e-4)*T2-0.223203e-3)*T2+0.29218256e-2)*T2+0.3989422819;
           Q1:=T*(((((-0.36594e-5*T2+0.1622e-4)*T2-0.398708e-4)*T2+0.1064741e-3)*T2-0.63904e-3)*T2+0.0374008364);
           TA1:=X-0.75*PI;
           BY1:=2.0/sqrt(X)*(P1*sin(TA1)+Q1*cos(TA1));
           SH1:=2.0/PI*(1.0+S/(X*X))+BY1
        end
End;


{main program}
BEGIN

    writeln; 
    write(' Please enter x: ');
    readln(X);
    writeln;
    writeln('   x          H1(x)             ');
    writeln('--------------------------------');
        
    STVH1(X, SH1);
        
    writeln(X:5:1, '  ',SH1);

    Readkey;
    DoneWinCrt

END.

{end of file mstvh1.pas}