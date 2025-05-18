{**************************************************************
!*       Purpose: This program computes Struve function       * 
!*                H0(x) using subroutine STVH0                *
!*       Input :  x   --- Argument of H0(x) ( x ò 0 )         *
!*       Output:  SH0 --- H0(x)                               *
!*       Example:                                             *
!*                   x          H0(x)                         *
!*                ----------------------                      *
!*                  0.0       .00000000                       * 
!*                  5.0      -.18521682                       *
!*                 10.0       .11874368                       *
!*                 15.0       .24772383                       *
!*                 20.0       .09439370                       *
!*                 25.0      -.10182519                       *
!* ---------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special    *
!*             Functions, jin.ece.uiuc.edu/routines/routines  *
!*             .html".                                        *
!*                                                            *
!*                       Pascal Release By J-P Moreau, Paris. *
!*                                 (www.jpmoreau.fr)          *
!*************************************************************} 
Uses WinCrt;

Var  X,SH0: Double;


Procedure STVH0(X: double; Var SH0: double);
{       =============================================
!       Purpose: Compute Struve function H0(x)
!       Input :  x   --- Argument of H0(x) ( x ò 0 )
!       Output:  SH0 --- H0(x)
!       ============================================= }
Label   15, 25;
Var     A0,BY0,P0,Q0,R,S,T,T2,TA0: double;
		K, KM: integer;
Begin	    
        S:=1.0;
        R:=1.0;
        if (X <:= 20.0) then
		begin
           A0:=2.0*X/PI;
           for K:=1 to 60 do
		   begin
              R:=-R*X/(2.0*K+1.0)*X/(2.0*K+1.0);
              S:=S+R;
              if (abs(R) < abs(S)*1.0e-12) then goto 15;
           end;
15:        SH0:=A0*S;
        end
        else
		begin
           KM:=Round(0.5*(X+1.0));
           if (X >= 50.0) then KM:=25;
           for K:=1 to KM do
		   begin
              R:=-R*Sqr((2.0*K-1.0)/X);
              S:=S+R;
              if (abs(R) < abs(S)*1.0e-12) then goto 25;
           end;
25:        T:=4.0/X;
           T2:=T*T;
           P0:=((((-.37043e-5*T2+.173565e-4)*T2-.487613e-4)*T2+.17343e-3)*T2-0.1753062e-2)*T2+.3989422793;
           Q0:=T*(((((.32312e-5*T2-0.142078e-4)*T2+0.342468e-4)*T2-0.869791e-4)*T2+0.4564324e-3)*T2-0.0124669441);
           TA0:=X-0.25*PI;
           BY0:=2.0/sqrt(X)*(P0*sin(TA0)+Q0*cos(TA0));
           SH0:=2.0/(PI*X)*S+BY0;
        end
End;


BEGIN
        
	writeln;
	write(' Please enter x: ');
        readln(X);
	writeln;
        writeln('   x        H0(x)     ');
        writeln('----------------------');
        
	STVH0(X, SH0);
        
	writeln(X:5:1,'  ',SH0);

	Readln;
	DoneWinCrt

END.

{end of file mstvh0.pas}