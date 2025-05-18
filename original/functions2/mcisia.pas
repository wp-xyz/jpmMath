{******************************************************************
!*      Purpose: This program computes the cosine and sine        *
!*               integrals using subroutine CISIA                 *
!*      Input :  x  --- Argument of Ci(x) and Si(x)               *
!*      Output:  CI --- Ci(x)                                     *
!*               SI --- Si(x)                                     *
!*      Example:                                                  * 
!*                   x         Ci(x)          Si(x)               *
!*                ------------------------------------            *
!*                  0.0     - ì             .00000000             *
!*                  5.0     -.19002975     1.54993124             *
!*                 10.0     -.04545643     1.65834759             *
!*                 20.0      .04441982     1.54824170             * 
!*                 30.0     -.03303242     1.56675654             *
!*                 40.0      .01902001     1.58698512             *
!* -------------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special        *
!*             Functions jin.ece.uiuc.edu/routines/routines.html" *
!*                                                                *
!*                           Pascal Release By J-P Moreau, Paris. *
!*                                     (www.jpmoreau.fr)          *
!*****************************************************************}
Program MCISIA;

Uses WinCrt1;

Var
      CI,SI,X: Double;

{calculate y power n}
FUNCTION Power(y:Double; n:Integer): Double;
var i: integer; result: double;
begin
  result :=1.0;
  if n=0 then begin
    Power:=result;
    exit;
  end
  else for i:=1 to n do result := y * result;
  Power:=result;
end;


Procedure CISIA(X: double; Var CI, SI: Double);

{       =============================================
!       Purpose: Compute cosine and sine integrals
!                Si(x) and Ci(x)  ( x ò 0 )
!       Input :  x  --- Argument of Ci(x) and Si(x)
!       Output:  CI --- Ci(x)
!                SI --- Si(x)
!       ============================================= }
Label   10, 20, fin;
Var     EL, EPS, P2, X2, XA, XA0, XA1, XCS, XF, XG, XG1, XG2, XM, XR, XS, XSS: Double;
        BJ: Array[1..100] of Double;
	K, M: Integer;
Begin
        P2:=1.570796326794897;
        EL:=0.5772156649015329;
        EPS:=1.0E-15;
        X2:=X*X;
        if (X = 0.0) then
        begin
           CI:=-1.0E+100;
           SI:=0.0
        end
        else if (X <= 16.0) then
        begin
           XR:=-0.25*X2;
           CI:=EL+Ln(X)+XR;
           for K:=2 to 40 do
           begin
              XR:=-0.5*XR*(K-1)/(K*K*(2*K-1))*X2;
              CI:=CI+XR;
              if (abs(XR) < abs(CI)*EPS) then goto 10
	   end;
10:        XR:=X;
           SI:=X;
           for K:=1 to 40 do
           begin
              XR:=-0.5*XR*(2*K-1)/K/(4*K*K+4*K+1)*X2;
              SI:=SI+XR;
              if (abs(XR) < abs(SI)*EPS) then goto fin;
	   end
        end
        else if X <= 32.0 then
        begin
           XM:=Int(47.2+0.82*X);
           M:=Round(XM);
           XA1:=0.0;
           XA0:=1.0E-100;
	   for K:=M Downto 1 do
           begin
              XA:=4.0*K*XA0/X-XA1;
              BJ[K]:=XA;
              XA1:=XA0;
              XA0:=XA
           end;
           XS:=BJ[1];

           K:=3;
20:        XS:=XS+2.0*BJ[K];
           Inc(K,2);
           if K<=M then goto 20;

           BJ[1]:=BJ[1]/XS;
           for K:=2 to M do   BJ[K]:=BJ[K]/XS;
           XR:=1.0;
           XG1:=BJ[1];
	   for K:=2 to M do
           begin
              XR:=0.25*XR*power(2.0*K-3.0,2)/((K-1.0)*power(2.0*K-1.0,2))*X;
              XG1:=XG1+BJ[K]*XR
           end;
           XR:=1.0;
           XG2:=BJ[1];
           for K:=2 to M do
           begin
              XR:=0.25*XR*power(2.0*K-5.0,2)/((K-1.0)*power(2.0*K-3.0,2))*X;
              XG2:=XG2+BJ[K]*XR
	   end;
           XCS:=cos(X/2.0);
           XSS:=sin(X/2.0);
           CI:=EL+Ln(X)-X*XSS*XG1+2*XCS*XG2-2*XCS*XCS;
           SI:=X*XCS*XG1+2*XSS*XG2-sin(X);
        end
        else
        begin
           XR:=1.0;
           XF:=1.0;
	   for K:=1 to 9 do
           begin
              XR:=-2.0*XR*K*(2*K-1)/X2;
              XF:=XF+XR
           end;
           XR:=1.0/X;
           XG:=XR;
           for K:=1 to 8 do
           begin
              XR:=-2.0*XR*(2*K+1)*K/X2;
              XG:=XG+XR
           end;
           CI:=XF*sin(X)/X-XG*cos(X)/X;
           SI:=P2-XF*cos(X)/X-XG*sin(X)/X
        end;
fin: End;


{main program}
BEGIN
        writeln;
        writeln('   x               Ci(x)                      Si(x)            ');
        writeln(' --------------------------------------------------------------');
       
	X:=0.0;
	CISIA(X, CI, SI);
        writeln('   0.0     -i                           0');

	X:=5.0;
        CISIA(X, CI, SI);
        writeln(' ',X:5:1,'     ',CI,'     ',SI);

        X:=10.0;
	Repeat
          CISIA(X, CI, SI);
          writeln(' ',X:5:1,'     ',CI,'     ',SI);
          X:=X+10.0
        Until X > 40.0;
	writeln;
        ReadKey;
        DoneWinCrt
END.

{end of file mcisia.pas}