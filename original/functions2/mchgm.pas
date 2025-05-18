{******************************************************************
!*      Purpose: This program computes the confluent              *
!*               hypergeometric function M(a,b,x) using           *
!*               subroutine CHGM                                  *
!*      Input  : a  --- Parameter                                 *
!*               b  --- Parameter ( b <> 0,-1,-2,... )            *
!*               x  --- Argument                                  *
!*      Output:  HG --- M(a,b,x)                                  *
!*      Example:                                                  *
!*                  a       b       x          M(a,b,x)           *
!*                -----------------------------------------       *
!*                 1.5     2.0    20.0     .1208527185D+09        *
!*                 4.5     2.0    20.0     .1103561117D+12        *
!*                -1.5     2.0    20.0     .1004836854D+05        *
!*                -4.5     2.0    20.0    -.3936045244D+03        *
!*                 1.5     2.0    50.0     .8231906643D+21        *
!*                 4.5     2.0    50.0     .9310512715D+25        *
!*                -1.5     2.0    50.0     .2998660728D+16        *
!*                -4.5     2.0    50.0    -.1807604439D+13        *
!* -------------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special        *
!*             Functions jin.ece.uiuc.edu/routines/routines.html" *
!*                                                                *
!*                           Pascal Release By J-P Moreau, Paris. *
!*                                    (www.jpmoreau.fr)           *
!*****************************************************************}
Uses WinCrt;

Var  A, B, X, HG: double;


FUNCTION Power(Y:double; X:double): double;
BEGIN
  IF Y<=0.0 THEN
  begin
    writeln(' Error in Power.');
    Power:=0.0
  end
  ELSE
    Power:=Exp(x*Ln(y))
END;


Procedure GAMMA(X: double; Var GA: double);
{       ==================================================
!       Purpose: Compute gamma function â(x)
!       Input :  x  --- Argument of â(x)
!                       ( x is not equal to 0,-1,-2,úúú)
!       Output:  GA --- â(x)
!       ================================================== }
Var        G: Array[1..26] of double;
	   GR, R, Z: double;
	   K, M, M1: integer;
Begin
	   if X = Round(X) then
	     if X > 0.0 then
             begin
              GA:=1.0;
              M1:=Round(X-1);
              for K:=2 to M1 do GA := GA * K
             end
             else
              GA:=1.0E+300
	   else
           begin
	     if abs(X) > 1.0 then
             begin
               Z:=abs(X);
               M:=Round(Z);
               R:=1.0;
               for K:=1 to M do R:=R*(Z-K);
               Z:=Z-M
             end
             else
               Z:=X;

             G[1] := 1.0; G[2] := 0.5772156649015329;
	     G[3] := -0.6558780715202538; G[4] := -0.420026350340952e-1;
	     G[5] := 0.1665386113822915;  G[6] := -0.421977345555443e-1;
	     G[7] := -0.96219715278770e-2; G[8] := 0.72189432466630e-2;
	     G[9] := -0.11651675918591e-2; G[10] := -0.2152416741149e-3;
	     G[11] := 0.1280502823882e-3;  G[12] := -0.201348547807e-4;
	     G[13] := -0.12504934821e-5;   G[14] := 0.11330272320e-5;
	     G[15] := -0.2056338417e-6;    G[16] := 0.61160950e-8;
	     G[17] := 0.50020075e-8;       G[18] := -0.11812746e-8;
	     G[19] := 0.1043427e-9;        G[20] := 0.77823e-11;
	     G[21] := -0.36968e-11;        G[22] := 0.51e-12;
	     G[23] := -0.206e-13;          G[24] := -0.54e-14;
             G[25] := 14e-14;              G[26] := 1e-15;

	     GR:=G[26];
             for K:=25 Downto 1 do  GR:=GR*Z+G[K];
             GA:=1.0/(GR*Z);
             if abs(X) > 1.0 then
             begin
               GA:=GA*R;
               if X < 0.0 then GA:=-PI/(X*GA*sin(PI*X))
             end
           end
End;


Procedure CHGM(A, B, X: double; Var HG: double);
{       ==================================================
!       Purpose: Compute confluent hypergeometric function
!                M(a,b,x)
!       Input  : a  --- Parameter
!                b  --- Parameter ( b <> 0,-1,-2,... )
!                x  --- Argument
!       Output:  HG --- M(a,b,x)
!       Routine called: GAMMA for computing â(x)
!       ================================================== }
Label   25, fin;
Var     A0,A1,HG1,HG2,R,R1,R2,RG,SUM1,SUM2,TA,TB,TBA,X0,XG,Y0,Y1: double;
	I, J, K, LA, M, N, NL: integer;
Begin
        A0:=A;
        A1:=A;
        X0:=X;
        HG:=0.0;
        if (B = 0.0) or (B = -abs(Round(B))) then
           HG:=1.0E+300
        else if (A = 0.0) or (X = 0.0) then
           HG:=1.0
        else if A = -1.0 then
           HG:=1.0-X/B
        else if A = B then
           HG:=exp(X)
        else if A-B = 1.0 then
           HG:=(1.0+X/B)*exp(X)
        else if (A = 1.0) and (B = 2.0) then
           HG:=(exp(X)-1.0)/X
        else if (A = Int(A)) and (A < 0.0) then
        begin
           M:=Round(-A);
           R:=1.0;
           HG:=1.0;
           for K:=1 to M do
           begin
              R:=R*(A+K-1.0)/K/(B+K-1.0)*X;
              HG:=HG+R
           end
        end;
        if HG <> 0.0 then goto fin;
        if X < 0.0 then
        begin
           A:=B-A;
           A0:=A;
           X:=abs(X)
        end;
        if (A < 2.0) then NL:=0;
        if (A >= 2.0) then
        begin
           NL:=1;
           LA:=Round(A);
           A:=A-LA-1.0
        end;
        for N:=0 to NL do
        begin
           if (A0 >= 2.0) then A:=A+1.0;
           if (X <= 30.0+abs(B)) or (A < 0.0) then
           begin
              HG:=1.0;
              RG:=1.0;
              for J:=1 to 500 do
              begin
                 RG:=RG*(A+J-1.0)/(J*(B+J-1.0))*X;
                 HG:=HG+RG;
                 if (abs(RG/(HG)) < 1.0E-15) then goto 25;
              end
           end
           else
           begin
              GAMMA(A,TA);
              GAMMA(B,TB);
              XG:=B-A;
              GAMMA(XG,TBA);
              SUM1:=1.0;
              SUM2:=1.0;
              R1:=1.0;
              R2:=1.0;
              for I:=1 to 8 do
              begin
                 R1:=-R1*(A+I-1.0)*(A-B+I)/(X*I);
                 R2:=-R2*(B-A+I-1.0)*(A-I)/(X*I);
                 SUM1:=SUM1+R1;
                 SUM2:=SUM2+R2
              end;
              HG1:=TB/TBA*power(X,-A)*cos(PI*A)*SUM1;
              HG2:=TB/TA*exp(X)*power(X,A-B)*SUM2;
              HG:=HG1+HG2
           end;
25:        if (N = 0) then Y0:=HG;
           if (N = 1) then Y1:=HG;
	end;
        if A0 >= 2.0 then
	    for I:=1 to LA-1 do
            begin
              HG:=((2.0*A-B+X)*Y1+(B-A)*Y0)/A;
              Y0:=Y1;
              Y1:=HG;
              A:=A+1.0
            end;
        if (X0 < 0.0) then HG:=HG*exp(X0);
        A:=A1;
        X:=X0;
fin: End;


BEGIN

	writeln;
        writeln('    a     b     x          M(a,b,x)');
        writeln(' ------------------------------------------');

        A:=1.5; B:=2.0; X:=20.0;
	CHGM(A,B,X,HG);
	writeln(' ',A:5:1,' ',B:5:1,' ',X:5:1,' ',HG);

        A:=4.5;
	CHGM(A,B,X,HG);
	writeln(' ',A:5:1,' ',B:5:1,' ',X:5:1,' ',HG);

        A:=-1.5;
	CHGM(A,B,X,HG);
	writeln(' ',A:5:1,' ',B:5:1,' ',X:5:1,' ',HG);

        A:=-4.5;
	CHGM(A,B,X,HG);
	writeln(' ',A:5:1,' ',B:5:1,' ',X:5:1,' ',HG);

        A:=1.5; X:=50.0;
	CHGM(A,B,X,HG);
	writeln(' ',A:5:1,' ',B:5:1,' ',X:5:1,' ',HG);

        A:=4.5;
	CHGM(A,B,X,HG);
	writeln(' ',A:5:1,' ',B:5:1,' ',X:5:1,' ',HG);

        A:=-1.5;
	CHGM(A,B,X,HG);
	writeln(' ',A:5:1,' ',B:5:1,' ',X:5:1,' ',HG);

        A:=-4.5;
	CHGM(A,B,X,HG);
	writeln(' ',A:5:1,' ',B:5:1,' ',X:5:1,' ',HG);

        Readln;

END.

{end of file mchgm.pas}