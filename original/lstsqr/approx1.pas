{*******************************************************
*  LEAST SQUARE APPROXIMATION OF A DISCRETE FUNCTION   *
*  USING ORTHOGONAL POLYNOMIALS                        *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
* (Approximate SIN(X) defined by 60 points, from x=0   *
*  to x=PI/2, with degree K=7).                        *
*                                                      *
*  I     COEFFICIENTS   STD DEVIATION                  *
*  0     0.63383395       0.14142136                   *
*  1     0.66275983       0.30570229                   *
*  2    -0.37764682       0.73926175                   *
*  3    -0.11371820       1.82234391                   * 
*  4     0.02861664       4.52658861                   *
*  5     0.00574936      11.29632598                   *
*  6    -0.00096142      28.29696331                   *
*  7    -0.00013770      71.13668785                   *
*                                                      *  
*  I     VARIABLE        EXACT R         APPROX. R     *
*  1     0.00000000      0.00000000     -0.00000003    *
*  2     0.31415927      0.30901699      0.30901698    *
*  3     0.62831853      0.58778525      0.58778525    *
*  4     0.94247780      0.80901699      0.80901700    *
*  5     1.25663706      0.95105652      0.95105651    *
*  6     1.57079633      1.00000000      0.99999997    *
*                                                      *
* ---------------------------------------------------- *    
* Ref.: From Numath Library By Tuan Dang Trong in      *
*       Fortran 77 [BIBLI 18].                         *
*                                                      *
*                  TPW Release By J-P Moreau, Paris.   *
*                         (www.jpmoreau.fr)            *
*******************************************************}
PROGRAM APPROX;

Uses  WinCrt;

Const NMAX = 100;

Type  pVEC = ^VEC;
      VEC = Array[0..NMAX] of Double;

Var   I,K,M,NW: Integer;
      X,Y,SIGMA,S,ECART,ALPHA,BETA: pVEC;
      W1,W2,W3,W4: pVEC;
      R,W,W0,DW: Double;


    {Auxiliary function}
    FUNCTION PRD(X,Y,Z:pVEC;M:Integer): Double;
    Var SUM: Double;
        I: Integer;
    Begin
      SUM:=0.0;
      For I:=1 to M do
        SUM:=SUM+X^[I]*Y^[I]/Sqr(Z^[I]);
      PRD:=SUM
    End;


    Procedure MCARRE(K,M:Integer; X,Y: pVEC; Var SIGMA,S,ALPHA,BETA,ECART,P1,P2,P3,P4:pVEC);
{-----------------------------------------------------------------------
!     LEAST SQUARES APPROXIMATION OF A FUNCTION F(X) DEFINED BY M POINTS
!     X(I), Y(I) BY USING ORTHOGONAL POLYNOMIALS
!-----------------------------------------------------------------------
!     INPUTS:
!     K   : DEGREE OF POLYNOMIALS
!     M   : NUMBER OF POINTS
!     X,Y : TABLES OF DIMENSION M TO STORE M ABSCISSAS AND
!           M ORDINATES OF GIVEN POINTS
!     SIGMA : TABLE OF DIMENSION M TO STORE THE STANDARD DEVIATIONS
!             OF VARIABLE Y
!     OUTPUTS:
!     S    : TABLE OF DIMENSION(0:K)
!     ALPHA: TABLE OF DIMENSION (K)
!     BETA : TABLE OF DIMENSION (0:K-1)
!     WORKING SPACE:
!     W1..W4: AUXILIARY TABLES OF DIMENSION (M)
!     NOTE:
!     COEFFICIENTS S,ALPHA,BETA ARE USED TO EVALUATE VALUE
!     AT POINT Z OF THE BEST POLYNOMIAL OF DEGREE K
!     BY USING FUNCTION P(K,S,ALPHA,BETA,Z) DESCRIBED BELOW.
!----------------------------------------------------------------------}
    Var I,L: Integer;
        OMEGA,T,W,WPR: Double;
    Begin
      For I:=1 to M do
      begin
        P1^[I]:=0.0;
        P2^[I]:=1.0
      end;
      W:=1.0*M;
      BETA^[0]:=0.0;
      For I:=0 to K-1 do
      begin
        OMEGA:=PRD(Y,P2,SIGMA,M);
        S^[I]:=OMEGA/W;
        T:=PRD(P2,P2,SIGMA,M)/Sqr(W);
        ECART^[I]:=SQRT(T);
        For L:=1 to M  do P4^[L]:=X^[L]*P2^[L];
        ALPHA^[I+1]:=PRD(P4,P2,SIGMA,M)/W;
        For L:=1 to M do
          P3^[L]:=(X^[L]-ALPHA^[I+1])*P2^[L]-BETA^[I]*P1^[L];
        WPR:=PRD(P3,P3,SIGMA,M);
        If I+1 <= K-1 then  BETA^[I+1]:=WPR/W;
        W:=WPR;
        For L:=1 to M do
        begin
          P1^[L]:=P2^[L];
          P2^[L]:=P3^[L]
        end
      end;
      OMEGA:=PRD(Y,P2,SIGMA,M);
      S^[K]:=OMEGA/W;
      T:=PRD(P2,P2,SIGMA,M)/Sqr(W);
      ECART^[K]:=SQRT(T)
    End;


    FUNCTION P(K:Integer; S,ALPHA,BETA:pVEC;X:Double): Double;
{---------------------------------------------------------------------
!     THIS FUNCTIONE ALLOWS EVALUATING VALUE AT POINT X OF A FUNCTION
!     F(X), APPROXIMATED BY A SYSTEM OF ORTHOGONAL POLYNOMIALS Pj(X)
!     THE COEFFICIENTS OF WHICH, ALPHA,BETA HAVE BEEN DETERMINED BY
!     LEAST SQUARES.
!--------------------------------------------------------------------}
    Var B,BPR,BSD: Double;
    Var I: Integer;
    Begin
      B:=S^[K];
      BPR:=S^[K-1]+(X-ALPHA^[K])*S^[K];
      For I:=K-2 Downto 0 do
      begin
        BSD:=S^[I]+(X-ALPHA^[I+1])*BPR-BETA^[I+1]*B;
        B:=BPR;
        BPR:=BSD
      end;
      P:=BPR
    End;


{main program}
BEGIN

      K:=7;
      M:=50;
      NW:=6;

      New(X); New(Y); New(SIGMA); New(S);
      New(ECART); New(ALPHA); New(BETA);
      New(W1); New(W2); New(W3); New(W4);

      For I:=1 to M do
      begin
        X^[I]:=PI*(I-1)/(2*(M-1));
        Y^[I]:=SIN(X^[I]);
        SIGMA^[I]:=1.0
      end;

      MCARRE(K,M,X,Y,SIGMA,S,ALPHA,BETA,ECART,W1,W2,W3,W4);

      writeln;
      WRITELN('   I    COEFFICIENTS    STD DEVIATION');
      For I:= 0 to K do
        writeln(I:4,'  ',S^[I]:12:8,'      ',ECART^[I]:12:8);
      writeln;
      WRITELN('   I    VARIABLE          EXACT R         APPROX. R');
      W0:=0.0;
      DW:=PI/10.0;
      For I:=1 to NW do
      begin
        W:=W0+(I-1)*DW;
        R:=P(K,S,ALPHA,BETA,W);
        writeln(I:4,'  ',W:12:8,'      ',SIN(W):12:8,'    ',R:12:8);
      end;
      ReadKey;
      Dispose(X); Dispose(Y); Dispose(SIGMA); Dispose(S);
      Dispose(ECART); Dispose(ALPHA); Dispose(BETA);
      Dispose(W1); Dispose(W2); Dispose(W3); Dispose(W4);
      DoneWinCrt

END. {OF MAIN PROGRAM}


{end of file approx.pas}