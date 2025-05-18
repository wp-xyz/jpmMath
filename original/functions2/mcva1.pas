{********************************************************************
!*      Purpose: This program computes a sequence of characteristic *
!*               values of Mathieu functions using subroutine CVA1  *
!*      Input :  m  --- Order of Mathieu functions                  *
!*               q  --- Parameter of Mathieu functions              *
!*               KD --- Case code                                   *
!*                      KD=1 for cem(x,q)  ( m = 0,2,4,...)         *
!*                      KD=2 for cem(x,q)  ( m = 1,3,5,...)         *
!*                      KD=3 for sem(x,q)  ( m = 1,3,5,...)         *
!*                      KD=4 for sem(x,q)  ( m = 2,4,6,...)         *
!*      Output:  CV^[I] --- Characteristic values; I = 1,2,3,...    *
!*               For KD=1, CV^[1], CV^[2], CV(3),..., correspond to *
!*               For KD=2, CV^[1], CV^[2], CV(3),..., correspond to *
!*               the characteristic values of cem for m = 1,3,5,..  *
!*               For KD=3, CV^[1], CV^[2], CV(3),..., correspond to *
!*               the characteristic values of sem for m = 1,3,5,..  *
!*               For KD=4, CV^[1], CV^[2], CV(3),..., correspond to *
!*               the characteristic values of sem for m = 0,2,4,..  *
!*                                                                  *
!*      Example: Mmax = 12,    q = 25.00                            *
!*                                                                  *
!*               Characteristic values of Mathieu functions         *
!*                                                                  *
!*                 m            a                  b                *
!*               ------------------------------------------         *
!*                 0      -40.256779547                             *
!*                 1      -21.314899691      -40.256778985          *
!*                 2       -3.522164727      -21.314860622          *
!*                 3       12.964079444       -3.520941527          *
!*                 4       27.805240581       12.986489953          *
!*                 5       40.050190986       28.062765899          *
!*                 6       48.975786716       41.801071292          *
!*                 7       57.534689001       55.002957151          *
!*                 8       69.524065166       69.057988351          *
!*                 9       85.076999882       85.023356505          *
!*                10      103.230204804      103.225680042          *
!*                11      123.643012376      123.642713667          *
!*                12      146.207690643      146.207674647          *
!* ---------------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special          *
!*             Functions jin.ece.uiuc.edu/routines/routines.html"   *
!*                                                                  *
!*                            Pascal Release By J-P Moreau, Paris.  *
!*                                     (www.jpmoreau.fr)            *
!*******************************************************************}
        PROGRAM MCVA1;
 Type   pVEC = ^VEC;
        VEC = Array[1..200] of Double;
        pVEC1 = ^VEC1;
        VEC1 = Array[1..500] of Double;

Var		
        CV1, CV2, CVE, CVS: pVec;
        J, KD, MMAX: Integer;
       	Q: Double;

{calculate y power n}
FUNCTION Power(y:Double; n:Integer): Double;
var i: integer; result: double;
begin
  result :=1.0;
  if n=0 then begin
    Power:=result;
    exit;
  end
  else for i:=1 to abs(n) do result := y * result;
  if n>0 then Power:=result else Power:=1.0/result
end;

Procedure CVA1(KD,M: Integer; Q: Double; CV:pVEC);
{       ============================================================
!       Purpose: Compute a sequence of characteristic values of
!                Mathieu functions
!       Input :  M  --- Maximum order of Mathieu functions
!                q  --- Parameter of Mathieu functions
!                KD --- Case code
!                       KD=1 for cem(x,q)  ( m = 0,2,4,... )
!                       KD=2 for cem(x,q)  ( m = 1,3,5,... )
!                       KD=3 for sem(x,q)  ( m = 1,3,5,... )
!                       KD=4 for sem(x,q)  ( m = 2,4,6,... )
!       Output:  CV^[I] --- Characteristic values; I = 1,2,3,...
!                For KD=1, CV^[1], CV^[2], CV(3),..., correspond to
!                the characteristic values of cem for m = 0,2,4,...
!                For KD=2, CV^[1], CV^[2], CV(3),..., correspond to
!                the characteristic values of cem for m = 1,3,5,...
!                For KD=3, CV^[1], CV^[2], CV(3),..., correspond to
!                the characteristic values of sem for m = 1,3,5,...
!                For KD=4, CV^[1], CV^[2], CV(3),..., correspond to
!                the characteristic values of sem for m = 0,2,4,...
!       ============================================================ }
Label   55, 60, 70;
Var     G, H: pVEC;
        D, E, F: pVEC1;
		EPS, S,T,T1, X1,XA,XB: Double;
		I,IC,ICM, J,K,K1, NM,NM1: Integer;
Begin
        New(G); New(H); New(D); New(E); New(F);
        EPS:=1E-14;
        ICM:=Round(M/2)+1;
        IF (KD = 4) then ICM:=(M Div 2);
        IF (Q = 0.0) THEN
           IF (KD = 1) THEN
              For IC:=1 to ICM do
                CV^[IC]:=4.0*(IC-1.0)*(IC-1.0)
           ELSE IF (KD <> 4) THEN
              For IC:=1 to ICM do
                CV^[IC]:=(2.0*IC-1.0)*(2.0*IC-1.0)
           ELSE
              For IC:=1 to ICM do
                CV^[IC]:=4.0*IC*IC
        ELSE
	begin
           NM:=Round(10+1.5*M+0.5*Q);
           E^[1]:=0.0;
           F^[1]:=0.0;
           IF (KD = 1) THEN
           begin
              D^[1]:=0.0;
              For I:=2 to NM do
			  begin
                 D^[I]:=4.0*(I-1)*(I-1);
                 E^[I]:=Q;
                 F^[I]:=Q*Q
              end;
              E^[2]:=SQRT(2.0)*Q;
              F^[2]:=2.0*Q*Q
           end
           ELSE IF KD <> 4 THEN
           begin
              D^[1]:=1.0+Power(-1,KD)*Q;
              For I:=2 to NM do
			  begin
                 D^[I]:=(2.0*I-1.0)*(2.0*I-1.0);
                 E^[I]:=Q;
                 F^[I]:=Q*Q
              end
           end
           ELSE
	       begin
              D^[1]:=4.0;
              For I:=2 to NM do
	          begin
                 D^[I]:=4.0*I*I;
                 E^[I]:=Q;
                 F^[I]:=Q*Q
              end
           end;

           XA:=D^[NM]+ABS(E^[NM]);
           XB:=D^[NM]-ABS(E^[NM]);
           NM1:=NM-1;
           For I:=1 to NM1 do
	       begin
              T:=ABS(E^[I])+ABS(E^[I+1]);
              T1:=D^[I]+T;
              IF (XA < T1) then XA:=T1;
              T1:=D^[I]-T;
              IF (T1 < XB) then XB:=T1
           end;
           For I:=1 to ICM do
	   begin
              G^[I]:=XA;
              H^[I]:=XB
           end;
           For K:=1 to ICM do
	   begin
              For K1:=K to ICM do
	      begin
                 IF (G^[K1] < G^[K]) THEN
	      	 begin
                    G^[K]:=G^[K1];
                    GOTO 55
                 end
              end;
55:           IF (K <> 1) AND (H^[K] < H^[K-1]) then H^[K]:=H^[K-1];
60:           X1:=(G^[K]+H^[K])/2.0;
              CV^[K]:=X1;
              IF (ABS((G^[K]-H^[K])/X1) < EPS) then GOTO 70;
              J:=0;
              S:=1.0;
              For I:=1 to NM do
	      begin
                 IF (S = 0.0) then S:=S+1E-30;
                 T:=F^[I]/S;
                 S:=D^[I]-T-X1;
                 IF (S < 0.0) then J:=J+1
              end;
              IF (J < K) THEN
                 H^[K]:=X1
              ELSE
	      begin
                 G^[K]:=X1;
                 IF (J >= ICM) THEN
                    G^[ICM]:=X1
                 ELSE
	      	 begin
                    IF (H^[J+1] < X1) then H^[J+1]:=X1;
                    IF (X1 < G^[J]) then G^[J]:=X1
                 end
              end;
              GOTO 60;
70:           CV^[K]:=X1
           end
        end;
        Dispose(G); Dispose(H); Dispose(D); Dispose(E); Dispose(F)
End;


{main program}
BEGIN
        New(CV1); New(CV2); New(CVE); New(CVS);
	Writeln;
        WRITE(' Please enter Mmax,q: ');
        READLN(MMAX, Q);
        Writeln;
        CVA1(1,MMAX,Q,CV1);
        CVA1(2,MMAX,Q,CV2);
        For J:=1 to (MMAX Div 2)+1 do
	begin
           CVE^[2*J-1]:=CV1^[J];
           CVE^[2*J]:=CV2^[J]
        end;
        CVA1(3,MMAX,Q,CV1);
        CVA1(4,MMAX,Q,CV2);
        For J:=1 to (MMAX Div 2)+1 do
	begin
           CVS^[2*J]:=CV1^[J];
           CVS^[2*J+1]:=CV2^[J]
        end;
        WRITELn(' Characteristic values of Mathieu functions');
        WRITELn;
        WRITELn('  m            a                   b');
        WRITELn('-------------------------------------------------');
        For J:=0 to MMAX do
	begin
           IF (J = 0) then WRITELn(' ',J,' ',CVE^[J+1]);
           IF (J <> 0) then WRITELn(' ',J,' ',CVE^[J+1],' ',CVS^[J+1])
        end;
	Readln;
END.

{end of file mcva1.pas}