{*****************************************************************************
*              RESOLUTION OF A SET OF NON-LINEAR EQUATIONS                   *
* -------------------------------------------------------------------------- *
* COPYRIGHT 1991, BY R.S. BAIN                                               *
* Originally wrapped by bain@luther.che.wisc.edu on Tue Mar 30 16:27:48 1993 *
*                                                                            *
* As of this writing, we do not know how to contact Rod Bain, the            * 
* author of NNES.                                                            *
*                                                                            *
* -- David M. Gay (dmg@bell-labs.com)                                        *
*    Bell Labs, Murray Hill                                                  *
*    8 February 1999                                                         *
*                                                                            *
*                                   Pascal Release 1.1 By J-P Moreau, Paris. *
*                                              (PART 3/3)                    *
*                                           (www.jpmoreau.fr)                *
*                                                                            *
* Release 1.1: Corrected bug in qrdcom                                       * 
*****************************************************************************}
Unit Unnes2;

Interface

Uses WinCrt, Utils, Fcn, Strings, Unnes, Unnes1;
 
  Procedure nnes (var absnew, cauchy, deuflh, geoms, linesr, newton, overch:boolean;
                  acptcr:integer; var itsclf, itsclx:integer; jactyp, jupdm:integer;
                  var maxexp:integer; maxit, maxns, maxqns, mgll, minqns, n, narmij,
                  niejev, njacch:integer; var njetot:integer; nunit:integer; var output:
                  integer; qnupdm:integer; var stopcr:integer; supprs:integer; var trmcod,
                  trupdm:integer; var alpha, confac, delta, delfac, epsmch, etafac, fcnnew,
                  fdtolj, ftol, lam0, mstpf, nsttol, omega, ratiof, sigma, stptol:double;
                  a:pMat; boundl, boundu, delf, fsave, ftrack, fvec, fvecc:pVec; h:pMat;
                  hhpi:pVec; jac, plee:pMat; rdiag, s, sbar, scalef, scalex, sn, ssdhat,
                  strack, vhat, xc, xplus, xsave: pVec; help: String);


Implementation

Procedure qform (n:integer; a:pMat; hhpi:pVec; jac:pMat);
{-----------------------------------------------------
!    FEB. 14, 1991
!
!    FORM Q^  FROM THE HOUSEHOLDER MATRICES STORED IN
!    MATRICES A AND HHPI AND STORE IT IN JAC.
!----------------------------------------------------}
Var
    tau: double;
    i,j,k: integer;
    V1,V2:pVec;
Begin
  New(V1); New(V2);
  For i:=1 to n do
    For j:=1 to n do
      jac^[i,j] := zero;
  For j:=1 to n do jac^[j,j] := one;
  For k:=1 to n-1 do
    IF hhpi^[k]<> zero THEN
      For j:=1 to n do
      begin
        For i:=k to n do
        begin
          V1^[i]:=a^[i,k]; V2^[i]:=jac^[i,j]
        end;
        tau := DOT_PRODUCT(n, V1, V2);
        tau := tau/hhpi^[k];
        For i:=k to n do 
          jac^[i,j] := jac^[i,j] - tau*a^[i,k]
      end;
  Dispose(V1); Dispose(V2)
End; {qform}

Procedure qrdcom (var qrsing:boolean; n:integer; epsmch:double;
                  a:pMat; hhpi, rdiag:pVec);
{---------------------------------------------------------------------
!    FEB. 23, 1992
!
!    THIS SUBROUTINE COMPUTES THE QR DECOMPOSITION OF THE MATRIX A.
!    THE DECOMPOSITION IS COMPLETED EVEN IF A SINGULARITY IS DETECTED
!    (WHEREUPON QRSING IS SET TO TRUE).
!--------------------------------------------------------------------}
Var
    eta, sigma, tau: double;
    i, j, k: integer;
    V1,V2:pVec;
Begin
  New(V1); New(V2);
  qrsing:=false;
  For k:=1 to n-1 do
  begin
    eta:=zero;
    For i:=k to n do eta:=MAX(eta, ABS(a^[i,k]));
    IF eta < epsmch THEN
    begin
      qrsing:=true;
      hhpi^[k]:=zero;
      rdiag^[k]:=zero
    end
    ELSE
    begin
      For j:= k to n do a^[j,k]:=a^[j,k]/eta;
      sigma:=zero;
      For j:= k to n do sigma := sigma + Sqr(a^[j,k]);
      sigma:=SIGN(SQRT(sigma),a^[k,k]);
      a^[k,k]:=a^[k,k] + sigma;
      hhpi^[k]:=sigma*a^[k,k];
      rdiag^[k]:=-eta*sigma;
      For j:=k+1 to n do
      begin
        For i:=k to n do
        begin
          V1^[i-k+1]:=a^[i,k]; V2^[i-k+1]:=a^[i,j];  {corected bug 04/05/07}
        end;
        tau := DOT_PRODUCT(n-k+1,V1,V2);
        tau:=tau/hhpi^[k];
        For i:=k to n do a^[i,j] := a^[i,j] - tau*a^[i,k]
      end
    end
  end;
  rdiag^[n]:=a^[n,n];
  IF ABS(rdiag^[n]) < epsmch then qrsing:=true;
  Dispose(V1); Dispose(V2)
End; {qrdcom}

Procedure utumul (nradec, ncadec, nraact, ncaact, nrbdec, ncbdec:integer;
                  amat, bmat: pMat);
{-----------------------------------------------------------------------------
!   FEB. 8, 1991
!
!   MATRIX MULTIPLICATION:   A^A:=B   WHERE A IS UPPER TRIANGULAR
!
!   VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4.
!
!   NRADEC IS NUMBER OF ROWS IN A DECLARED
!   NCADEC IS NUMBER OF COLUMNS IN A DECLARED
!   NRAACT IS THE LIMIT FOR THE 1ST INDEX IN A
!   NCAACT IS THE LIMIT FOR THE 2ND INDEX IN A
!   NRBDEC IS NUMBER OF ROWS IN B DECLARED
!   NCBDEC IS NUMBER OF COLUMNS IN B DECLARED
!
!   MODIFIED VERSION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES MACKINNON,
!   QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!-----------------------------------------------------------------------------}
Var 
    i, j, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, ncc32r,
    nend: integer;  SUM: double;
Begin

{  FIND ENTRY IN MATRIX B }

For i:=1 to ncaact do
begin

{  FIND NUMBER OF GROUPS OF SIZE 32, 16...}

  nend:=IMIN(i,nraact);

  ncc32:=nend Div 32;
  ncc32r:=nend - 32*ncc32;
  ncc16:=ncc32r Div 16;
  ncc16r:=ncc32r - 16*ncc16;
  ncc8:=ncc16r Div 8;
  ncc8r:=ncc16r - 8*ncc8;
  ncc4:=ncc8r Div 4;
  ncc4r:=ncc8r - 4*ncc4;

  For j:=i to ncaact do
  begin
    sum:=zero;
    k:=0;
    IF ncc32 > 0 THEN
      For kk:=1 to ncc32 do
      begin
        k:=k+32;
        sum:=sum + amat^[k-31,i]*amat^[k-31,j]+amat^[k-30,i]*amat^[k-30,
            j]+amat^[k-29,i]*amat^[k-29,j]+amat^[k-28,i]*amat^[k-28,j]+
            amat^[k-27,i]*amat^[k-27,j]+amat^[k-26,i]*amat^[k-26,j]+
            amat^[k-25,i]*amat^[k-25,j]+amat^[k-24,i]*amat^[k-24,j];
        sum:=sum + amat^[k-23,i]*amat^[k-23,j]+amat^[k-22,i]*amat^[k-22,
            j]+amat^[k-21,i]*amat^[k-21,j]+amat^[k-20,i]*amat^[k-20,j]+
            amat^[k-19,i]*amat^[k-19,j]+amat^[k-18,i]*amat^[k-18,j]+
            amat^[k-17,i]*amat^[k-17,j]+amat^[k-16,i]*amat^[k-16,j];
        sum:=sum + amat^[k-15,i]*amat^[k-15,j]+amat^[k-14,i]*amat^[k-14,
            j]+amat^[k-13,i]*amat^[k-13,j]+amat^[k-12,i]*amat^[k-12,j]+
            amat^[k-11,i]*amat^[k-11,j]+amat^[k-10,i]*amat^[k-10,j]+
            amat^[k-9,i]*amat^[k-9,j]+amat^[k-8,i]*amat^[k-8,j];
        sum:=sum + amat^[k-7,i]*amat^[k-7,j]+amat^[k-6,i]*amat^[k-6,j]+
            amat^[k-5,i]*amat^[k-5,j]+amat^[k-4,i]*amat^[k-4,j]+amat^[k-3,
            i]*amat^[k-3,j]+amat^[k-2,i]*amat^[k-2,j]+amat^[k-1,i]*
            amat^[k-1,j]+amat^[k,i]*amat^[k,j]
      end;

    IF ncc16 > 0 THEN
      For kk:=1 to ncc16 do
      begin
        k:=k+16;
        sum:=sum + amat^[k-15,i]*amat^[k-15,j]+amat^[k-14,i]*amat^[k-14,
            j]+amat^[k-13,i]*amat^[k-13,j]+amat^[k-12,i]*amat^[k-12,j]+
            amat^[k-11,i]*amat^[k-11,j]+amat^[k-10,i]*amat^[k-10,j]+
            amat^[k-9,i]*amat^[k-9,j]+amat^[k-8,i]*amat^[k-8,j];
        sum:=sum + amat^[k-7,i]*amat^[k-7,j]+amat^[k-6,i]*amat^[k-6,j]+
            amat^[k-5,i]*amat^[k-5,j]+amat^[k-4,i]*amat^[k-4,j]+amat^[k-3,
            i]*amat^[k-3,j]+amat^[k-2,i]*amat^[k-2,j]+amat^[k-1,i]*
            amat^[k-1,j]+amat^[k,i]*amat^[k,j]
      end;

    IF ncc8 > 0 THEN
      For kk:=1 to ncc8 do
      begin
        k:=k+8;
        sum:=sum + amat^[k-7,i]*amat^[k-7,j]+amat^[k-6,i]*amat^[k-6,j]+
            amat^[k-5,i]*amat^[k-5,j]+amat^[k-4,i]*amat^[k-4,j]+amat^[k-3,
            i]*amat^[k-3,j]+amat^[k-2,i]*amat^[k-2,j]+amat^[k-1,i]*
            amat^[k-1,j]+amat^[k,i]*amat^[k,j]
      end;

    IF ncc4 > 0 THEN
      For kk:=1 to ncc4 do
      begin
        k:=k+4;
        sum:=sum + amat^[k-3,i]*amat^[k-3,j]+amat^[k-2,i]*amat^[k-2,j] +
            amat^[k-1,i]*amat^[k-1,j]+amat^[k,i]*amat^[k,j]
      end;

    IF ncc4r > 0 THEN
      For kk:=1 to ncc4r do
      begin
        k:=k+1;
        sum:=sum + amat^[k,i]*amat^[k,j]
      end;

    bmat^[i,j]:=sum;
    IF i <> j then  bmat^[j,i]:=bmat^[i,j]
  end
end

End; {utumul}

Procedure rtrmul (n:integer; a, h:pMat; rdiag:pVec);
{--------------------------------------------------------
!    SEPT. 4, 1991
!
!    FIND R^R FOR QR-DECOMPOSED JACOBIAN.
!
!    R IS STORED IN STRICT UPPER TRIANGLE OF A AND RDIAG.
!-------------------------------------------------------}
Var 
    i:integer;
    wv1:pVec;
Begin
  New(wv1);
{ TEMPORARILY REPLACE DIAGONAL OF R IN A (A IS RESTORED LATER). }
  For i:=1 to n do
  begin
    wv1^[i]:=a^[i,i];
    a^[i,i]:=rdiag^[i]
  end;
  utumul (n, n, n, n, n, n, a, h);
  For i:=1 to n do  a^[i,i]:=wv1^[i];
  Dispose(wv1)
End; {rtrmul}

Procedure onenrm (var abort, pertrb:boolean; n, nunit, output:integer;
                  epsmch:double; var h1norm:double; h:pMat; scalex:pVec);
{----------------------------------------------------------------------------
!    FEB. 23, 1992
!
!    FIND 1-NORM OF H MATRIX IF PERTURBATION IS DESIRED AND PERTURB DIAGONAL.
!---------------------------------------------------------------------------}
Label return;
Var
    sqrtep, temp:double;
    i,j:integer;
Begin
  sqrtep:=SQRT(epsmch);
  IF output > 4 THEN
  begin
    Line0;
    Msg('       DIAGONAL OF MATRIX H (= JAC^JAC) BEFORE BEING PERTURBED:');
    Line0;
    For i:=1 to n do
      Writeln(fp_out,'   *          H(',i:3,',',i:3,') = ',h^[i,i]:12:3,'                                     *')
  end;
  h1norm:=zero;
  For j:=1 to n do h1norm:=h1norm + ABS(h^[1,j])/scalex^[j];
  h1norm:=h1norm/scalex^[1];
  For i:=2 to n do
  begin
    temp:=zero;
    For j:=1 to i do temp:=temp + ABS(h^[j,i])/scalex^[j];
    For j:=i+1 to n do temp:=temp + ABS(h^[i,j])/scalex^[j];
    h1norm:=MAX(h1norm,temp/scalex^[i])
  end;
  IF output > 4 THEN
  begin
    Line0;
    Writeln(fp_out,'   *       1-NORM OF MATRIX H: ',h1norm:11:3,'                               *')
  end;
  IF h1norm < epsmch THEN
  begin
    IF output > 0 THEN
    begin
      Line1; Line0;
      Msg('    PROGRAM FAILS AS 1-NORM OF JACOBIAN IS TOO SMALL.');
      Line0; Line1
    end;
    abort:=true;
    goto RETURN
  end
  ELSE
  begin
{   PERTURB DIAGONAL OF MATRIX H - USE THIS TO FIND "SN" }
    pertrb:=true;
    For i:=1 to n do
      h^[i,i]:=h^[i,i] + SQRT(1.0*n)*sqrtep*h1norm*scalex^[i]*scalex^[i];
    IF output > 4 THEN
    begin
      Line0; Line0;
      Msg('    PERTURBED H MATRIX:');
      matprt (n, n, h)
    end;
    abort := FALSE
  end;

Return: End; {onenrm}


Procedure nstpun (var abort, linesr:boolean; overch:boolean; var overfl, qrsing,
                  sclfch: boolean; sclxch:boolean; itnum, maxexp, n, nunit, output:
                  integer; epsmch:double; a:pMat; delf, fvecc:pVec; h: pMat;
                  hhpi:pVec; jac: pMat; rdiag, scalef, scalex, sn: pVec);
{--------------------------------------------------------------------------------
!   FEB. 23, 1992
!
!   THIS SUBROUTINE FINDS THE NEWTON STEP.
!
!   IF THE JACOBIAN IS DETECTED AS SINGULAR OR IF THE ESTIMATED CONDITION
!   NUMBER IS TOO HIGH (GREATER THAN EPSMCH^(-2/3) ) THEN H=J^J IS FORMED
!   AND THE DIAGONAL IS PERTURBED BY ADDING SQRT(N*EPSMCH)*H1NORM*SCALEX(I)^2
!   TO THE CORRESPONDING ELEMENT.   A CHOLESKY DECOMPOSITION IS PERFORMED ON
!   THIS MODIFIED MATRIX PRODUCING A PSEUDO-NEWTON STEP.
!   NOTE: THIS PROCEDURE MAY BE BE BYPASSED FOR ILL-CONDITIONED JACOBIANS
!   BY SETTING THE LOGICAL VARIABLE BYPASS TO TRUE IN THE DRIVER.
!
!   ABORT    IF THE 1-NORM OF MATRIX H BECOMES TOO SMALL
!            ALTHOUGH BUT NOT AT A SOLUTION
!   BYPASS   ALLOWS BYPASSING OF THE SPECIAL TREATMENT FOR
!            BADLY CONDITIONED JACOBIANS
!   QRSING   INDICATES SINGULAR JACOBIAN DETECTED.
!-------------------------------------------------------------------------------}
Label return;
Var
    connum, h1norm, maxadd, maxffl, scalfi, scalxj, sqrtep: double;
    wv1: pVec;
    pertrb: boolean;
    i,j,newstm: integer;
Begin
  New(wv1);
  newstm := 0;

  abort:=false;
  overfl:=false;
  pertrb:=false;
  sqrtep:=SQRT(epsmch);

  IF output > 3 THEN
  begin
    Line0; Line0;
    Msg('    SOLUTION OF LINEAR SYSTEM FOR NEWTON STEP, SN')
  end;

{  STORE (POSSIBLY SCALED) JACOBIAN IN MATRIX A }

  IF (NOT sclfch) THEN
    matcop (n, n, n, n, n, n, jac, a)
  ELSE
    For i:=1 to n do
      IF scalef^[i] <> one THEN
      begin
        scalfi:=scalef^[i];
        For j:=1 to n do a^[i,j]:=jac^[i,j]*scalfi
       {For j:=1 to n do a^[j,i]:=jac^[i,j]*scalfi}
      end
      ELSE
        For j:=1 to n do a^[i,j]:=jac^[i,j];
       {For j:=1 to n do a^[j,i]:=jac^[i,j];}

{  SCALED JACOBIAN IS PRINTED ONLY IF AT LEAST ONE SCALING
   FACTOR IS NOT 1.0 }

  IF (output > 4) AND (sclfch) THEN
  begin
    Line0; Line0;
    Msg('       SCALED JACOBIAN MATRIX:');
    matprt (n, n, a)
  end;

{   QR DECOMPOSITION OF (POSSIBLY SCALED) JACOBIAN }

  qrdcom (qrsing, n, epsmch, a, hhpi, rdiag);

{   SAVE MATRIX A FOR BACK SUBSTITUTION TO CHECK DEUFLHARDS'S SECOND
    STEP ACCEPTANCE CRITERION IN LINE SEARCH OR TRUST REGION METHOD }

  matcop (n, n, n, n, n, n, a, h);

  IF (output > 4) AND (n > 1) AND (NOT matsup) THEN
  begin
    Line0;
    IF (NOT sclfch) THEN
      Msg('       QR DECOMPOSITION OF JACOBIAN MATRIX:')
    ELSE
      Msg('       QR DECOMPOSITION OF SCALED JACOBIAN MATRIX:');
    matprt (n, n, a);
    Line0;
    Msg('           DIAGONAL OF R          PI FACTORS FROM QR DECOMPOSITION');
    Line0;
    For i:=1 to n-1 do
      if Abs(rdiag^[i])<100000 then
        Writeln(fp_out,'   *       RDIAG(',i:3,') = ',rdiag^[i]:11:3,'        HHPI(',i:3,') = ',hhpi^[i]:12:3,'        *')
      else
        Writeln(fp_out,'   *       RDIAG(',i:3,') = ',rdiag^[i]:11:-3,'        HHPI(',i:3,') = ',hhpi^[i]:12:3,'        *');
    if Abs(rdiag^[n])<100000 then
      Writeln(fp_out,'   *       RDIAG(',n:3,') = ',rdiag^[n]:11:3,'        HHPI(',n:3,') = ',hhpi^[n]:12:3,'        *')
    else
      Writeln(fp_out,'   *       RDIAG(',n:3,') = ',rdiag^[n]:11:-3,'        HHPI(',n:3,') = ',hhpi^[n]:12:3,'        *');
    Line0;
    IF itnum = 1 THEN
    begin
      Msg('       NOTE: R IS IN STRICT UPPER TRIANGLE OF MATRIX A PLUS RDIAG.');
      Line0;
      Msg('           THE COLUMNS OF THE LOWER TRIANGLE OF MATRIX A PLUS');
      Msg('           THE ELEMENTS OF VECTOR HHPI FORM THE HOUSEHOLDER');
      Msg('           MATRICES WHICH, WHEN MULTIPLIED TOGETHER, FORM Q.')
    end
  end;

{  ESTIMATE CONDITION NUMBER IF (SCALED) JACOBIAN IS NOT SINGULAR }

  IF (NOT bypass) AND (NOT qrsing) AND (n > 1) THEN
  begin
    IF (sclxch) THEN

{     SET UP FOR CONDITION NUMBER ESTIMATOR - SCALE R WRT X'S } 

      For j:=1 to n do
      begin
        IF scalex^[j] <> one THEN
        begin
          scalxj:=scalex^[j];
          rdiag^[j]:=rdiag^[j]/scalxj;
          For i:=1 to j-1 do a^[i,j]:=a^[i,j]/scalxj
        end
      end;

     condno (overch, overfl, maxexp, n, nunit, output, connum, a, rdiag);

{    IF OVERFLOW DETECTED IN CONDITION NUMBER ESTIMATOR, ASSIGN QRSING AS TRUE }

  IF (overfl) then qrsing:=true;

{    NOTE: OVERFL SWITCHED TO FALSE BEFORE FORMATION OF H LATER }
  end
  ELSE

{    ASSIGN DUMMY TO CONNUM FOR SINGULAR JACOBIAN UNLESS N=1 }

    IF n = 1 THEN
      connum:=one
    ELSE
      connum:=zero;

{    MATRIX H=JAC^JAC IS FORMED IN THREE CASES:
        1)  THE (SCALED) JACOBIAN IS SINGULAR
        2)  THE CONDITION NUMBER IS TOO HIGH AND THE OPTION TO BYPASS THE
            PERTURBATION OF THE JACOBIAN IS NOT BEING USED.
        3)  REQUESTED BY USER THROUGH NEWSTM.  }

IF (qrsing) OR ((NOT bypass) AND (connum > one/Power1(sqrtep,1.3330))) OR (newstm = 77) THEN
begin

{       FORM H=(DF*JAC]^(DF*JAC) WHERE DF=DIAG(SCALEF).  USE
        PREVIOUSLY COMPUTED QR DECOMPOSITION OF (SCALED) JACOBIAN
        WHERE R IS STORED IN THE UPPER TRIANGLE OF A AND RDIAG.  }

  IF (overch) THEN
  begin
    overfl:=false;
    ataov (overfl, maxexp, n, nunit, output, jac, h, scalef)
  end
  ELSE
    rtrmul (n, a, h, rdiag);

  IF output > 3 THEN
  begin
    Line0;
    IF (qrsing) AND (NOT overfl) THEN
      Msg('       SINGULAR JACOBIAN DETECTED: JACOBIAN PERTURBED.')
    ELSE
{     NOTE: IF OVERFL IS TRUE THEN QRSING MUST BE TRUE }
      IF (overfl) THEN
      begin
        Msg('       POTENTIAL OVERFLOW DETECTED IN CONDITION NUMBER ESTIMATOR.');
        Msg('       MATRIX "ASSIGNED" AS SINGULAR AND JACOBIAN PERTURBED.')
      end
      ELSE
        Writeln(fp_out,'   *       CONDITION NUMBER TOO HIGH: ',connum:12:3,' JACOBIAN PERTURBED.','     *')
  end;

  overfl:=false;

  IF newstm <> 77 THEN
  begin 
{   FIND 1-NORM OF H MATRIX AND PERTURB DIAGONAL }
    onenrm (abort, pertrb, n, nunit, output, epsmch, h1norm, h, scalex);
    IF (abort) then goto RETURN
  end;

{  CHOLESKY DECOMPOSITION OF H MATRIX - MAXFFL=0 IMPLIES
   THAT H IS KNOWN TO BE POSITIVE DEFINITE.  }

  maxffl:=zero;

  cholde (n, maxadd, maxffl, sqrtep, h, a);

  IF output > 4 THEN
  begin
    Line0; Line0;
    Msg('    CHOLESKY DECOMPOSITION OF H MATRIX:');
    matprt (n, n, a)
  end;

{   FIND NEWTON STEP FROM CHOLESKY DECOMPOSITION.  IF THE DIAGONAL HAS
    BEEN PERTURBED THEN THIS IS NOT THE ACTUAL NEWTON STEP BUT
    ONLY AN APPORXIMATION THEREOF.  }
  For i:=1 to n do wv1^[i]:=-delf^[i];

  chsolv (overch, overfl, maxexp, n, nunit, output, a, wv1, sn);

  overfl:=false;

  IF output > 3 THEN
    IF (NOT sclxch) THEN
    begin
      Line0;
      IF (NOT pertrb) THEN
        Msg('   NEWTON STEP FROM CHOLESKY DECOMPOSITION')
      ELSE
        Msg('   APPROXIMATE NEWTON STEP FROM PERTURBED JACOBIAN');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *       SN(',i:3,') = ',sn^[i]:12:3,'                                        *')
    end
    ELSE
    begin
      Line0;
      IF (NOT pertrb) THEN
        Msg('   NEWTON STEP FROM CHOLESKY DECOMPOSITION    IN SCALED UNITS:')
      ELSE
        Msg('   APPROXIMATE NEWTON STEP FROM PERTURBED JACOBIAN   IN SCALED UNITS:');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *       SN(',i:3,') = ',sn^[i]:12:3,'            SN(',i:3,') = ',scalex^[i]*sn^[i]:12:3,'       *')
    end;

{   SET QRSING TO TRUE SO THAT THE CORRECT MATRIX FACTORIZATION IS USED
    IN THE BACK-CALCULATION OF SBAR FOR DEUFLHARD RELAXATION FACTOR
    INITIALIZATION IN LINE SEARCH (ONLY MATTERS WHEN JACOBIAN IS ILL-
    CONDITIONED BUT NOT SINGULAR).   }

  qrsing:=true
end
ELSE
begin
  IF (output > 3) AND (n > 1) THEN
  begin
    IF (NOT bypass) AND (connum <= one/Power1(sqrtep,1.333)) THEN
    begin
      Line0;
      Writeln(fp_out,'   *      CONDITION NUMBER ACCEPTABLE: ',connum:9:2,', JACOBIAN NOT PERTURBED   *')
    end;
    IF (bypass) AND (connum > one/Power1(sqrtep,1.333)) THEN
    begin
      Line0;
      Writeln(fp_out,'   *      CONDITION NUMBER HIGH: ',connum:9:2,', JACOBIAN NOT PERTURBED AS      *');
      Msg('      BYPASS IS TRUE.')
    end
  end;

{   NOTE: HERE SN STORES THE R.H.S. - IT IS OVERWRITTEN }

  For i:=1 to n do sn^[i] := -fvecc^[i]*scalef^[i];

  IF (NOT bypass) AND (sclxch) THEN

{   R WAS SCALED BEFORE THE CONDITION NUMBER ESTIMATOR -
    THIS CONVERTS IT BACK TO THE UNSCALED FORM. }
    For j:=1 to n do
      IF scalex^[j] <> one THEN
      begin
        scalxj:=scalex^[j];
        rdiag^[j]:=rdiag^[j]*scalxj;
        For i:= 1 to j-1 do a^[i,j]:=a^[i,j]*scalxj
      end;

{   ACCEPTABLE CONDITION NUMBER - USE BACK SUBSTITUTION TO
    FIND NEWTON STEP FROM QR DECOMPOSITION. }

  qrsolv (overch, overfl, maxexp, n, nunit, output, a, hhpi, rdiag, sn);
  overfl:=false;

  IF output > 3 THEN
    IF (NOT sclxch) THEN
    begin
      Line0;
      Msg('       NEWTON STEP FROM QR DECOMPOSITION:');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *       SN(',i:3,') = ',sn^[i]:12:3,'                                          *')
    end
    ELSE
    begin
      Line0;
      Msg('       NEWTON STEP FROM QR DECOMPOSITION IN SCALED UNITS:');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *       SN(',i:3,') = ',sn^[i]:12:3,'            SN(',i:3,') = ',scalex^[i]*sn^[i]:12:3,
                       '         *')
    end;

{   TRANSFORM MATRICES FOR SUBSEQUENT CALCULATIONS IN TRUST
    REGION METHODS (A IS STORED ABOVE IN H).  }

  IF (NOT linesr) THEN
    For i:=1 to n do
    begin
      a^[i,i]:=rdiag^[i];
      For j:=1 to i-1 do a^[i,j]:=a^[j,i]
    end
end;

return: Dispose(wv1)
End; {nstpun}

Procedure nstpfa (var abort, linesr, overch, overfl, qrsing, restrt, sclfch,
                  sclxch:boolean; var itnum, maxexp, n, newstm, nunit, output:integer;
                  epsmch:double; a:pMat; delf, fvecc:pVec; h:pMat;  hhpi:pVec;
                  jac:pMat; rdiag, scalef, scalex, sn: pVec);
{--------------------------------------------------------------------------------
!
!    FEB. 23, 1992
!
!    THIS SUBROUTINE FINDS THE NEWTON STEP.
!
!    IF THE JACOBIAN IS DETECTED AS SINGULAR OR IF THE ESTIMATED CONDITION
!    NUMBER IS TOO HIGH ^[GREATER THAN EPSMCH**^[-2/3]] THEN H::=J^J IS FORMED
!    AND THE DIAGONAL IS PERTURBED BY ADDING SQRT^[N*EPSMCH]*H1NORM*SCALEX^[I]**2
!    TO THE CORRESPONDING ELEMENT.  A CHOLESKY DECOMPOSITION IS PERFORMED ON
!    THIS MODIFIED MATRIX PRODUCING A PSEUDO-NEWTON STEP.
!
!    IF THE CONDITION NUMBER IS SMALL THEN THE NEWTON STEP, SN,
!    IS FOUND DIRECTLY BY BACK SUBSTITUTION.
!
!    ABORT    IF THE 1-NORM OF MATRIX H BECOMES TOO SMALL
!             ALTHOUGH BUT NOT AT A SOLUTION
!    BYPASS   ALLOWS BYPASSING OF THE SPECIAL TREATMENT FOR
!             BADLY CONDITIONED JACOBIANS
!    QRSING   INDICATES SINGULAR JACOBIAN DETECTED
!--------------------------------------------------------------------------------}
Var
    pertrb: boolean;
    connum, h1norm, maxadd, maxffl, scalfi, scalxj, sqrtep, SUM: double;
    wv1: pVec;
    i,j: integer;
Begin
  New(wv1);
  abort:=false;
  overfl:=false;
  sqrtep:=SQRT(epsmch);

IF (restrt) THEN
begin
  IF output > 3 THEN
  begin
    Line0; Line0;
    Msg('    SOLUTION OF LINEAR SYSTEM FOR NEWTON STEP, SN')
  end;

{   STORE (POSSIBLY SCALED) JACOBIAN IN MATRIX A }

  IF (NOT sclfch) THEN
    matcop (n, n, n, n, n, n, jac, a)
  ELSE
    For i:=1 to n do
      IF scalef^[i] <> one THEN
      begin
        scalfi:=scalef^[i];
        For j:=1 to n do a^[i,j]:=jac^[i,j]*scalfi
      end
      ELSE
        For j:=1 to n do a^[i,j]:=jac^[i,j];

{ SCALED JACOBIAN IS PRINTED ONLY IF AT LEAST ONE SCALING FACTOR IS NOT 1.0 }

  IF (output > 4) AND (sclfch) THEN
  begin
    Line0; Line0;
    Msg('       SCALED JACOBIAN MATRIX:');
    matprt (n, n, a)
  end;

{   QR DECOMPOSITION OF (POSSIBLY SCALED) JACOBIAN }

  qrdcom (qrsing, n, epsmch, a, hhpi, rdiag);

  IF (output > 4) AND (n > 1) AND (NOT matsup) THEN
  begin
    Line0;
    IF (NOT sclfch) THEN
      Msg('       QR DECOMPOSITION OF JACOBIAN MATRIX:')
    ELSE
      Msg('       QR DECOMPOSITION OF SCALED JACOBIAN MATRIX:');
    matprt (n, n, a);
    Line0;
    Msg('          DIAGONAL OF R             PI FACTORS FROM QR DECOMPOSITION');
    Line0;
    For i:=1 to n-1 do
      Writeln(fp_out,'   *       RDIAG(',i:3,') = ',rdiag^[i]:11:3,'        HHPI(',i:3,') = ',hhpi^[i]:12:3,'      *');
    Writeln(fp_out,'   *       RDIAG(',n:3,') = ',rdiag^[n]:11:3,'                                      *');
    Line0;
    IF itnum = 1 THEN
    begin
      Msg('      NOTE: R IS IN STRICT UPPER TRIANGLE OF MATRIX A PLUS RDIAG');
      Line0;
      Msg('           THE COLUMNS OF THE LOWER TRIANGLE OF MATRIX A PLUS');
      Msg('           THE ELEMENTS OF VECTOR HHPI FORM THE HOUSEHOLDER');
      Msg('           MATRICES WHICH, WHEN MULTIPLIED TOGETHER, FORM Q.')
    end
  end;

{   FORM THE ACTUAL Q^ MATRIX FROM THE HOUSEHOLDER TRANSFORMATIONS STORED
    IN THE LOWER TRIANGLE OF A AND THE FACTORS IN HHPI: STORE IT IN JAC }

  qform (n, a, hhpi, jac);

{   COMPLETE THE UPPER TRIANGULAR R MATRIX BY REPLACING THE
    DIAGONAL OF A.  THE QR DECOMPOSITION IS NOW AVAILABLE }

  For i:=1 to n do a^[i,i]:=rdiag^[i];
end
ELSE
begin

{ USING UPDATED FACTORED FORM OF JACOBIAN - CHECK FOR SINGULARITY }

  qrsing:=false;
  For i:=1 to n do IF a^[i,i] = zero then qrsing:=true

end;

{  ESTIMATE CONDITION NUMBER IF JACOBIAN IS NOT SINGULAR }

IF (NOT bypass) AND (NOT qrsing) AND (n > 1) THEN
begin
  IF (sclxch) THEN

{   SET UP FOR CONDITION NUMBER ESTIMATION - SCALE R WRT X'S }

    For j:=1 to n do
      IF scalex^[j] <> one THEN
      begin
        scalxj:=scalex^[j];
        rdiag^[j]:=rdiag^[j]/scalxj;
        For i:=1 to j-1 do a^[i,j]:=a^[i,j]/scalxj
      end;

  condno (overch, overfl, maxexp, n, nunit, output, connum, a, rdiag);

{   UNSCALE R IF IT WAS SCALED BEFORE THE CALL TO CONDNO }

  IF (sclxch) THEN
    For j:=1 to n do
      IF scalex^[j] <> one THEN
      begin
        scalxj:=scalex^[j];
        rdiag^[j]:=rdiag^[j]*scalxj;
        For i:=1 to j-1 do a^[i,j]:=a^[i,j]*scalxj
      end;

{   IF OVERFLOW DETECTED IN CONDITION NUMBER ESTIMATOR ASSIGN
    QRSING AS TRUE SO THAT THE JACOBIAN WILL BE PERTURBED. }

  IF (overfl) then qrsing:=true

{   NOTE: OVERFL SWITCHED TO FALSE BEFORE FORMATION OF H }

end
ELSE
  IF n = 1 THEN
    connum:=one
  ELSE
    connum:=zero;

{   MATRIX H:=JAC^JAC IS FORMED IN TWO CASES:
      1)  THE ^[SCALED] JACOBIAN IS SINGULAR
      2)  THE CONDITION NUMBER IS TOO HIGH AND THE OPTION TO BYPASS
          THE PERTURBATION OF THE JACOBIAN IS NOT BEING USED
      3)  REQUESTED BY THE USER THROUGH NEWSTM.  }

IF (qrsing) OR ((NOT bypass) AND (connum > one/Power1(sqrtep,1.333))) OR (newstm = 77) THEN
begin

{    FORM H=(DF*JAC)^(DF*JAC) WHERE DF=DIAG(SCALEF).  USE
     PREVIOUSLY COMPUTED QR DECOMPOSITION OF (SCALED) JACOBIAN
     WHERE R IS STORED IN THE UPPER TRIANGLE OF A AND RDIAG. }

  overfl:=false;
  IF (overch) THEN
    ataov (overfl, maxexp, n, nunit, output, jac, h, scalef)
  ELSE
    rtrmul (n, a, h, rdiag);

  IF output > 3 THEN
  begin
    Line0;
    IF (qrsing) AND (NOT overfl) THEN
      Msg('       SINGULAR JACOBIAN DETECTED: JACOBIAN PERTURBED.')
    ELSE
    begin
      IF (overfl) THEN
      begin
        Msg('       POTENTIAL OVERFLOW DETECTED IN CONDITION NUMBER ESTIMATOR,');
        Msg('       MATRIX "ASSIGNED" AS SINGULAR AND JACOBIAN PERTURBED.')
      end
      ELSE
        Writeln(fp_out,'   *       CONDITION NUMBER TOO HIGH: ',connum:12:3,' JACOBIAN PERTURBED.   *')
    end
  end;
  overfl:=false;

  IF newstm <> 77 THEN
{   FIND 1-NORM OF H MATRIX AND PERTURB DIAGONAL }
    onenrm (abort, pertrb, n, nunit, output, epsmch, h1norm, h, scalex);

{   CHOLESKY DECOMPOSITION OF H MATRIX - MAXFFL=0 INDICATES
    THAT H IS KNOWN TO BE POSITIVE DEFINITE. }

  maxffl:=zero;
  cholde (n, maxadd, maxffl, sqrtep, h, a);

  IF output > 4 THEN
  begin
    Line0; Line0;
    Msg('    CHOLESKY DECOMPOSITION OF H MATRIX:');
    matprt (n, n, a)
  end;

{   FIND NEWTON STEP FROM CHOLESKY DECOMPOSITION.  IF THE
    DIAGONAL HAS BEEN PERTURBED THEN THIS IS NOT THE ACTUAL
    NEWTON STEP BUT ONLY AN APPROXIMATION THEREOF. }

  For i:=1 to n do wv1^[i] := -delf^[i];
  chsolv (overch, overfl, maxexp, n, nunit, output, a, wv1, sn);
  overfl:=false;

  IF output > 3 THEN
    IF (NOT sclxch) THEN
    begin
      IF (NOT pertrb) THEN
        Msg('     NEWTON STEP FROM CHOLESKY DECOMPOSITION')
      ELSE
        Msg('     APPROXIMATE NEWTON STEP FROM PERTURBED JACOBIAN');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *       SN(',i:3,') = ',sn^[i]:12:3,'                                        *')
    end
    ELSE
    begin
      Line0;
      IF (NOT pertrb) THEN
        Msg('     NEWTON STEP FROM CHOLESKY DECOMPOSITION IN SCALED UNITS')
      ELSE
        Msg('     APPROXIMATE NEWTON STEP FROM PERTURBED JACOBIAN   IN SCALED UNITS');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *       SN(',i:3,') = ',sn^[i]:12:3,'   SN(',i:3,') = ',scalex^[i]*sn^[i]:12:3,'               *')
    end;

{     SET QRSING TO TRUE SO THAT THE CORRECT MATRIX FACTORIZATION
      IS USED IN THE BACK-CALCULATION OF SBAR FOR DEUFLHARD
      RELAXATION FACTOR INITIALIZATION.  }

  qrsing:=true
end
ELSE
begin
  IF (output > 3) AND (n > 1) THEN
  begin
    IF (NOT bypass) AND (connum <= one/Power1(sqrtep,1.333)) THEN
    begin
      Line0;
      Writeln(fp_out,'   *     CONDITION NUMBER ACCEPTABLE: ',connum:9:2,', JACOBIAN NOT PERTURBED  *')
    end;
    IF (bypass) AND (connum > one/Power1(sqrtep,1.333)) THEN
    begin
      Line0;
      Writeln(fp_out,'   *     CONDITION NUMBER HIGH: ',connum:9:2,', JACOBIAN NOT PERTURBED AS     *');
      Msg('     BYPASS IS TRUE.')
    end
  end;
  For i:=1 to n do
  begin
    sum:=zero;
    For j:=1 to n do sum:=sum - jac^[i,j]*scalef^[j]*fvecc^[j];
    sn^[i]:=sum
  end;

  rsolv (overch, overfl, maxexp, n, nunit, output, a, rdiag, sn);
  overfl:=false;

  IF output > 3 THEN
  begin
    IF (NOT sclxch) THEN
    begin
      Line0;
      Msg('       NEWTON STEP FROM QR DECOMPOSITION');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *       SN(',i:3,') = ',sn^[i]:12:3,'                                        *')
    end
    ELSE
    begin
      Line0;
      Msg('       NEWTON STEP FROM QR DECOMPOSITION IN SCALED UNITS');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *       SN(',i:3,') = ',sn^[i]:12:3,'   SN(',i:3,') = ',scalex^[i]*sn^[i]:12:3,'               *')
    end
  end;

{ TRANSFORM MATRICES FOR SUBSEQUENT CALCULATIONS IN TRUST REGION METHOD }

  IF (NOT linesr) THEN
    For i:=2 to n do
      For j:=1 to i-1 do
        a^[i,j]:=a^[j,i]
end;
Dispose(wv1)
End; {nstpfa}

Procedure title (cauchy, deuflh, geoms, linesr, newton, overch: boolean;
                 acptcr, contyp, itsclf, itsclx, jactyp, jupdm, maxit, maxns,
                 maxqns, mgll, minqns, n, narmij, ninitn, njacch, nunit,
                 output, qnupdm, stopcr, trupdm:integer; alpha, confac, delfac,
                 delta, epsmch, etafac, fcnold, ftol, lam0, maxstp, mstpf,
                 nsttol, omega, ratiof, sigma, stptol:double; boundl, boundu,
                 fvecc, scalef, scalex, xc: pVec);
{--------------------------------------------------------
!    APR. 13, 1991
!
!    THIS PROCEDURE WRITES TITLE AND RECORDS PARAMETERS.
!-------------------------------------------------------}
Label 780,return;
Var
    i: integer;
Begin
  IF output < 2 then goto RETURN;
  Writeln(fp_out); Writeln(fp_out); Writeln(fp_out);
  Line1; Line1; Line0;
  Msg('                           N N E S');
  Line0;
  Msg('       NONMONOTONIC NONLINEAR EQUATION SOLVER VERSION 1.05');
  Line0;
  Msg('                 COPYRIGHT 1991, BY R.S. BAIN');
  Line0; Line1; Line1;
  Writeln(fp_out); Writeln(fp_out);

  IF output < 3 then GOTO 780;

  Line1; Line1;

IF (newton) THEN
begin
  IF jupdm = 0 THEN
    Msg('  METHOD: NEWTON (NO LINE SEARCH)')
  ELSE IF jupdm = 1 THEN
    Msg('  METHOD: QUASI-NEWTON (NO LINE SEARCH) USING BROYDEN UPDATE')
  ELSE IF jupdm = 2 THEN
    Msg('  METHOD: QUASI-NEWTON (NO LINE SEARCH) USING LEE AND LEE UPDATE');
  Line0;
  IF overch THEN
    Msg('  OVERLOW CHECKING IN USE.')
  ELSE
    Msg('  OVERFLOW CHECKING NOT IN USE.');
  IF jactyp = 0 THEN
    IF njacch > 0 THEN
      Writeln(fp_out,'   *  ANALYTICAL JACOBIAN USED, CHECKED NUMERICALLY, NJACCH: ',njacch:5,'         *')
    ELSE
      Msg('  ANALYTICAL JACOBIAN USED; NOT CHECKED.')
  ELSE IF jactyp = 1 THEN
    Msg('  JACOBIAN ESTIMATED USING FORWARD DIFFERENCES.')
  ELSE IF jactyp = 2 THEN
    Msg('  JACOBIAN ESTIMATED USING BACKWARD DIFFERENCES.')
  ELSE
    Msg('  JACOBIAN ESTIMATED USING CENTRAL DIFFERENCES.');
  Line0; Line1
end
ELSE
begin
  IF (linesr) THEN
  begin
    Line0;
    IF (deuflh) THEN
      Msg('  DEUFLHARD RELAXATION FACTOR INITIALIZATION IN EFFECT.')
    ELSE
      Msg('  DEUFLHARD RELAXATION FACTOR INITIALIZATION NOT IN EFFECT.')
  end
  ELSE
  begin
    IF etafac = one THEN
      Msg('  METHOD: TRUST REGION USING SINGLE DOGLEG STEPS.')
    ELSE
      Msg('  METHOD: TRUST REGION USING DOUBLE DOGLEG STEPS.');
    Line0;
    IF (cauchy) THEN
      Msg('  INITIAL STEP CONSTRAINED BY SCALED CAUCHY STEP.')
    ELSE
      Msg('  INITIAL STEP CONSTRAINED BY SCALED NEWTON STEP.')
  end;
  IF (geoms) THEN
    Msg('  METHOD: GEOMETRIC SEARCH.')
  ELSE
    Msg('  METHOD: SEARCH BASED ON SUCCESSIVE MINIMIZATIONS.');
  IF (overch) THEN
    Msg('  OVERLOW CHECKING IN USE.')
  ELSE
    Msg('  OVERFLOW CHECKING NOT IN USE.');
  IF jupdm = 0 THEN
    Msg('  NO QUASI-NEWTON UPDATE USED.');
  IF jupdm = 1 THEN
    IF qnupdm = 0 THEN
      Msg('  BROYDEN QUASI-NEWTON UPDATE OF UNFACTORED JACOBIAN.')
    ELSE
      Msg('  BROYDEN QUASI-NEWTON UPDATE OF FACTORED JACOBIAN.')
  ELSE IF jupdm = 2 THEN
    IF qnupdm = 0 THEN
      Msg('  LEE AND LEE QUASI-NEWTON UPDATE OF UNFACTORED JACOBIAN.')
    ELSE
      Msg('  LEE AND LEE QUASI-NEWTON UPDATE OF FACTORED JACOBIAN.');
  IF jactyp = 0 THEN
    IF njacch > 0 THEN
      Writeln(fp_out,'   *  ANALYTICAL JACOBIAN USED, CHECKED NUMERICALLY, NJACCH: ',njacch:5,'         *')
    ELSE
      Msg('  ANALYTICAL JACOBIAN USED; NOT CHECKED.')
  ELSE IF jactyp = 1 THEN
    Msg('  JACOBIAN ESTIMATED USING FORWARD DIFFERENCES.')
  ELSE IF jactyp = 2 THEN
    Msg('  JACOBIAN ESTIMATED USING BACKWARD DIFFERENCES.')
  ELSE
    Msg('  JACOBIAN ESTIMATED USING CENTRAL DIFFERENCES.');
  IF (NOT linesr) THEN
  begin
    Line0;
    IF (trupdm = 0) AND (jupdm > 0) THEN
      Msg('  TRUST REGION UPDATED USING POWELL STRATEGY.')
    ELSE
      Msg('  TRUST REGION UPDATED USING DENNIS AND SCHNABEL STRATEGY.');
  end;
  Line0; Line1; Line0;
  IF itsclf <> 0 THEN
  begin
    Writeln(fp_out,'   *  ADAPTIVE FUNCTION SCALING STARTED AT ITERATION: ..........',itsclf:6,'     *');
    Line0
  end;
  IF itsclx <> 0 THEN
  begin
    Writeln(fp_out,'   *  ADAPTIVE VARIABLE SCALING STARTED AT ITERATION: ..........',itsclx:6,'     *');
    Line0
  end;
  IF (linesr) THEN
    IF jupdm = 0 THEN
      Writeln(fp_out,'   *  MAXIMUM NUMBER OF STEPS IN LINE SEARCH, MAXNS: ...........',maxns:6,'     *')
    ELSE
    begin
      Writeln(fp_out,'   *  MAXIMUM NUMBER OF NEWTON LINE SEARCH STEPS, MAXNS: .......',maxns:6,'     *');
      Writeln(fp_out,'   *  MAXIMUM NUMBER OF QUASI-NEWTON LINE SEARCH STEPS, MAXQNS: ',maxqns:6,'     *')
    end
  ELSE
    IF jupdm = 0 THEN
      Writeln(fp_out,'   *  MAXIMUM NUMBER OF TRUST REGION UPDATES, MAXNS: ...........',maxns:6,'     *')
    ELSE
    begin
      Writeln(fp_out,'   *  MAXIMUM NO. OF NEWTON TRUST REGION UPDATES, MAXNS: .......',maxns:6,'     *');
      Writeln(fp_out,'   *  MAXIMUM NO. OF QUASI-NEWTON TRUST REGION UPDATES, MAXQNS: ',maxqns:6,'     *')
    end;
  IF narmij < maxit THEN
  begin
    Line0;
    Writeln(fp_out,'   *  NUMBER OF OBJECTIVE FUNCTION VALUES COMPARED, MGLL: ......',mgll:6,'     *')
  end;
  IF jupdm > 0 THEN
  begin
    IF narmij = maxit Then Line0;
    Writeln(fp_out,'   *  MINIMUM NUMBER OF STEPS BETWEEN JACOBIAN UPDATES, MINQNS: ',minqns:6,'     *');
    Writeln(fp_out,'   *  NUMBER OF NON-QUASI-NEWTON STEPS AT START, NINITN: .......',ninitn:6,'     *')
  end;
  Writeln(fp_out,'   *  NUMBER OF ARMIJO STEPS AT START, NARMIJ: .................',narmij:6,'     *')
end;
Line0;
IF stopcr = 3 THEN
  Writeln(fp_out,'   *  FUNCTION AND STEP SIZE STOPPING CRITERIA, STOPCR: ........',stopcr:6,'     *')
ELSE IF stopcr = 12 THEN
  Writeln(fp_out,'   *  FUNCTION OR STEP SIZE STOPPING CRITERIA, STOPCR: .........',stopcr:6,'     *')
ELSE IF stopcr = 1 THEN
  Writeln(fp_out,'   *  STEP SIZE STOPPING CRITERION, STOPCR: ....................',stopcr:6,'     *')
ELSE
  Writeln(fp_out,'   *  FUNCTION STOPPING CRITERION, STOPCR: .....................',stopcr:6,'     *');
IF (NOT newton) THEN
begin
  Line0;
  IF acptcr = 12 THEN
    Writeln(fp_out,'   *  FUNCTION AND STEP SIZE ACCEPTANCE CRITERIA, ACPTCR: ......',acptcr:6,'     *')
  ELSE IF acptcr = 2 THEN
    Writeln(fp_out,'   *  STEP SIZE ACCEPTANCE CRITERION, ACPTCR: ..................',acptcr:6,'     *')
  ELSE
    Writeln(fp_out,'   *  FUNCTION ACCEPTANCE CRITERION, ACPTCR: ...................',acptcr:6,'     *');
  IF contyp <> 0 THEN
  begin
    Line0;
    Writeln(fp_out,'   *  CONSTRAINTS IN USE, CONTYP: ..............................',contyp:6,'     *')
  end
end;
Line0; Line1; Line0;
Writeln(fp_out,'   *  ESTIMATED MACHINE EPSILON, EPSMCH: ......',epsmch,'     *');
Line0;
Writeln(fp_out,'   *  FACTOR TO ESTABLISH MAXIMUM STEP SIZE, MSTPF : .......',mstpf:10:3,'     *');
Writeln(fp_out,'   *  CALCULATED MAXIMUM STEP SIZE, MAXSTP: ................',maxstp:10:3,'     *');
IF (NOT linesr) THEN
begin
  IF delta < zero THEN
  begin
    Line0;
    Msg('  INITIAL TRUST REGION NOT PROVIDED.')
  end
  ELSE
  begin
    Line0;
    Writeln(fp_out,'   *  INITIAL TRUST REGION SIZE, DELTA: ....................',delta:10:3,'     *')
  end;
  IF etafac < one THEN
    Writeln(fp_out,'   *  FACTOR TO SET DIRECTION OF TRUST REGION STEP, ETAFAC: ... ',etafac:6:4,'     *');
  Writeln(fp_out,'   *  TRUST REGION UPDATING FACTOR, DELFAC: ................',delfac:10:3,'     *')
end;
IF (NOT newton) THEN
begin
  Line0;
  Writeln(fp_out,'   *  FACTOR IN OBJECTIVE FUNCTION COMPARISON, ALPHA: ......',alpha:10:5,'     *');
  IF (linesr) AND (NOT newton) THEN
  begin
    Line0;
    Writeln(fp_out,'   *  REDUCTION FACTOR FOR RELAXATION FACTOR, SIGMA: .......',sigma:10:3,'     *')
  end;
  IF jupdm <> 0 THEN
  begin
    Line0;
    Writeln(fp_out,'   *  REDUCTION REQUIRED IN OBJ. FUNCTION FOR QN STEP, RATIOF:  ',ratiof:6:4,'     *')
  end;
  IF jupdm = 2 THEN
  begin
    Line0;
    Writeln(fp_out,'   *  FACTOR IN LEE AND LEE UPDATE, OMEGA: .....................',omega:6:4,'     *')
  end
end;
Line0;
IF stopcr <> 2 THEN
begin
  Writeln(fp_out,'   *  STOPPING TOLERANCE FOR STEP SIZE, STPTOL (E-10): .....',1E10*stptol:10:3,'     *');
  Writeln(fp_out,'   *  STOPPING TOLERANCE FOR NEWTON STEP, NSTTOL (E-10): ...',1E10*nsttol:10:3,'     *')
end;
IF stopcr <> 1 THEN
  Writeln(fp_out,'   *  STOPPING TOLERANCE FOR OBJECTIVE FUNCTION, FTOL(E-5)..',1E5*ftol:10:3,'     *');
IF (linesr) AND (NOT newton) AND (lam0 < one) THEN
begin
  Line0;
  Writeln(fp_out,'   *  INITIAL LAMBDA IN LINE SEARCH, LAM0: .................',lam0:10:3,'     *')
end;
IF contyp > 0 THEN
begin
  Line0;
  Writeln(fp_out,'   *  FACTOR TO ENSURE STEP WITHIN CONSTRAINTS, CONFAC: ........',confac:6:4,'     *')
end;
Line0; Line1; Line0;
Msg('  SCALING FACTORS');
Line0;
Msg('      COMPONENT VALUES                   FUNCTION VALUES');
Line0;
For i:=1 to n do
  Writeln(fp_out,'   *  SCALEX(',i:3,') = ',scalex^[i]:10:3,'      SCALEF(',i:3,') = ',scalef^[i]:10:3,'               *');
IF contyp > 0 THEN
begin
  Line0; Line1; Line0;
  Msg('  LOWER AND UPPER BOUNDS');
  Line0;
  Msg('      LOWER BOUNDS                          UPPER BOUNDS');
  Line0;
  For i:=1 to n do
    Writeln(fp_out,'   *  BOUNDL(',i:3,') = ',boundl^[i]:10:3,'      BOUNDU(',i:3,') = ',boundu^[i]:10:3,'               *')
end;
Line0;

780:IF output = 2 Then Line1;
  Line1; Line0;
  Msg('    INITIAL ESTIMATES               INITIAL FUNCTION VALUES');
  Line0;
  For i:=1 to n do
    if abs(fvecc^[i]) < 100000 then
      Writeln(fp_out,'   *  X(',i:3,') = ',xc^[i]:12:3,'             F(',i:3,') = ',fvecc^[i]:12:3,'              *')
    else
      Writeln(fp_out,'   *  X(',i:3,') = ',xc^[i]:12:3,'             F(',i:3,') = ',fvecc^[i]:12:-3,'              *');
  Line0;
  if Abs(fcnold) < 100000 then
    Writeln(fp_out,'   *  INITIAL OBJECTIVE FUNCTION VALUE = ',fcnold:12:3,'                      *')
  else
    Writeln(fp_out,'   *  INITIAL OBJECTIVE FUNCTION VALUE = ',fcnold:12:-3,'                      *');
  Line0; Line1; Line1;

return: End; {title}

Procedure qmin (nunit, output:integer; delfts:double; var delta:double;
                deltaf, stplen:double);
{-----------------------------------------------------------------------
!    FEB. 9, 1991
!
!    SET THE NEW TRUST REGION SIZE, DELTA, BASED ON A QUADRATIC
!    MINIMIZATION WHERE DELTA IS THE INDEPENDENT VARIABLE.
!
!    DELTAF IS THE DIFFERENCE IN THE SUM-OF-SQUARES OBJECTIVE
!    FUNCTION VALUE AND DELFTS IS THE DIRECTIONAL DERIVATIVE IN
!    THE DIRECTION OF THE CURRENT STEP, S, WHICH HAS STEP LENGTH STPLEN.
!----------------------------------------------------------------------}
Const point1 = 0.1;
Var
    deltmp: double;
Begin

  IF deltaf-delfts <> zero THEN
  begin

{  CALCULATE DELTA WHERE MINIMUM WOULD OCCUR - DELTMP.
   THIS IS PROVISIONAL AS IT MUST BE WITHIN CERTAIN LIMITS TO BE ACCEPTED }

    deltmp:=-delfts*stplen/(two*(deltaf-delfts));
    IF output > 4 THEN
    begin
      Line0;
      Writeln(fp_out,'   *       TEMPORARY DELTA FROM QUADRATIC MINIMIZATION: ',deltmp:12:3,'     *');
      Writeln(fp_out,'   *                            VERSUS CURRENT DELTA: ',delta:12:3,'       *')
    end;

{   REDUCE DELTA DEPENDING ON THE MAGNITUDE OF DELTMP.
    IT MUST BE WITHIN (0.1*DELTA,0.5*DELTA) TO BE ACCEPTED -
    OTHERWISE THE NEAREST ENDPOINT OF THE INTERVAL IS USED }

    IF deltmp < point1*delta THEN
    begin
      delta:=point1*delta;
      IF output > 4 THEN
      begin
        Line0;
        Msg('       NEW DELTA SET TO 0.1 CURRENT DELTA')
      end
    end
    ELSE IF deltmp > delta/two THEN
    begin
      delta:=delta/two;
      IF output > 4 THEN
      begin
        Line0;
        Msg('       NEW DELTA SET TO 0.5 CURRENT DELTA')
      end
    end
    ELSE
    begin
      delta:=deltmp;
      IF output > 4 THEN
      begin
        Line0;
        Msg('       NEW DELTA SET TO DELTMP')
      end
    end
  end
  ELSE
  begin
    IF output > 4 THEN
    begin
      Line0;
      Msg('       TO AVOID OVERFLOW NEW DELTA SET TO 0.5 CURRENT DELTA.')
    end;
    delta:=delta/two
  end

End; {qmin}


Procedure trstup (var geoms, newtkn, overch, overfl, qrsing, sclfch, sclxch:boolean;
                  var acpcod, acpstr, acptcr, contyp, isejac, jupdm, maxexp,
                  mgll, mnew, n, narmij, nfunc, notrst, nunit, output,
                  qnupdm, retcod, trupdm:integer; var alpha, confac, delfac, delstr,
                  delta:double; epsmch:double; var fcnmax, fcnnew, fcnold, fcnpre, maxstp,
                  newlen, newmax, powtau, rellen, stptol:double;  a, astore: pMat;
                  boundl, boundu, delf, fplpre, ftrack, fvec, fvecc, hhpi: pVec; jac:pMat;
                  rdiag, rhs, s, sbar, scalef, scalex, strack, xc, xplpre, xplus:pVec);
{-------------------------------------------------------------
!    FEB. 28, 1992
!
!    THIS SUBROUTINE CHECKS FOR ACCEPTANCE OF A TRUST REGION
!    STEP GENERATED BY THE DOUBLE DOGLEG METHOD.
!------------------------------------------------------------}
Label return;
Const
    pt5 = 0.5; threeq = 0.75; onept1 = 1.1;
Var
    convio: boolean;
    i, j, k: integer;
    delfpr, delfts, deltaf, dlftsm, dmult, powlam, powmu, ratio, ratiom,
    sbrnrm, sp, ss, stplen, SUM: double;
    V1, wv3: pVec; temp1,temp2: pMat;
Begin
  New(temp1); New(temp2); New(V1); New(wv3);

{ NOTE: ACCEPTANCE CODE, ACPCOD, IS 0 ON ENTRANCE TO TRSTUP}

  convio:=false;
  overfl:=false;
  retcod := 0;

IF output > 3 THEN
begin
  Line0; Line0;
  IF (NOT sclfch) AND (NOT sclxch) THEN
    Msg('    TRUST REGION UPDATING')
  ELSE
    Msg('    TRUST REGION UPDATING (ALL X''S AND F''S IN UNSCALED UNITS)');
  Line0
end;

{ CHECK TO MAKE SURE "S" IS A DESCENT DIRECTION - FIND DIRECTIONAL
  DERIVATIVE AT CURRENT XC USING S GENERATED BY DOGLEG SUBROUTINE. }

innerp (overch, overfl, maxexp, n, n, n, nunit, output, delfts, delf, s);

IF output > 3 THEN
begin
  Line0;
  Writeln(fp_out,'   *       INNER PRODUCT OF DELF AND S, DELFTS: ........',delfts:13:4,'      *')
end;
IF delfts > zero THEN
begin
  IF output > 3 THEN
  begin
    Line0;
    Msg('       DIRECTIONAL DERIVATIVE POSITIVE, SEARCH DIRECTION REVERSED.')
  end;
  For i:=1 to n do s^[i] := -s^[i]
end;

{   FIND MAXIMUM OBJECTIVE FUNCTION VALUE AND MAXIMIUM STEP
    LENGTH FOR NONMONOTONIC SEARCH.  THIS HAS TO BE DONE ONLY
    ONCE DURING EACH ITERATION (WHERE NOTRST=1). }

IF notrst = 1 THEN
begin
  newmax:=newlen;
  fcnmax:=fcnold;
  IF isejac > narmij THEN
    IF (isejac < narmij+mgll) THEN
      For j:=1 to mnew do
      begin
        fcnmax:=MAX(fcnmax, ftrack^[j-1]);
        newmax:=MAX(newmax, strack^[j-1])
      end
    ELSE
      For j:=0 to mnew do
      begin
        fcnmax:=MAX(fcnmax, ftrack^[j]);
        newmax:=MAX(newmax, strack^[j])
      end
end;

{     TEST TRIAL POINT - FIND XPLUS AND TEST FOR CONSTRAINT
      VIOLATIONS IF CONTYP DOES NOT EQUAL 0.  }

For  i:=1 to n do
begin
  wv3^[i]:=-one;

{ WV3 IS A MARKER FOR "VIOLATORS" - IT CHANGES TO 1 OR 2 }

  xplus^[i]:=xc^[i] + s^[i];
  IF contyp > 0 THEN
    IF xplus^[i] < boundl^[i] THEN
    begin
      convio:=true;
      wv3^[i]:=one
    end
    ELSE IF xplus^[i] > boundu^[i] THEN
    begin
      convio:=true;
      wv3^[i]:=two
    end
end;

{   IF CONSTRAINT IS VIOLATED... }

IF (convio) THEN
begin
  IF output > 3 THEN
  begin
    Line0;
    Msg('       CONSTRAINT VIOLATED:');
    Line0;
    Msg('      TRIAL ESTIMATES (VIOLATIONS MARKED):');
    Line0;
    For i:=1 to n do
      IF wv3^[i] > zero THEN
        Writeln(fp_out,'   *            XPLUS(',i:3,') = ',xplus^[i]:12:3,'  *                             *')
      ELSE
        Writeln(fp_out,'   *            XPLUS(',i:3,') = ',xplus^[i]:12:3,'                                *')
  end;

{     FIND STEP WITHIN CONSTRAINED REGION.

      FIND THE RATIO OF THE DISTANCE FROM THE (I)TH COMPONENT TO ITS
      CONSTRAINT TO THE LENGTH OF THE PROPOSED STEP, XPLUS(I)-XC(I).
      MULTIPLY THIS BY CONFAC (DEFAULT 0.95) TO ENSURE THE NEW STEP STAYS
      WITHIN THE ACCEPTABLE REGION UNLESS XC IS CLOSE TO THE BOUNDARY
      (RATIO <= 1/2).   IN SUCH CASES A FACTOR OF 0.5*CONFAC IS USED.

      NOTE: ONLY THE VIOLATING COMPONENTS ARE REDUCED.  }

  ratiom:=one;

{     RATIOM STORES THE MINIMUM VALUE OF RATIO }

  For i:=1 to n do
  begin
    IF wv3^[i] = one THEN
      ratio:=(boundl^[i]-xc^[i])/s^[i]
    ELSE IF wv3^[i] = two THEN
      ratio:=(boundu^[i]-xc^[i])/s^[i];

    IF wv3^[i] > zero THEN
    begin
{   NOTE: RATIO IS STORED IN WV3 FOR OUTPUT ONLY }

      wv3^[i]:=ratio;

      ratiom:=MIN(ratiom,ratio);
      IF ratio > pt5 THEN
        xplus^[i]:=xc^[i] + confac*ratio*s^[i]
      ELSE
{   WITHIN BUFFER ZONE }
        xplus^[i]:=xc^[i] + confac*ratio*s^[i]/two;
      s^[i]:=xplus^[i] - xc^[i]
    end
  end;

  IF output > 3 THEN
  begin
    Line0; Line0;
    Msg('       NEW S AND XPLUS VECTORS (WITH RATIOS FOR VIOLATIONS)');
    Line0;
    Msg('       NOTE: RATIOS ARE RATIO OF LENGTH TO BOUNDARY FROM CURRENT');
    Msg('       X VECTOR TO MAGNITUDE OF CORRESPONDING PROPOSED STEP.');
    Line0;
    For i:=1 to n do
      IF wv3^[i] < zero THEN
        Writeln(fp_out,'   *       S(',i:3,') = ',s^[i]:12:3,'    XPLUS(',i:3,') = ',xplus^[i]:12:3,'            *')
      ELSE
        Writeln(fp_out,'   *       S(',i:3,') = ',s^[i]:12:3,'    XPLUS(',i:3,') = ',xplus^[i]:12:3,' ',wv3^[i]:9:3,' *');
    Line0;
    Writeln(fp_out,'   *       MINIMUM OF RATIOS, RATIOM: ',ratiom:12:3,'                       *')
  end;

{   THE NEW POINT, XPLUS, IS NOT NECESSARILY IN A DESCENT DIRECTION.
    CHECK DIRECTIONAL DERIVATIVE FOR MODIFIED STEP, DLFTSM. }

  innerp (overch, overfl, maxexp, n, n, n, nunit, output, dlftsm, delf, s);

  IF output > 3 THEN
  begin
    Line0;
    Writeln(fp_out,'   *       INNER PRODUCT OF DELF AND MODIFIED S, DL FTSM: ',dlftsm:12:3,'     *')
  end;

{   IF DLFTSM IS POSITIVE REDUCE TRUST REGION.  IF NOT, TEST NEW POINT }

  IF dlftsm > zero THEN
  begin
    delta:=confac*ratiom*delta;
    retcod:=8;
    goto RETURN
  end
end; {if convio}

{ CONSTRAINTS NOT (OR NO LONGER) VIOLATED - TEST NEW POINT }

  fcn1 (overfl, n, fvec, xplus);
  nfunc:=nfunc+1;

{ IF OVERFLOW AT NEW POINT REDUCE TRUST REGION AND RETURN }

IF (overfl) THEN
begin

{   IF THE OVERFLOW COMES AS A RESULT OF INCREASING DELTA WITHIN THE
    CURRENT ITERATION (IMPLYING DELSTR IS POSITIVE) AND DIVIDING DELTA
    BY 10 WOULD PRODUCE A DELTA WHICH IS SMALLER THAN THAT AT THE
    STORED POINT, THEN USE STORED POINT AS THE UPDATED RESULT. }

  IF delstr > delta/ten THEN
  begin
    For i:=1 to n do
    begin
      temp1^[i,1]:=xplpre^[i]; temp2^[i,1]:=xplus^[i];
    end;
    matcop (n, n, 1, 1, n, 1, temp1, temp2);
    For i:=1 to n do xplus^[i]:=temp2^[i,1];
    For i:=1 to n do
    begin
      temp1^[i,1]:=fplpre^[i]; temp2^[i,1]:=fvec^[i];
    end;
    matcop (n, n, 1, 1, n, 1, temp1, temp2);
    For i:=1 to n do fvec^[i]:=temp2^[i,1];
    acpcod:=acpstr;
    delta:=delstr;
    fcnnew:=fcnpre;
    retcod:=1
  end
  ELSE
  begin
    delta:=delta/ten;
    retcod:=9
  end;
  goto RETURN
end;

{ NO OVERFLOW IN RESIDUAL VECTOR }

IF output > 3 THEN
begin
  Line0;
  Msg('          TRIAL ESTIMATES                 FUNCTION VALUES');
  Line0;
  For i:=1 to n do
    Writeln(fp_out,'   *      XPLUS(',i:3,') = ',xplus^[i]:12:3,'       FVEC(',i:3,') = ',fvec^[i]:12:3,'       *');
end;

{   IF NO OVERFLOW WITHIN RESIDUAL VECTOR FIND OBJECTIVE FUNCTION }

  fcnevl (overfl, maxexp, n, nunit, output, epsmch, fcnnew, fvec, scalef);

{ IF OVERFLOW IN OBJECTIVE FUNCTION EVALUATION REDUCE TRUST REGION AND RETURN }

IF (overfl) THEN
begin
{   IF THE OVERFLOW COMES AS A RESULT OF INCREASING DELTA WITHIN THE
    CURRENT ITERATION (SO THAT DELSTR IS POSITIVE) AND DIVIDING DELTA
    BY 10 WOULD PRODUCE A DELTA WHICH IS SMALLER THAN THAT AT THE STORED
    POINT THEN USE STORED POINT AS THE UPDATED RESULT.  }

  IF delstr > delta/ten THEN
  begin
    For i:=1 to n do
    begin
      temp1^[i,1]:=xplpre^[i]; temp2^[i,1]:=xplus^[i];
    end;
    matcop (n, n, 1, 1, n, 1, temp1, temp2);
    For i:=1 to n do xplus^[i]:=temp2^[i,1];
    For i:=1 to n do
    begin
      temp1^[i,1]:=fplpre^[i]; temp2^[i,1]:=fvec^[i];
    end;
    matcop (n, n, 1, 1, n, 1, temp1, temp2);
    For i:=1 to n do fvec^[i]:=temp2^[i,1];
    acpcod:=acpstr;
    delta:=delstr;
    fcnnew:=fcnpre;
    retcod:=2
  end
  ELSE
  begin
    delta:=delta/ten;
    retcod:=10
  end;
  goto RETURN
end
ELSE
begin

{ NO OVERFLOW AT TRIAL POINT - COMPARE OBJECTIVE FUNCTION TO FCNMAX }

  IF output > 3 THEN
  begin
    Line0;
    IF (NOT sclfch) THEN
      Writeln(fp_out,'   *       OBJECTIVE FUNCTION AT XPLUS, FCNNEW: .........',fcnnew:12:4,'      *')
    ELSE
      Writeln(fp_out,'   *       SCALED OBJECTIVE FUNCTION AT XPLUS, FCNNEW: ..',fcnnew:12:4,'      *');
    Writeln(fp_out,'   *       COMPARE TO FCNMAX+ALPHA*DELFTS: ..............',(fcnmax + alpha*delfts):12:4,'      *')
  end
end;

{   IF ACPTCR=12 CHECK SECOND DEUFLHARD STEP ACCEPTANCE TEST
    BY FINDING 2-NORM OF SBAR.  THERE ARE FOUR POSSIBILITIES
    DEPENDING ON WHETHER THE JACOBIAN IS OR IS NOT SINGULAR
    AND WHETHER QNUPDM IS 0 OR 1 }

IF acptcr = 12 THEN
begin
  IF (qrsing) THEN
  begin
{   FORM -J^F AS RIGHT HAND SIDE - METHOD DEPENDS ON
    WHETHER QNUPDM EQUALS 0 OR 1 (UNFACTORED OR FACTORED) }

    IF qnupdm = 0 THEN
    begin
{   UNSCALED JACOBIAN IS IN MATRIX JAC }

      For k:=1 to n do wv3^[k] := -fvec^[k]*scalef^[k]*scalef^[k];
      mmulv(n, jac, wv3, rhs)
    end
    ELSE
    begin
{   R IN UPPER TRIANGLE OF A PLUS RDIAG AND Q^ IN JAC
    FROM QR DECOMPOSITION OF SCALED JACOBIAN. }
      For i:=1 to n do
      begin
        sum:=zero;
        For j:=1 to n do sum:=sum - jac^[i,j]*fvec^[j]*scalef^[j];
        wv3^[i]:=sum
      end;
      rhs^[1]:=rdiag^[1]*wv3^[1];
      For j:=2 to n do
      begin
        For k:=1 to j-1 do V1^[k]:=a^[k,j];
        sum:=DOT_PRODUCT(j-1, V1, wv3);
        rhs^[j]:=sum + rdiag^[j]*wv3^[j]
      end
    end;
    chsolv (overch, overfl, maxexp, n, nunit, output, a, rhs, sbar)
  end
  ELSE
  begin

{   RIGHT HAND SIDE IS -FVEC }

    IF (qnupdm = 0) OR (jupdm = 0) THEN
    begin

{   QR DECOMPOSITION OF SCALED JACOBIAN STORED IN ASTORE }

      For k:=1 to n do sbar^[k] := -fvec^[k]*scalef^[k];
      qrsolv (overch, overfl, maxexp, n, nunit, output, astore, hhpi, rdiag, sbar)
    end
    ELSE
    begin

{     SET UP RIGHT HAND SIDE - MULTIPLY -FVEC BY Q^
      (STORED IN JAC).  RHS IS A WORK VECTOR ONLY HERE }

      For k:=1 to n do wv3^[k] := -fvec^[k]*scalef^[k];
      avmul (n, n, n, n, jac, wv3, sbar);
      rsolv (overch, overfl, maxexp, n, nunit, output, a, rdiag, sbar)
    end
  end;

{   NORM OF (SCALED) SBAR IS NEEDED FOR SECOND ACCEPTANCE TEST }

  For k:=1 to n do wv3^[k]:=scalex^[k]*sbar^[k];
  twonrm (overfl, maxexp, n, epsmch, sbrnrm, wv3);

  IF output > 4 THEN
  begin
    Line0;
    IF (NOT sclxch) THEN
    begin
      Msg('          DEUFLHARD SBAR VECTOR');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *          SBAR(',i:3,') = ',sbar^[i]:12:3,'                                   *')
    end
    ELSE
    begin
      Msg('          DEUFLHARD SBAR VECTOR              IN SCALED X UNITS       *');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *          SBAR(',i:3,') = ', sbar^[i]:12:3,'        SBAR(',i:3,') = ',
                (scalex^[i]*sbar^[i]):12:3,'   *')
    end
  end;

  IF output > 3 THEN
  begin
    Line0;
    IF (NOT sclxch) THEN
      Writeln(fp_out,'   *          VALUE OF SBRNRM AT XPLUS: .................',sbrnrm:12:4,'    *')
    ELSE
      Writeln(fp_out,'   *          VALUE OF SCALED SBRNRM AT XPLUS: ..........',sbrnrm:12:4,'    *');
    Writeln(fp_out,'   *          NEWMAX: ...................................',newmax:12:4,'    *')
  end;

  IF sbrnrm < newmax then acpcod:=2;

{ FUNCTION VALUE ACCEPTANCE IS ALSO CHECKED REGARDLESS
  OF WHETHER SECOND STEP ACCEPTANCE CRITERION WAS MET }

end; {if acptor=12}

{ ESTABLISH DELTAF FOR USE IN COMPARISON TO PREDICTED
  CHANGE IN OBJECTIVE FUNCTION, DELFPR, LATER }

  deltaf:=fcnnew - fcnold;

IF fcnnew >= fcnmax + alpha*delfts THEN
begin

{   FAILURE OF FIRST STEP ACCEPTANCE TEST. TEST LENGTH OF
    STEP TO ENSURE PROGRESS IS STILL BEING MADE }

  rellen:=zero;
  For i:=1 to n do
    rellen:=MAX(rellen, ABS(s^[i])/MAX((ABS(xplus^[i])),one/scalex^[i]));

  IF rellen < stptol THEN
  begin

{ NO PROGRESS BEING MADE - RETCOD = 7 STOPS PROGRAM }
    For i:=1 to n do
    begin
      temp1^[i,1]:=xc^[i]; temp2^[i,1]:=xplus^[i];
    end;
    matcop (n, n, 1, 1, n, 1, temp1, temp2);
    For i:=1 to n do xplus^[i]:=temp2^[i,1];
    retcod:=7;
    goto RETURN
  end
  ELSE
  begin

{   FAILURE OF STEP BY OBJECTIVE FUNCTION CRITERION.
    ESTABLISH A NEW DELTA FROM EITHER SIMPLE DIVISION
    BY DELFAC OR BY FINDING THE MINIMUM OF A QUADRATIC MODEL }

    IF (geoms) THEN
      delta:=delta/delfac
    ELSE
    begin
{   FIRST FIND LENGTH OF TRUST REGION STEP }
      For i:=1 to n do wv3^[i]:=s^[i]*scalex^[i];
      twonrm (overfl, maxexp, n, epsmch, stplen, wv3);
      qmin (nunit, output, delfts, delta, deltaf, stplen)
    end;

    IF (delta < delstr) THEN
    begin

{     IF DELTA HAS BEEN INCREASED AT THIS ITERATION AND THE
      DELTA FROM QMIN IS LESS THAN THE DELTA AT THE PREVIOUSLY
      ACCEPTED (STORED) POINT THEN RETURN TO THAT POINT AND
      ACCEPT IT AS THE UPDATED ITERATE. }

      For i:=1 to n do
      begin
        temp1^[i,1]:=xplpre^[i]; temp2^[i,1]:=xplus^[i];
      end;
      matcop (n, n, 1, 1, n, 1, temp1, temp2);
      For i:=1 to n do xplus^[i]:=temp2^[i,1];
      For i:=1 to n do
      begin
        temp1^[i,1]:=fplpre^[i]; temp2^[i,1]:=fvec^[i];
      end;
      matcop (n, n, 1, 1, n, 1, temp1, temp2);
      For i:=1 to n do fvec^[i]:=temp2^[i,1];
      acpcod:=acpstr;
      delta:=delstr;
      fcnnew:=fcnpre;
      retcod:=3;
      goto RETURN
    end;

{     IF THE SECOND ACCEPTANCE TEST HAS BEEN PASSED RETURN
      WITH NEW TRUST REGION AND CONTINUE ON TO NEXT ITERATION;
      OTHERWISE TRY A NEW STEP WITH REDUCED DELTA. }

    IF acpcod = 2 THEN
      retcod:=4
    ELSE
{     FAILURE OF FIRST STEP ACCEPTANCE TEST }
      retcod:=11;
    goto RETURN
  end
end
ELSE
begin

{     OBJECTIVE FUNCTION MEETS FIRST ACCEPTANCE CRITERION.
      IN NONMONOTONIC SEARCHES IT MAY BE GREATER THAN THE
      PREVIOUS OBJECTIVE FUNCTION VALUE - CONSIDER THIS CASE FIRST. }

  IF deltaf >= alpha*delfts THEN
  begin
{   AN ACCEPTABLE STEP HAS BEEN FOUND FOR THE NONMONOTONIC SEARCH
    BUT THE OBJECTIVE FUNCTION VALUE IS NOT A "DECREASE" FROM THE
    PREVIOUS ITERATION (ACTUALLY IT MIGHT BE BETWEEN ZERO AND
    ALPHA*DELFTS).  ACCEPT STEP BUT REDUCE DELTA. }
    delta:=delta/delfac;
    retcod:=5;
    IF acpcod = 2 THEN
      acpcod:=12
    ELSE
      acpcod:=1;
    goto RETURN
  end;

{   COMPARE DELTAF TO DELTAF PREDICTED, DELFPR, TO DETERMINE NEXT
    TRUST REGION SIZE.  NOTE: DELTAF MUST BE LESS THAN ALPHA*DELFTS
    (IN ESSENCE NEGATIVE) TO HAVE REACHED THIS POINT IN TRSTUP.
    R IS IN UPPER TRIANGLE OF MATRIX A SO THE FOLLOWING CODE FINDS:

    DELFPR = DELF^S + 1/2 S^J^JS = DELF^S + 1/2 S^R^RS }

  uvmul (n, n, n, n, a, s, wv3);
  delfpr := delfts + DOT_PRODUCT(n, wv3, wv3) / two;
  IF output > 4 THEN
  begin
    Line0;
    Writeln(fp_out,'   *       PREDICTED CHANGE IN OBJECTIVE FUNCTION, DELFPR:',delfpr:12:3,'   *');
    Writeln(fp_out,'   *          ACTUAL CHANGE IN OBJECTIVE FUNCTION, DELTAF:',deltaf:12:3,'   *')
  end;

  IF ((retcod <= 6) AND (ABS(delfpr-deltaf) <= ABS(deltaf)/ten)) OR ((deltaf <= delfts) AND (NOT newtkn) AND (NOT convio))
         AND (delstr = zero) THEN
  begin
    IF MIN(newlen,maxstp)/delta > onept1 THEN
    begin
{           PROMISING STEP - INCREASE TRUST REGION.
            STORE CURRENT POINT.  }

      For i:=1 to n do
      begin
        temp1^[i,1]:=xplus^[i]; temp2^[i,1]:=xplpre^[i];
      end;
      matcop (n, n, 1, 1, n, 1, temp1, temp2);
      For i:=1 to n do xplpre^[i]:=temp2^[i,1];
      For i:=1 to n do
      begin
        temp1^[i,1]:=fvec^[i]; temp2^[i,1]:=fplpre^[i];
      end;
      matcop (n, n, 1, 1, n, 1, temp1, temp2);
      For i:=1 to n do fplpre^[i]:=temp2^[i,1];
      delstr:=delta;
      fcnpre:=fcnnew;

{     IF NONMONOTONIC STEPS ARE BEING USED EXPAND TRUST
      REGION TO NEWLEN, OTHERWISE EXPAND BY DELFAC. }

      IF isejac > narmij THEN
        delta:=MIN(newlen,maxstp)
      ELSE
        delta:=MIN(delfac*delta,maxstp);
      retcod:=12;
      IF acpcod = 2 THEN
        acpstr:=12
      ELSE
        acpstr:=1;
      acpcod:=0
    end
    ELSE
    begin
      retcod:=0;
      IF acpcod = 2 THEN
        acpcod:=12
      ELSE
        acpcod:=1;
    end;
    goto RETURN
  end
  ELSE
  begin

{   CHANGE TRUST REGION SIZE DEPENDING ON DELTAF AND DELFPR }

    retcod:=6;
    IF acpcod = 2 THEN
      acpcod:=12
    ELSE
      acpcod:=1;
    IF deltaf >= delfpr/ten THEN
    begin
      delta:=delta/delfac;
      IF output > 3 THEN
      begin
        Line0;
        Msg('    CHANGE IN F, DELTAF, IS > 0.1 DELFPR - REDUCE DELTA.')
      end
    end
    ELSE IF (trupdm = 0) AND (jupdm > 0) THEN
    begin

{     POWELL'S UPDATING SCHEME - FIND JAC S FIRST }

      IF qnupdm =0 THEN
      begin
{       UNSCALED JACOBIAN IN JAC }
        For i:=1 to n do rhs^[i] := s^[i]*scalef^[i];
        avmul (n, n, n, n, jac, rhs, wv3)
      end
      ELSE
      begin
{       MULTIPLY BY R FIRST }
        uvmul (n, n, n, n, a, s, rhs);
{       THEN Q (IN JAC^) }
        mmulv (n, jac, rhs, wv3)
      end;
      dmult:=delfpr/ten - deltaf;
      sp:=zero;
      ss:=zero;
      For k:=1 to n do
      begin
        wv3^[k]:=wv3^[k] + fvecc^[k];
        sp:=sp + ABS(fvec^[k]*(fvec^[k] - wv3^[k]));
        ss:=ss + (fvec^[k]-wv3^[k])*(fvec^[k] - wv3^[k])
      end;
      IF sp + SQRT(sp*sp + dmult*ss) < epsmch THEN
        powlam:=ten
      ELSE
        powlam:=one + dmult/(sp + SQRT(sp*sp + dmult*ss));
      powlam:=SQRT(powlam);
      powmu:=MIN(MIN(delfac,powlam), powtau);
      powtau:=powlam/powmu;
      IF output > 3 THEN
      begin
        Line0;
        Msg('      FACTORS IN POWELL UPDATING SCHEME:');
        Line0;
        Writeln(fp_out,'   *      LAMBDA: ',powlam:12:3,'    MU: ',powmu:12:3,'    TAU: ', powtau:12:3,'  *');
        Msg('       DELTA IS MINIMUM OF MU*DELTA AND MAXSTP.')
      end;
      delta:=MIN(powmu*delta, maxstp)
    end
    ELSE
    begin
      IF deltaf < threeq*delfpr THEN
      begin
        delta:=MIN(delfac*delta, maxstp);
        IF output > 3 THEN
        begin
          Line0;
          Msg('       CHANGE IN F, DELTAF, IS LESS THAN 0.75 X PREDICTED.');
          Writeln(fp_out,'   *       DELTA INCREASED TO: ',delta:12:3,'                              *')
        end
      end
      ELSE
        IF output > 3 THEN
        begin
          Line0;
          Msg('       DELTAF BETWEEN 0.1 AND 0.75 DELFPR - LEAVE DELTA UNCHANGED')
        end
    end
  end
end; {if fcnnew...}

return: Dispose(temp1); Dispose(temp2); Dispose(V1); Dispose(wv3)
End; {trstup}

Procedure rcdprt (nunit, retcod:integer; delta, rellen, stptol:double);
{-------------------------------------------------------------------------
!    FEB. 14, 1991
!
!    DESCRIBE MEANING OF RETURN CODES, RETCOD, FROM TRUST REGION UPDATING.
!------------------------------------------------------------------------}
Begin

  Line0;
  Writeln(fp_out,'   *       RETCOD, FROM TRUST REGION UPDATING:',retcod:5,'                      *');
  Line0;
  IF retcod = 1 THEN
  begin
    Msg('       PROMISING STEP FOUND; DELTA HAS BEEN INCREASED TO NEWLEN BUT');
    Msg('       BECAUSE OF OVERFLOWS IN THE FUNCTION VECTOR(S) IN SUBSEQUENT');
    Msg('       STEP(S) THE PROJECTED DELTA IS LESS THAN THAT AT THE ALREADY');
    Msg('       SUCCESSFUL STEP - RETURN TO SUCCESSFUL STEP AND ACCEPT AS');
    Msg('       NEW POINT.')
  end
  ELSE IF retcod = 2 THEN
  begin
    Msg('       PROMISING STEP FOUND; DELTA HAS BEEN INCREASED TO NEWLEN BUT');
    Msg('       BECAUSE OF OVERFLOWS IN THE OBJECTIVE FUNCTION IN SUBSEQUENT');
    Msg('       STEP(S) THE PROJECTED DELTA IS LESS THAN THAT AT THE ALREADY');
    Msg('       SUCCESSFUL STEP - RETURN TO SUCCESSFUL STEP AND ACCEPT AS');
    Msg('       NEW POINT.')
  end
  ELSE IF retcod = 3 THEN
  begin
    Msg('       PROMISING STEP FOUND; DELTA HAS BEEN INCREASED TO NEWLEN BUT');
    Msg('       BECAUSE OF SUBSEQUENT FAILURES IN THE STEP ACCEPTANCE TEST^[S)');
    Msg('       THE PROJECTED DELTA IS LESS THAN THAT AT THE ALREADY');
    Msg('       SUCCESSFUL STEP - RETURN TO SUCCESSFUL STEP AND ACCEPT.')
  end
  ELSE IF retcod = 4 THEN
    Msg('       STEP ACCEPTED BY STEP SIZE CRITERION ONLY - DELTA REDUCED')
  ELSE IF retcod = 5 THEN
  begin
    Msg('       STEP ACCEPTED - NEW FUNCTION VALUE > PREVIOUS =>');
    Msg('       REDUCE TRUST REGION.')
  end
  ELSE IF retcod = 6 THEN
    Msg('       STEP ACCEPTED - DELTA CHANGED AS DETAILED ABOVE.')
  ELSE IF retcod = 7 THEN
  begin
    Msg('       NO PROGRESS MADE: RELATIVE STEP SIZE IS TOO SMALL');
    Writeln(fp_out,'   *       REL. STEP SIZE, RELLEN = ',rellen:12:3,', STPTOL = ',stptol:12:3,'   *')
  end
  ELSE IF retcod = 8 THEN
  begin
    Msg('       POINT MODIFIED BY CONSTRAINTS NOT A DESCENT DIRECTION');
    Msg('       DELTA REDUCED TO CONFAC*RATIOM*DELTA.')
  end
  ELSE IF retcod = 9 THEN
    Msg('       OVERFLOW DETECTED IN FUNCTION VECTOR - DELTA REDUCED.')
  ELSE IF retcod = 10 THEN
    Msg('       OVERFLOW IN OBJECTIVE FUNCTION - DELTA REDUCED.')
  ELSE IF retcod = 11 THEN
  begin
    Msg('       STEP NOT ACCEPTED - REDUCE TRUST REGION SIZE BY MINIMIZATION');
    Msg('       OF QUADRATIC MODEL IN STEP DIRECTION.')
  end
  ELSE
    Msg('       PROMISING STEP - INCREASE DELTA TO NEWLEN AND TRY A NEW STEP.');
  Line0;
  Writeln(fp_out,'   *       DELTA ON RETURN FROM TRUST REGION UPDATING: ',delta:11:3,'       *')
End; {rcdprt}


Procedure nnes (var absnew, cauchy, deuflh, geoms, linesr, newton, overch:boolean;
                 acptcr:integer; var itsclf, itsclx:integer; jactyp, jupdm:integer;
                 var maxexp:integer; maxit, maxns, maxqns, mgll, minqns, n, narmij,
                 niejev, njacch:integer; var njetot:integer; nunit:integer; var output:
                 integer; qnupdm:integer; var stopcr:integer; supprs:integer; var trmcod,
                 trupdm:integer; var alpha, confac, delta, delfac, epsmch, etafac, fcnnew,
                 fdtolj, ftol, lam0, mstpf, nsttol, omega, ratiof, sigma, stptol:double;
                 a:pMat; boundl, boundu, delf, fsave, ftrack, fvec, fvecc:pVec; h:pMat;
                 hhpi:pVec; jac, plee:pMat; rdiag, s, sbar, scalef, scalex, sn, ssdhat,
                 strack, vhat, xc, xplus, xsave: pVec; help: String);
{----------------------------------------------------------------------------
!    FEB. 28, 1992
!    *** Solve a non-linear system - Main procedure ***
!---------------------------------------------------------------------------}
Label 310, 320, 330, return;
Const pt01 = 0.01;
Var
    beta, caulen, delstr, delta0, fcnmax, fcnmin, fcnold, fcnpre,
    fstore, maxstp, newlen, newmax, powtau, rellen, sbrnrm, sqrtz,
    stpmax, xnorm: double;
    wv1, wv2, wv4: pVec;
    acpcod, acpstr, contyp, countr, i, isejac, itemp, itstr, j,
    maxlin, maxtrs, mnew, mold, nac1, nac2, nac12, nfail, nfestr,
    nfetot, nosubt, notrst, outtmp, retcod: integer;
    abort, checkj, frstdg, instop, jacerr, newtkn, overfl, qnfail,
    qrsing, restrt, savest, sclfch, sclxch: boolean;
    newstm: integer;
Begin

  newstm:=0;

  New(wv1); New(wv2); New(wv4);

{ PRINT HELP IF REQUESTED }
  IF help <> 'NONE' THEN
  begin
    {olhelp (nunit, help); }
    goto RETURN
  end;

  qnfail:=false;
  restrt:=true;
  savest:=false;
  countr:=0;
  isejac:=0;
  mnew:=0;
  nac1:=0;
  nac2:=0;
  nac12:=0;
  nfunc:=1;
  njetot:=0;
  outtmp:=output;
  trmcod:=0;
  retcod := 0;

{ ESTABLISH INITIAL FUNCTION VALUE AND CHECK FOR STARTING ESTIMATE WHICH
  IS A SOLUTION.  ALSO, CHECK FOR INCOMPATIBILITIES IN INPUT PARAMETERS. }

  initch (instop, linesr, newton, overfl, sclfch, sclxch, acptcr, contyp,
          jactyp, jupdm, maxexp, n, nunit, output, qnupdm, stopcr,
          trupdm, epsmch, fcnold, ftol, boundl, boundu, fvecc, scalef,
          scalex, xc);

{ IF A FATAL ERROR IS DETECTED IN INITCH RETURN TO MAIN PROGRAM.
  NOTE: SOME INCOMPATIBILITIES ARE CORRECTED WITHIN INITCH AND EXECUTION
  CONTINUES.   WARNINGS ARE GENERATED WITHIN INITCH.  }

IF (instop) then goto RETURN;

{  ESTABLISH MAXIMUM STEP LENGTH ALLOWED (USUALLY THIS IS MUCH
   LARGER THAN ACTUAL STEP SIZES - IT IS ONLY TO PREVENT
   EXCESSIVELY LARGE STEPS).  THE FACTOR MSTPF CONTROLS THE
   MAGNITUDE OF MAXSTP AND IS SET BY THE USER (DEFAULT=1000). }

  maxst (overfl, maxexp, n, nunit, output, epsmch, maxstp, mstpf, scalex, xc);

{  WRITE TITLE AND RECORD PARAMETERS }

  title (cauchy, deuflh, geoms, linesr, newton, overch, acptcr, contyp,
         itsclf, itsclx, jactyp, jupdm, maxit, maxns, maxqns, mgll,
         minqns, n, narmij, niejev, njacch, nunit, output, qnupdm, stopcr,
         trupdm, alpha, confac, delfac, delta, epsmch, etafac, fcnold,
         ftol, lam0, maxstp, mstpf, nsttol, omega, ratiof, sigma, stptol,
         boundl, boundu, fvecc, scalef, scalex, xc);

{ INITIALIZE FTRACK AND STRACK VECTORS (FTRACK STORES (TRACKS) THE FUNCTION
  VALUES FOR THE NONMONOTONIC COMPARISON AND STRACK, SIMILARLY, STORES THE
  LENGTH OF THE NEWTON STEPS - WHICH ARE USED IN CONJUNCTION WITH
  DEUFLHARD'S SECOND ACCEPTANCE CRITERION).  }

  For j:=0 to mgll-1 do
  begin
    ftrack^[j]:=zero;
    strack^[j]:=zero
  end;

{ MAIN ITERATIVE LOOP - MAXIT IS SPECIFIED BY USER.

  ITNUM COUNTS OVERALL ITERATIONS.
  ISEJAC COUNTS ITERATIONS SINCE LAST EXPLICIT JACOBIAN
  EVALUATION IF A QUASI-NEWTON METHOD IS BEING USED.  }

For itnum:=1 to maxit do
begin

{ SUPPRESS OUTPUT IF DESIRED - USED IF DETAILED OUTPUT IS DESIRED
  FOR LATER ITERATIONS ONLY (DEFAULT VALUE FOR SUPPRS IS 0 -
  I.E. NO SUPPRESSION).  }

  IF itnum < supprs THEN
    output:=3
  ELSE
    output:=outtmp;

  IF output > 2 THEN
  begin
    Line1; Line0;
    Writeln(fp_out,'   *  ITERATION NUMBER: ',itnum:5,'                                              *')
  end;

{  UPDATE ITERATION COUNTER, ISEJAC, FOR QUASI-NEWTON METHODS
    ONLY (I.E. JUPDM > 0).  }

  IF (jupdm > 0) AND (restrt) AND (NOT newton) THEN
    isejac:=1
  ELSE
    isejac:=isejac+1;

  IF (output > 4) AND (jupdm > 0) AND (NOT newton) THEN
  begin
    Line0;
    IF (restrt) THEN
      IF itnum > niejev THEN Msg('    RESTRT IS TRUE, ISEJAC SET TO 1.')
    ELSE
      Writeln(fp_out,'   *  # OF ITERATIONS SINCE EXPLICIT JACOBIAN, ISEJAC, INCREASED TO ',isejac:4,'*')
  end;

{   WHEN AN EXPLICIT JACOBIAN IS BEING USED IN QUASI-NEWTON METHODS THEN
    MAXNS STEPS ARE ALLOWED IN THE LINE SEARCH.   FOR STEPS BASED ON A
    QUASI-NEWTON APPROXIMATION MAXQNS ARE ALLOWED.  SIMILARLY MAXNS AND
    MAXQNS STEPS ARE ALLOWED IN TRUST REGION METHODS RESPECTIVELY.
    THIS IS AN ATTEMPT TO AVOID AN EXCESSIVE NUMBER OF FUNCTION
    EVALUATIONS IN A DIRECTION WHICH WOULD NOT LEAD TO A SIGNIFICANT
    REDUCTION. }

  IF (NOT newton) THEN
  begin
    IF (linesr) THEN
      IF jupdm > 0 THEN
        IF isejac = 1 THEN

{   JACOBIAN UPDATED EXPLICITLY }
          maxlin:=maxns

        ELSE

{   QUASI-NEWTON UPDATE }
          maxlin:=maxqns
      ELSE
        maxlin:=maxns
    ELSE
      IF jupdm > 0 THEN
        IF isejac = 1 THEN

{   JACOBIAN UPDATED EXPLICITLY }
          maxtrs:=maxns

        ELSE

{   QUASI-NEWTON UPDATE }
          maxtrs:=maxqns

      ELSE
        maxtrs:=maxns;

    IF (jupdm > 0) AND (output > 4) THEN
    begin
      Line0;
      IF (linesr) THEN
        Writeln(fp_out,'   *    MAXLIN SET TO :',maxlin:5,'                                                *')
      ELSE
        Writeln(fp_out,'   *    MAXTRS SET TO: ',maxtrs:5,'                                                *')
    end
  end;

{   ESTABLISH WHETHER JACOBIAN IS TO BE CHECKED NUMERICALLY -
    NJACCH ESTABLISHES THE NUMBER OF ITERATIONS FOR WHICH JACOBIAN
    CHECKING IS DESIRED.  IF CHECKJ IS TRUE, A FORWARD DIFFERENCE NUMERICAL
    APPROXIMATION OF THE JACOBIAN IS COMPARED TO THE ANALYTICAL VERSION.
    STORE THE NUMBER OF FUNCTION EVALUATIONS SO THAT THESE "EXTRA" ARE NOT
    INCLUDED IN OVERALL STATISTICS.  }

  IF jactyp = 0 THEN
    IF itnum > njacch THEN checkj:=false
    ELSE
    begin
      checkj:=true;
      nfestr:=nfetot
    end
  ELSE
    checkj:=false;

{   EVALUATE JACOBIAN AT FIRST STEP OR IF NO QUASI-NEWTON
    UPDATE IS BEING USED (RESTRT IS FALSE ONLY IN QUASI-
    NEWTON METHODS WHEN THE QUASI-NEWTON UPDATE IS BEING USED). }

  IF (restrt) OR (jupdm = 0) THEN
  begin

{   IF MORE THAN ONE DAMPED NEWTON STEP HAS BEEN REQUESTED AT THE START,
    IDENTIFY THIS AS THE REASON FOR EXPLICIT JACOBIAN EVALUATION }

    IF (jupdm > 0) AND (itnum <= niejev) AND (itnum > 1) AND (output > 4) THEN
    begin
      Line0;
      Msg('    AS ITNUM <= NIEJEV JACOBIAN EVALUATED EXPLICITLY.')
    end;

{     OTHERWISE JUPDM IS 0 OR RESTRT IS TRUE AND ITNUM IS > NIEJEV }

    IF (jupdm > 0) AND (itnum > 1) AND (output > 4) AND (itnum > niejev) THEN
    begin
      Line0;
      Msg('    RESTRT IS TRUE - JACOBIAN EVALUATED EXPLICITLY.')
    end;

{   NOTE: MATRIX H IS USED HERE TO HOLD THE FINITE DIFFERENCE
    ESTIMATION USED IN CHECKING THE ANALYTICAL JACOBIAN IF CHECKJ
    IS TRUE.   VECTORS WV1 AND, FOR CENTRAL DIFFERENCES, WV2
    TEMPORARILY HOLD THE FINITE DIFFERENCE FUNCTION EVALUATIONS. }

    jacobi (checkj, jacerr, overfl, jactyp, n, nunit, output, epsmch,
            fdtolj, boundl, boundu, fvecc, wv1, wv2, jac, h, scalex, xc);

    njetot:=njetot + 1;

{   RETURN IF ANALYTICAL AND NUMERICAL JACOBIANS DON'T AGREE
    (APPLICABLE ONLY IF CHECKJ IS TRUE AND A DISCREPANCY IS FOUND
    WITHIN SUBROUTINE JACOBI). A WARNING IS GIVEN FROM WITHIN JACOBI }

    IF (jacerr) then goto RETURN;

{   RESET TOTAL NUMBER OF FUNCTION EVALUATIONS TO NEGLECT
    THOSE USED IN CHECKING ANALYTICAL JACOBIAN. }

    IF (checkj) then nfetot:=nfestr;

    IF jupdm > 0 THEN
    begin 
{   FCNMIN IS THE MINIMUM OF THE OBJECTIVE FUNCTION FOUND SINCE
    THE LAST EXPLICIT JACOBIAN EVALUATION.  IT IS USED TO ESTABLISH
    WHICH STEP THE PROGRAM RETURNS TO WHEN A QUASI-NEWTON STEP FAILS. }

      fcnmin:=fcnold;

{   POWTAU IS THE TAU FROM POWELL'S TRUST REGION UPDATING SCHEME
    (USED WHEN JUPDM=1).  IT IS RESET TO 1.0 AT EVERY EXPLICIT
    JACOBIAN EVALUATION IN QUASI-NEWTON METHODS.  }

      powtau:=one;

{   AT EVERY EXPLICIT JACOBIAN EVALUATION EXCEPT POSSIBLY THE FIRST,
    IN QUASI-NEWTON METHODS,  A NEW TRUST REGION IS CALCULATED
    INTERNALLY USING EITHER THE NEWTON STEP OR THE CAUCHY STEP AS
    SPECIFIED BY THE USER IN LOGICAL VARIABLE "CAUCHY" IN
    SUBROUTINE DELCAU.  THIS IS FORCED BY SETTING DELTA TO A
    NEGATIVE NUMBER.  }

      IF itnum > 1 then delta:=-one;

{   RESET COUNTER FOR NUMBER OF FAILURES OF RATIO TEST (IS
    FCNNEW/FCNOLD < RATIOF?) SINCE LAST EXPLICIT JACOBIAN UPDATE }

      nfail:=0;
      IF (output > 4) AND (itnum > niejev) THEN
      begin
        Line0;
        Msg('    NUMBER OF FAILURES OF RATIO TEST, NFAIL, SET BACK TO 0.')
      end;

{     ESTABLISH A NEW MAXIMUM STEP LENGTH ALLOWED }

      IF itnum > 1 then maxst (overfl, maxexp, n, nunit, output, epsmch, maxstp, mstpf, scalex, xc);

{     SET "P" MATRIX TO IDENTITY FOR LEE AND LEE UPDATE (JUPDM=2). }

      IF (jupdm = 2) AND (itnum >= niejev) THEN
        For j:=1 to n do
        begin
          For i:=1 to n do plee^[i,j]:=zero;
          plee^[j,j]:=one
        end
    end
  end;

  IF ((restrt) OR (qnupdm = 0)) AND (output > 4) AND (NOT matsup) THEN
  begin
{   WRITE JACOBIAN MATRIX }
    Line0;
    Msg('    JACOBIAN MATRIX:');
    matprt (n, n, jac)
  end;

{   ESTABLISH SCALING MATRICES IF DESIRED (ITSCLF=0 => NO ADAPTIVE
    SCALING, WHILE ITSCLF > 0 MEANS ADAPTIVE SCALING STARTS AFTER THE
    ^[ITSCLF]TH ITERATION).  SIMILARLY FOR COMPONENT SCALING...
    NOTE: SCALING FACTORS ARE UPDATED ONLY WHEN THE JACOBIAN IS UPDATED
    EXPLICITLY IN QUASI-NEWTON METHODS.

    FUNCTION SCALING.   }

  IF (restrt) AND (itsclf > 0) AND (itsclf <= itnum) THEN
  begin

    ascalf (n, epsmch, fvecc, jac, scalef);

    IF output > 4 THEN
    begin
      Line0;
      Msg('    FUNCTION SCALING VECTOR:');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *       SCALEF(',i:3,') = ',scalef^[i]:12:3,'                                    *')
    end;

{   RECALCULATE OBJECTIVE FUNCTION WITH NEW SCALING FACTORS.
    THIS AVOIDS PREMATURE FAILURES IF THE CHANGE IS SCALING FACTORS
    WOULD MAKE THE PREVIOUS OBJECTIVE FUNCTION VALUE SMALLER.  }

    fcnevl (overfl, maxexp, n, nunit, output, epsmch, fcnold, fvecc, scalef)
  end;

{   COMPONENT SCALING }

  IF (restrt) AND (itsclx > 0) AND (itsclx <= itnum) THEN
  begin

    ascalx (n, epsmch, jac, scalex);

    IF output > 4 THEN
    begin
      Line0;
      Msg('       COMPONENT SCALING VECTOR:');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *       SCALEX(',i:3,') = ',scalex^[i]:12:3,'                                    *')
    end
  end;

{   FIND GRADIENT OF 1/2 FVECC^FVECC (NOT USED IF DECOMPOSED
    MATRIX IS UPDATED IN WHICH CASE THE GRADIENT IS FOUND
    WITHIN THAT SUBROUTINE - CALL IS MADE FOR OUTPUT ONLY). }

  gradf (overch, overfl, restrt, sclfch, sclxch, jupdm, maxexp, n,
         nunit, output, qnupdm, delf, fvecc, jac, scalef, scalex);

{   FIND NEWTON STEP USING QR DECOMPOSITION.

    IF JUPDM = 0 OR QNUPDM = 0 THEN THE UNFACTORED FORM
    FOR CALCULATING THE NEWTON STEP IS USED.  }

  IF (jupdm = 0) OR (qnupdm = 0) THEN

{       NEWTON STEP - UNFACTORED FORM BEING UPDATED }

    nstpun(abort, linesr, overch, overfl, qrsing, sclfch, sclxch,
           itnum, maxexp, n, nunit, output, epsmch, a, delf, fvecc, h,
           hhpi, jac, rdiag, scalef, scalex, sn)
  ELSE

{   NEWTON STEP - FACTORED FORM BEING UPDATED }

    nstpfa (abort, linesr, overch, overfl, qrsing, restrt, sclfch,
            sclxch, itnum, maxexp, n, newstm, nunit, output, epsmch, a,
            delf, fvecc, h, hhpi, jac, rdiag, scalef, scalex, sn);

{   RUN IS ABORTED IF THE JACOBIAN BECOMES ESSENTIALLY ALL ZEROS. }

  IF (abort) then goto return;

{  CHECK FOR CONVERGENCE ON LENGTH OF NEWTON STEP IF
   STOPCR = 1, 12 OR 3. }

  IF (stopcr <> 2) THEN
  begin
    stpmax:=zero;
    For i:=1 to n do
    begin
      stpmax:=MAX(stpmax, ABS(sn^[i])/MAX(xc^[i], scalex^[i]));
      wv1^[i]:=scalex^[i]*xc^[i]
    end;
    twonrm (overfl, maxexp, n, epsmch, xnorm, wv1);
    IF stpmax <= nsttol*(one+xnorm) THEN
    begin
      trmcod:=1;

{   IF STOPCR=3 THEN OBJECTIVE FUNCTION VALUE MUST BE
    DETERMINED AS WELL - OTHERWISE A SOLUTION HAS BEEN FOUND. }
      IF stopcr <> 3 THEN GOTO 330;

{   NOTE: STATEMENT 202 PRECEDES CONVERGENCE CHECKING SUBROUTINE }
    end
  end;

{   FIND LENGTH OF (SCALED) NEWTON STEP, NEWLEN }

  For i:=1 to n do wv1^[i]:=scalex^[i]*sn^[i];
  twonrm (overfl, maxexp, n, epsmch, newlen, wv1);

{   FOR ITERATIONS AFTER THE ARMIJO STEPS HAVE BEEN COMPLETED
    AT THE BEGINNING (IN OTHER WORDS THE MONOTONIC STEPS)
    STORE THE FUNCTION AND, POSSIBLY, THE NEWTON STEP LENGTHS
    IN THE FTRACK AND STRACK VECTORS, RESPECTIVELY. }

  IF isejac >= narmij THEN
  begin
    IF isejac = 1 THEN
    begin
      strack^[0]:=newlen;

{     NEWMAX IS USED TO KEEP A BOUND ON THE ENTRIES IN
      THE STRACK VECTOR. }

      newmax:=newlen
    end
    ELSE
      strack^[countr]:=MIN(newmax,newlen);

{     THE OBJECTIVE FUNCTION VALUE IS STORED EVEN IF IT IS
      GREATER THAN ANY PRECEEDING FUNCTION VALUE. }

    ftrack^[countr]:=fcnold;

{     WRITE FTRACK AND STRACK VECTORS IF DESIRED.  SINCE ONLY THE LAST
      MGLL VALUES ARE NEEDED THE COUNTER CIRCULATES THROUGH THE VECTOR
      CAUSING ONLY THE MGLL MOST RECENT VALUES TO BE KEPT.
      NOTE: THESE VECTORS ARE NOT APPLICABLE IF NEWTON'S METHOD IS BEING
            USED.  }

    IF (NOT newton) AND (output > 4) THEN
    begin
      Line0; Line0;

{     IF ONLY THE FUNCTION VALUE ACCEPTANCE TEST IS BEING USED,
      THUS ACPTCR=1, THEN ONLY THE FTRACK VECTOR IS APPLICABLE. }

      IF acptcr = 1 THEN
      begin
        Writeln(fp_out,'   *    CURRENT FTRACK VECTOR; LATEST CHANGE: ELEMENT ',countr:4,'                 *');
        Line0;
        For j:=0 to mgll-1 do
          if abs(ftrack^[j]) < 100000 then
            Writeln(fp_out,'   *       FTRACK(',j:3,') = ',ftrack^[j]:11:3,'                                       *')
          else
            Writeln(fp_out,'   *       FTRACK(',j:3,') = ',ftrack^[j]:11:-3,'                                       *');
      end
      ELSE
      begin
{       BOTH THE FUNCTION VALUE AND THE STEP SIZE
        ACCEPTANCE TESTS ARE BEING USED, ACPTCR=12. }

        Writeln(fp_out,'   *    CURRENT FTRACK AND STRACK VECTORS; LATEST CHANGE: ELEMENT',countr:4,'      *');
        Line0;
        For j:=0 to mgll-1 do
          if abs(ftrack^[j]) < 100000 then
            Writeln(fp_out,'   *       FTRACK(',j:3,') = ',ftrack^[j]:11:3,'  STRACK(',j:3,') = ',strack^[j]:11:3,
            '            *')
          else
            Writeln(fp_out,'   *       FTRACK(',j:3,') = ',ftrack^[j]:11:-3,'  STRACK(',j:3,') = ',strack^[j]:11:3,
            '            *')
      end
    end;

{   UPDATE COUNTING INTEGER, COUNTR. RECYCLE IF COUNTR
    HAS REACHED MGLL-1.  }

    IF countr = mgll-1 THEN
      countr:=0
    ELSE
      countr:=countr+1
  end;

{ RESET STEP ACCEPTANCE CODE AND DELSTR }
  acpcod:=0;

  IF (NOT linesr) then delstr:=zero;

{   RESET QNFAIL TO FALSE TO AVOID PREMATURE STOPPING }

  qnfail:=false;

  IF (linesr) THEN
  begin

{   THE MAIN LINE SEARCH IS CALLED IN SUBROUTINE LINE }

    Line (abort, absnew, deuflh, geoms, newton, overch, overfl, qnfail,
          qrsing, restrt, sclfch, sclxch, acpcod, acptcr, contyp,
          isejac, itnum, jupdm, maxexp, maxlin, mgll, mnew, n, narmij,
          nfunc, nunit, output, qnupdm, stopcr, trmcod, alpha, confac,
          epsmch, fcnmax, fcnnew, fcnold, lam0, maxstp, newlen, sbrnrm,
          sigma, a, boundl, boundu, delf, ftrack, fvec, h, hhpi, jac,
          rdiag, wv1, s, sbar, scalef, scalex, sn, strack, xc, xplus);

    IF (abort) THEN   {Line not ok}
    begin
      IF output > 0 THEN
      begin
        Line0; Line1
      end;
      goto RETURN
    end;

{   NOTE: 201 PRECEDES PRINTING OF ITERATION RESULTS }

    IF (newton) then GOTO 320
  end
  ELSE
  begin 
{         TRUST REGION METHOD

          ESTABLISH INITIAL TRUST REGION SIZE, DELTA, AND/OR
          FIND LENGTH OF SCALED DESCENT STEP, CAULEN.  }

    delcau (cauchy, overch, overfl, isejac, maxexp, n, nunit, output,
            beta, caulen, delta, epsmch, maxstp, newlen, sqrtz, a, delf, scalex);

    frstdg:=true;

    IF output > 3 THEN
    begin
      Line0; Line0;
      Msg('    SUMMARY OF TRUST REGION METHOD USING (DOUBLE) DOGLEG STEP.')
    end;

{         MAIN INTERNAL LOOP FOR TRUST REGION METHOD.

          THE TRUST REGION SIZE IS STORED FOR COMPARISON
          LATER TO SET THE PARAMETER POWTAU USED IN POWELL'S
          TRUST REGION UPDATING SCHEME (QUASI-NEWTON, TRUPDM=1). }

    delta0:=delta;

    For notrst:=1 to maxtrs do
    begin

      dogleg (frstdg, newtkn, overch, overfl, maxexp, n, notrst, nunit,
              output, beta, caulen, delta, etafac, newlen, sqrtz, delf,
              s, scalex, sn, ssdhat, vhat);

{       NOTE: WV1 AND WV4 HOLD THE COMPONENT AND RESIDUAL VECTOR
              RESPECTIVELY FOR A TRIAL POINT WHICH HAS BEEN FOUND
              TO BE ACCEPTABLE WHILE THE TRUST REGION IS EXPANDED
              AND A NEW TRIAL POINT TESTED.
              WV2 AND WV3 ARE WORK VECTORS.
              H IS CALLED ASTORE INSIDE TRSTUP. }

      trstup (geoms, newtkn, overch, overfl, qrsing, sclfch, sclxch,
              acpcod, acpstr, acptcr, contyp, isejac, jupdm, maxexp,
              mgll, mnew, n, narmij, nfunc, notrst, nunit, output,
              qnupdm, retcod, trupdm, alpha, confac, delfac, delstr,
              delta, epsmch, fcnmax, fcnnew, fcnold, fcnpre, maxstp,
              newlen, newmax, powtau, rellen, stptol, a, h, boundl,
              boundu, delf, wv1, ftrack, fvec, fvecc, hhpi, jac, rdiag,
              wv2, s, sbar, scalef, scalex, strack, xc, wv4, xplus);

      IF ((output > 4) OR (retcod = 7)) AND (output > 2) then
        rcdprt (nunit, retcod, delta, rellen, stptol);

{          IF NO PROGRESS WAS BEING MADE (RETCOD=7) IN A QUASI-NEWTON
           STEP RETRY WITH AN EXPLICIT JACOBIAN EVALUATION. }

      IF (retcod = 7) AND (NOT restrt) then qnfail:=true;

{          RETURN CODE LESS THAN 8 EXITS FROM TRUST REGION LOOP }

      IF retcod < 8 then GOTO 310

    end;

{       IF NO SUCCESSFUL STEP FOUND IN A QUASI-NEWTON STEP,
        RETRY WITN AN EXPLICIT JACOBIAN EVALUATION. }

    IF (NOT restrt) then qnfail:=true

  end;

310: IF (NOT linesr) THEN
  begin 
    IF delta < delta0 then powtau:=one;
    delta:=MAX(delta, 1E-10)
  end;

{   IF RETCOD=7 AND STOPCR=2, RESET STOPPING CRITERION TO
    AVOID HANGING IN TRUST REGION METHOD. (RETCOD=7 MEANS
    THE RELATIVE STEP LENGTH WAS LESS THAN STPTOL). }

  IF (NOT linesr) AND (retcod = 7) AND (stopcr = 2) AND (NOT qnfail) THEN stopcr:=12;

{   RETAIN NUMBER OF STEPS ACCEPTED BY EACH CRITERION FOR PERFORMANCE EVALUATION. }

  IF (NOT newton) THEN
    IF acpcod = 1 THEN
      nac1:=nac1+1
    ELSE IF acpcod = 2 THEN
      nac2:=nac2+1
    ELSE IF acpcod = 12 THEN
      nac12:=nac12+1;

{ *** PRINT RESULTS OF ITERATION *** }

320:IF output > 2 THEN

    nersl (newton, restrt, sclfch, sclxch, acpcod, jupdm, n, nunit, output, fcnnew, fvec, xplus);

{   CHECK FOR CONVERGENCE.  STATEMENT 202 IS USED IF THE
    STEP SIZE OF THE NEWTON STEP IS FOUND TO BE WITHIN
    THE SPECIFIED TOLERANCE AND STOPCR IS 1 OR 12 }

330:IF (NOT qnfail) THEN
    begin
{   IF QNFAIL IS TRUE THE QUASI-NEWTON SEARCH FAILED TO
    FIND A SATISFACTORY STEP - SINCE THE JACOBIAN IS TO
    BE RE-EVALUATED AVOID PREMATURE STOPPAGES IN NESTOP }

    nestop (absnew, linesr, newton, sclfch, sclxch, acptcr, itnum, n,
            nac1, nac2, nac12, nfunc, njetot, nunit, output, stopcr,
            trmcod, fcnnew, ftol, nsttol, stpmax, stptol, fvec, scalef,
            scalex, xc, xplus)
    end;

{   IF THE TERMINATION CODE, TRMCOD, IS GREATER THAN 0 THEN
    CONVERGENCE HAS BEEN REACHED }

  IF trmcod > 0 then goto RETURN;    {normal exit}

{   QUASI-NEWTON UPDATING - JUPDM > 0 }

  IF jupdm > 0 THEN
  begin

{   QNFAIL MEANS A FAILURE IN THE QUASI-NEWTON SEARCH.
    RE-EVALUATE JACOBIAN AND TRY A DAMPED NEWTON STEP.
    MAXLIN IS CHANGED FROM MAXQNS TO MAXNS OR MAXTRS IS
    CHANGED, SIMILARLY, AT THE START OF THE LOOP }

    IF (qnfail) THEN
      IF output > 4 THEN
      begin
        Line0;
        Msg('       FAILURE IN QUASI-NEWTON SEARCH: QNFAIL IS TRUE.')
      end

    ELSE IF (NOT newton) THEN
    begin

{     QNFAIL IS FALSE - THE STEP HAS BEEN ACCEPTED }

      IF output > 4 THEN
      begin
        Line0;
        Writeln(fp_out,'   *       FCNNEW= ',fcnnew:11:3,'  FCNOLD= ',fcnold:11:3,'  RATIO= ',(fcnnew/fcnold):11:3,'  *')
      end;

      IF fcnnew/fcnold > ratiof THEN
      begin
{     STEP ACCEPTED BUT NOT A SIGNIFICANT IMPROVEMENT }
        nfail:=nfail + 1;
        IF output > 4 THEN
        begin
          Line0;
          Writeln(fp_out,'   *       RATIO > RATIOF SO NFAIL INCREASED TO: ',nfail:5,'                   *');
          Line0
        end
      end
      ELSE
      begin

{       STEP ACCEPTED WITH A SIGNIFICANT IMPROVEMENT }

        IF fcnnew/fcnold > pt01 THEN
          nosubt:=0
        ELSE
          nosubt:=1;

{       ITEMPT IS USED LOCALLY FOR OUTPUT CONTROL }

        itemp:=nfail;
        nfail:=IMAX(nfail-nosubt,0);
        IF output > 4 THEN
        begin
          Line0;
          IF itemp = nfail THEN
            Writeln(fp_out,'   *       NFAIL STAYS AT: ',nfail:5,'                                         *')
          ELSE
            Writeln(fp_out,'   *       NFAIL CHANGED TO: ',nfail:5,'                                       *')
        end
      end;

{     SAVE THE RESULTS FOR RESTART IF A FAILURE IN THE
      QUASI-NEWTON METHOD OCCURS - ESSENTIALLY THIS
      FINDS THE BEST POINT SO FAR }

      IF (isejac = 1) OR (nfail = 1) OR ((nfail <= minqns) AND (fcnnew/fcnmin < one)) THEN
      begin
        savest:=true;
        fcnmin:=fcnnew
      end;

      IF (savest) THEN
      begin
        savest:=false;
        fstore:=fcnnew;
        IF output > 4 THEN
        begin
          Line0;
          Msg('       STEP IS SAVED');
          Line0;
          Msg('       SAVED COMPONENT AND FUNCTION VALUES.');
          Line0
        end;
        For  i:=1 to n do
        begin
          xsave^[i]:=xplus^[i];
          fsave^[i]:=fvec^[i];
          IF output > 4 THEN
            Writeln(fp_out,'   *       XSAVE(',i:3,') = ',xsave^[i]:12:3,'      FSAVE(',i:3,') = ',fsave^[i]:12:3,'      *')
        end;
        IF output > 4 then Line0;
        itstr:=itnum
      end
    end;

{   NOTE: IF QNFAIL IS TRUE THEN NFAIL CANNOT HAVE INCREASED IMPLYING
          THAT NFAIL CANNOT NOW BE GREATER THAN MINQNS }

    IF (qnfail) OR (nfail > minqns) THEN
    begin

{     RESTART FROM BEST POINT FOUND SO FAR }

      restrt:=true;
      IF output > 4 THEN
      begin
        Line0; Line0;
        Msg('       RESTRT IS TRUE.')
      end;
      For j:=0 to mgll-1 do ftrack^[j]:=zero;
      IF acptcr = 12 THEN
        For j:=0 to mgll-1 do strack^[j]:=zero;
      countr:=0;
      isejac:=0;
      mnew:=0;
      trmcod:=0;
      fcnold:=fstore;
      IF output > 4 THEN
      begin
        Line0;
        Writeln(fp_out,'   *       RETURN TO ITERATION:',itstr:5,'   WHERE:                            *');
        Line0
      end;
      For i:=1 to n do
      begin
        xc^[i]:=xsave^[i];
        fvecc^[i]:=fsave^[i];
        IF output > 4 THEN
          Writeln(fp_out,'   *       XC(',i:3,') = ',xc^[i]:12:3,'   FVECC(',i:3,') = ',fvecc^[i]:12:3,'            *')
      end;
      IF output > 4 then Line0
    end
    ELSE
      IF itnum >= niejev then restrt:=false
  end;

{   UPDATE JACOBIAN IF DESIRED:
    QNUPDM = 0  ==> ACTUAL JACOBIAN BEING UPDATED
    QNUPDM = 1  ==> FACTORED JACOBIAN BEING UPDATED }

  IF (NOT restrt) THEN
  begin
    IF qnupdm = 0 THEN
    begin
      IF jupdm = 1 THEN

{     USE BROYDEN UPDATE }

        broyun (overfl, maxexp, n, nunit, output, epsmch, fvec, fvecc, jac, scalex, xc, xplus)

      ELSE IF jupdm = 2 THEN

{     USE LEE AND LEE UPDATE }

        llun (overch, overfl, isejac, maxexp, n, nunit, output, epsmch,
                      omega, fvec, fvecc, jac, plee, s, scalex, xc, xplus);

{     THE FACTORED FORM OF THE JACOBIAN IS UPDATED }

      IF jupdm = 1 THEN

        broyfa (overch, overfl, sclfch, sclxch, maxexp, n, nunit,
                output, epsmch, a, delf, fvec, fvecc, jac, rdiag,
                s, scalef, scalex, wv1, wv2, xc, xplus)

      ELSE IF jupdm = 2 THEN

        llfa (overch, overfl, sclfch, sclxch, isejac, maxexp, n, nunit,
              output, epsmch, omega, a, delf, fvec, fvecc, jac, plee,
              rdiag, s, scalef, scalex, xc, xplus)
    end

  end;

{    UPDATE CURRENT VALUES - RESET TRMCOD TO ZERO.
     UPDATE M "VECTOR" (ACTUALLY ONLY THE LATEST VALUE IS NEEDED). }

  mold:=mnew;
  IF isejac < narmij THEN
    mnew:=0
  ELSE
    mnew:=IMIN(mold+1,mgll-1);

  IF (jupdm = 0) OR ((jupdm > 0) AND (NOT restrt)) THEN

    update (mnew, mold, n, trmcod, fcnnew, fcnold, fvec, fvecc, xc, xplus)

end; {*** itnum main loop ***}

IF output > 0 THEN
begin
  Line1; Line0;
  Writeln(fp_out,'   *  NO SOLUTION FOUND AFTER',(itnum-1):6,' ITERATION(S).                          *');
  Line0;
  Msg('       FINAL ESTIMATES               FINAL FUNCTION VALUES             *');
  Line0;
  For i:=1 to n do
    Writeln(fp_out,'   *    X(',i:3,') = ',xplus^[i]:15:6,'         F(',i:3,') = ',fvec^[i]:15:6,'          *');
  Line0;
  Writeln(fp_out,'   *  FINAL OBJECTIVE FUNCTION VALUE = ',fcnnew:15:6,'                     *');
  Line0; Line1
end;

Return: Dispose(wv1); Dispose(wv2); Dispose(wv4)
End; {nnes}

END.

{end of file unnes2.pas}