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
*                                  Pascal Release 1.1 By J-P Moreau, Paris.  *
*                                              (PART 2/3)                    *
*                                           (www.jpmoreau.fr)                *
*                                                                            *
* Release 1.1: corrected bug in gradf.                                       *
*****************************************************************************}
Unit Unnes1;

Interface

Uses WinCrt, Utils, Fcn, Strings, Unnes;

Var  delfts, lambda, mu, newmax, norm, nrmpre: double;


     Procedure dogleg (var frstdg, newtkn, overch, overfl:boolean; var maxexp, n, notrst, nunit,
                       output:integer; var beta:double; caulen:double; var delta:double; etafac,
                       newlen, sqrtz: double; delf, s, scalex, sn, ssdhat, vhat: pVec);
     Procedure gradf (overch:boolean; var overfl, restrt, sclfch,sclxch:boolean;
                      jupdm, maxexp,n, nunit, output, qnupdm:integer; delf, fvecc:pVec;
                      jac: pMat; scalef, scalex: pVec);
     Procedure initch (var instop, linesr, newton, overfl, sclfch, sclxch:boolean;
                       acptcr, contyp, jactyp, jupdm, maxexp, n, nunit, output, qnupdm,
                       stopcr, trupdm:integer; var epsmch, fcnold:double; ftol:double;
                       boundl, boundu,fvecc, scalef, scalex, xc: pVec);
     Procedure jacobi (checkj:boolean; var jacerr, overfl:boolean; jactyp, n, nunit, output:integer;
                       epsmch:double; var fdtolj:double; boundl, boundu, fvecc, fvecj1, fvecj2:pVec;
                       jac, jacfdm: pMat; scalex, xc: pVec);
     Procedure line (var abort, absnew, deuflh, geoms, newton, overch, overfl,
                     qnfail, qrsing, restrt, sclfch, sclxch:boolean; var acpcod,
                     acptcr, contyp, isejac, itnum, jupdm, maxexp, maxlin, mgll, mnew,
                     n, narmij, nfunc, nunit, output, qnupdm, stopcr, trmcod:integer;
                     var alpha, confac, epsmch, fcnmax, fcnnew, fcnold, lam0,
                     maxstp, newlen, sbrnrm, sigma:double;  a:pMat; boundl, boundu, delf,
                     ftrack, fvec:pVec; h:pMat; hhpi:pVec; jac:pMat; rdiag, rhs, s, sbar,
                     scalef, scalex, sn, strack, xc, xplus: pVec);
     Procedure llfa (var overch, overfl, sclfch, sclxch:boolean; isejac, maxexp, n, nunit,
                     output:integer; var epsmch, omega:double; a:pMat; delf, fvec, fvecc:
                     pVec; jac, plee:pMat; rdiag, s, scalef, scalex, xc, xplus: pVec);
     Procedure llun (overch:boolean; var overfl:boolean; isejac, maxexp, n, nunit,
                     output:integer; epsmch,omega:double; fvec, fvecc:pVec; jac, plee:
                     pMat; s, scalex, xc, xplus: pVec);
     Procedure matprt (nrowpr, ncolpr: integer; a: pMat);
     Procedure maxst (var overfl:boolean; maxexp, n, nunit, output:integer;
                      epsmch:double; var maxstp:double; mstpf:double;
                      scalex, xc: pVec);
     Procedure nersl (newton, restrt, sclfch, sclxch:boolean; acpcod, jupdm,
                      n, nunit, output:integer; fcnnew:double; fvec, xplus:pVec);
     Procedure nestop (var absnew, linesr, newton, sclfch, sclxch:boolean; var
                       acptcr, itnum, n, nac1, nac2, nac12, nfunc, njetot:integer;
                       nunit, output, stopcr:integer; var trmcod:integer; var
                       fcnnew:double; ftol, nsttol:double; var stpmax:double;
                       stptol:double; fvec, scalef, scalex, xc, xplus: pVec);

Implementation

Procedure dogleg (var frstdg, newtkn, overch, overfl:boolean; var maxexp, n, notrst, nunit,
                  output:integer; var beta:double; caulen:double; var delta:double; etafac,
                  newlen, sqrtz: double; delf, s, scalex, sn, ssdhat, vhat: pVec);
{--------------------------------------------------------------------------
!    FEB. 23, 1992
!
!    THIS SUBROUTINE FINDS A TRUST REGION STEP USING THE
!    (DOUBLE) DOGLEG METHOD.
!-------------------------------------------------------------------------}
Label return;
Var
    delfts, eta, factor, gamma, lambda, temp, tempv, zeta: double;
    i:integer;
Begin

  overfl:=false;
  eta := one;

  IF output > 3 THEN
  begin
    Line0;
    Writeln(fp_out,'   *    TRUST REGION STEP:',notrst:6,'  TRUST REGION LENGTH, DELTA:',delta:11:3,'  *');
    Line0;
    Writeln(fp_out,'   *    LENGTH OF NEWTON STEP, NEWLEN: ', newlen:11:3,'                        *')
  end;

{   CHECK FOR NEWTON STEP WITHIN TRUST REGION - IF SO USE NEWTON STEP }

  IF newlen <= delta THEN
  begin
    For i:=1 to n do s^[i]:=sn^[i];
    newtkn:=true;
    temp:=delta;
    delta:=newlen;
    IF output > 3 THEN
    begin
      Line0;
      Writeln(fp_out,'   *       NEWTON STEP WITHIN ACCEPTABLE RANGE ( <= THAN DELTA)           *');
      IF temp = delta THEN
        Writeln(fp_out,'   *       DELTA STAYS AT LENGTH OF NEWTON STEP: ', delta:11:3, '              *')
      ELSE
        Writeln(fp_out,'   *       DELTA SET TO LENGTH OF NEWTON STEP.                            *');
      Line0;
      Writeln(fp_out,'   *       FULL NEWTON STEP ATTEMPTED.                                    *')
    end;
    goto return
  end
  ELSE
  begin
{   NEWTON STEP NOT WITHIN TRUST REGION - APPLY ^[DOUBLE)
    DOGLEG PROCEDURE.  ^[IF ETAFAC EQUALS 1.0 THEN THE SINGLE
    DOGLEG PROCEDURE IS BEING APPLIED).  }

    newtkn:=false;
    IF (frstdg) THEN
    begin
{   SPECIAL SECTION FOR FIRST DOGLEG STEP - CALCULATES
    CAUCHY POINT ^[MINIMIZER OF MODEL FUNCTION IN STEEPEST
    DESCENT DIRECTION OF THE OBJECTIVE FUNCTION).  }

      frstdg:=false;
      IF output > 4 THEN
      begin
        Line0;
        IF etafac = one THEN
          Writeln(fp_out,'   *       FIRST SINGLE DOGLEG STEP                                       *')
        ELSE
          Writeln(fp_out,'   *       FIRST DOUBLE DOGLEG STEP                                       *');
        Line0;
        Writeln(fp_out,'   *          SCALED CAUCHY STEP.                                         *');
        Line0
      end;

{   NOTE: BETA AND SQRTZ WERE CALCULATED IN SUBROUTINE DELCAU }

      zeta:=sqrtz*sqrtz;

{   FIND STEP TO CAUCHY POINT }

      factor:=-(zeta/beta);
      For i:=1 to n do
      begin
        ssdhat^[i]:=factor*(delf^[i]/scalex^[i]);
        IF output > 4 THEN
          Writeln(fp_out,'   *           SSDHAT(',i:3,') = ',ssdhat^[i]:12:3,'                                  *');
      end;
      innerp (overch, overfl, maxexp, n, n, n, nunit, output, delfts, delf, sn);
      overfl:=false;

{   PROTECT AGAINST (RARE) CASE WHEN CALCULATED DIRECTIONAL
    DERIVATIVE EQUALS ZERO.  }

      IF delfts <> zero THEN
      begin
{   STANDARD EXECUTION }
        gamma:=(zeta/ABS(delfts))*(zeta/beta);
        eta:=etafac + (one-etafac)*gamma
      end
      ELSE
      begin
        IF (output > 1) AND (NOT wrnsup) THEN
        begin
          Line0;
          Writeln(fp_out,'   *    WARNING: DELFTS=0,  ETA SET TO 1.0 TO AVOID DIVISION BY ZERO      *')
        end;
        eta:=one
      end;
      IF output > 4 THEN
      begin
        Line0;
        Writeln(fp_out,'   *          ETA = ', eta:11:3, '                                           *');
        Line0;
        Writeln(fp_out,'   *        VHAT VECTOR     VHAT(I) = ETA*SN(I)*SCALEX(I) - SSDHAT(I)      *');
        Line0
      end;
      For i:=1 to n do
      begin
        vhat^[i]:=eta*scalex^[i]*sn^[i] - ssdhat^[i];
        IF output > 4 THEN
          Writeln(fp_out,'   *            VHAT(',i:3,') = ', vhat^[i]:12:3,'                                  *')
      end
    end; {if frstdg...}
  end; {if newlen...}

{ ETA*NEWLEN <= DELTA MEANS TAKE STEP IN NEWTON DIRECTION TO TRUST REGION BOUNDARY }

  IF eta*newlen <= delta THEN
  begin
    IF output > 4 THEN
    begin
      Line0;
      Writeln(fp_out,'   *          ETA*NEWLEN <= DELTA     S(I) = (DELTA/NEWLEN)*SN(I)          *')
    end;
    For i:=1 to n do s^[i]:=(delta/newlen)*sn^[i]
  end
  ELSE
  begin
{   DISTANCE TO CAUCHY POINT >:= DELTA MEANS TAKE STEP IN
    STEEPEST DESCENT DIRECTION TO TRUST REGION BOUNDARY.  }
    IF caulen >= delta THEN
    begin
      IF output > 4 THEN
      begin
        Line0;
        Writeln(fp_out,'   *          CAULEN >= DELTA,  S(I)=(DELTA/CAULEN)*(SSDHAT(I)/SCALEX(I))  *')
      end;
      For i:=1 to n do s^[i] := (delta/caulen)*(ssdhat^[i]/scalex^[i])
    end
    ELSE
    begin
{     TAKE (DOUBLE) DOGLEG STEP }

      innerp (overch, overfl, maxexp, n, n, n, nunit, output, temp, ssdhat, vhat);
      innerp (overch, overfl, maxexp, n, n, n, nunit, output, tempv, vhat, vhat);
      overfl:=false;
      lambda:=(-temp + SQRT(temp*temp - tempv*(caulen*caulen-delta*delta)))/tempv;
      IF output > 4 THEN
      begin
        Line0;
        Writeln(fp_out,'   *          S(I) = (SSDHAT(I)+LAMBDA*VHAT(I))/SCALEX(I)                  *');
        Line0;
        Writeln(fp_out,'   *          WHERE LAMBDA = ', lambda:12:3, '                                 *')
      end;
      For i:=1 to n do s^[i]:=(ssdhat^[i] + lambda*vhat^[i])/scalex^[i]
    end
  end;
  IF output > 3 THEN
  begin
    Line0; Line0;
    Writeln(fp_out,'   *          REVISED STEP FROM SUBROUTINE DOGLEG:                           *');
    Line0;
    For i:=1 to n do
      Writeln(fp_out,'   *           S(',i:3,') = ', s^[i]:12:3,'                                      *')
  end;
return: End; {dogleg}


Procedure gradf (overch:boolean; var overfl, restrt, sclfch,sclxch:boolean;
                jupdm, maxexp, n, nunit, output, qnupdm:integer; delf, fvecc:pVec;
                jac: pMat; scalef, scalex: pVec);
{----------------------------------------------------------------------------
!    FEB. 23, 1992
!
!    THIS SUBROUTINE COMPUTES THE GRADIENT OF THE FUNCTION
!
!           F=1/2{SCALEF*FVECC)^(SCALEF*FVECC)
!
!    WHICH IS USED AS THE OBJECTIVE FUNCTION FOR MINIMIZATION.
!
!    NOTE: WHEN THE FACTORED FORM OF THE JACOBIAN IS UPDATED IN QUASI-NEWTON
!          METHODS THE GRADIENT IS UPDATED AS WELL IN THE SAME SUBROUTINE -
!          IT IS PRINTED HERE THOUGH.  IN THESE CASES QNUPDM > 0.
!---------------------------------------------------------------------------}
Var
    i,j: integer; wv1:pVec;
Begin
  New(wv1);
  IF (restrt) OR (jupdm = 0) OR (qnupdm = 0) THEN
  begin

{   GRADIENT NOT ALREADY UPDATED:    FIND DELF = J^F }
    For i:=1 to n do
      IF (sclfch) THEN
        wv1^[i]:=fvecc^[i]*scalef^[i]*scalef^[i]
      ELSE
        wv1^[i]:=fvecc^[i];

    IF (overch) THEN
{   check each entry individually }
      atvov (overfl, maxexp, n, nunit, output, jac, wv1, delf)
    ELSE
      mmulvt(n,jac,wv1,delf);  {04/05/07 - replaced mmulv by mmulvt}

{   PRINT GRADIENT VECTOR, DELF }

    IF output > 3 THEN
      IF (NOT sclxch) THEN
      begin
        Line0;
        IF (NOT sclfch) THEN
          Writeln(fp_out,'   *    GRADIENT OF OBJECTIVE FUNCTION:                                    *')
        ELSE                                                                            
          Writeln(fp_out,'   *    GRADIENT OF SCALED OBJECTIVE FUNCTION                              *');
        Line0;
        For i:=1 to n do
          if abs(delf^[i]) < 100000.0 then
            Writeln(fp_out,'   *    DELF(',i:3,') = ',delf^[i]:12:3,'                                           *')
          else
            Writeln(fp_out,'   *    DELF(',i:3,') = ',delf^[i]:12:-3,'                                           *')
      end
      ELSE
      begin
        Line0;
        Writeln(fp_out,'   *    GRADIENT OF OBJECTIVE FUNCTION IN SCALED X UNITS                   *');
        Line0;
        For i:=1 to n do
          if abs(scalef^[i]*scalef^[i]*delf^[i]/ scalex^[i]) < 100000.0 then
            Writeln(fp_out,'   *    DELF(',i:3,') = ', delf^[i]:12:3,'         DELF(',i:3,') = ',
                    (scalef^[i]*scalef^[i]*delf^[i]/ scalex^[i]):12:3,'        *')
          else
            Writeln(fp_out,'   *    DELF(',i:3,') = ', delf^[i]:12:3,'         DELF(',i:3,') = ',
                    (scalef^[i]*scalef^[i]*delf^[i]/ scalex^[i]):12:-3,'        *')
      end
  end;

End; {gradf}

Procedure initch (var instop, linesr, newton, overfl, sclfch, sclxch:boolean;
                  acptcr, contyp, jactyp, jupdm, maxexp, n, nunit, output, qnupdm,
                  stopcr, trupdm:integer; var epsmch, fcnold:double; ftol:double;
                  boundl, boundu,fvecc, scalef, scalex, xc: pVec);
{-----------------------------------------------------------------------------
!    AUG. 27, 1991
!
!    THIS SUBROUTINE FIRST CHECKS TO SEE IF N IS WITHIN THE ACCEPTABLE RANGE.
!
!    THE SECOND CHECK IS TO SEE IF THE INITIAL ESTIMATE IS
!    ALREADY A SOLUTION BY THE FUNCTION VALUE CRITERION, FTOL.
!
!    THE THIRD CHECK IS MADE TO SEE IF THE NEWTON OPTION IS BEING
!    USED WITH THE LINE SEARCH.  IF NOT A WARNING IS GIVEN AND
!    THE LINE SEARCH OPTION IS INVOKED.
!
!    THE FOURTH CHECK IS TO ENSURE APPLICABILITY OF SELECTED
!    VALUES FOR INTEGER CONSTANTS.
!
!    THE FIFTH CHECK IS TO WARN THE USER IF INITIAL ESTIMATES ARE NOT WITHIN
!    THE RANGES SET BY THE BOUNDL AND BOUNDU VECTORS.
!    CONTYP IS CHANGED FROM 0 TO 1 IF ANY BOUND HAS BEEN SET BY THE USER
!
!    THE SIXTH CHECK ENSURES BOUNDL^[I) < BOUNDU^[I) FOR ALL I.
!---------------------------------------------------------------------------}
Label 130, return;
Var
  frster:boolean;
  temp1, temp2:double;
  i:integer;
Begin
  instop:=false;
  temp1:=-Power(ten,maxexp);
  temp2:=Power(ten,maxexp);

{ CHECK FOR N IN RANGE }

  IF n <= 0 THEN
  begin
    instop:=true;
    Line1; Line0;
    Writeln(fp_out,'   *  N IS OUT OF RANGE - RESET TO POSITIVE INTEGER                      *');
    Line0; Line1
  end;

{ CHECK FOR SCALING FACTORS POSITIVE }

  frster:=true;
  sclfch:=false;
  sclxch:=false;
  For i:=1 to n do
  begin
    IF scalef^[i] <= zero THEN
    begin
      IF (frster) THEN
      begin
        instop:=true;
        frster:=false;
        Line1
      end;
      Line0;
      Writeln(fp_out,'   *       SCALEF(',i:3,') = ', scalef^[i]:12:3,'    SHOULD BE POSITIVE              *')
    end;
    IF scalef^[i] <> one then sclfch:=true;
    IF scalex^[i] <= zero THEN
    begin
      IF (frster) THEN
      begin
        instop:=true;
        frster:=false;
        Line1
      end;
      Line0;
      Writeln(fp_out,'   *       SCALEX(',i:3,') = ', scalex^[i]:12:3,'    SHOULD BE POSITIVE              *')
    end;
    IF scalex^[i] <> one then sclxch:=true
  end;
  IF (NOT frster) THEN
  begin
    Line0; Line1
  end;

{ EVALUATE INITIAL RESIDUAL VECTOR AND OBJECTIVE FUNCTION AND
  CHECK TO SEE IF THE INITIAL GUESS IS ALREADY A SOLUTION.  }

  fcn1(overfl, n, fvecc, xc);

{ NOTE: NUMBER OF LINE SEARCH FUNCTION EVALUATIONS, NFUNC,
  INITIALIZED AT 1 WHICH REPRESENTS THIS EVALUATION.  }

  IF (overfl) THEN
  begin
    Line1; Line0;
    Writeln(fp_out,'   *       OVERFLOW IN INITIAL FUNCTION VECTOR EVALUATION.               *');
    Line0; Line1;
    instop:=true;
    goto RETURN
  end;

  fcnevl(overfl, maxexp, n, nunit, output, epsmch, fcnold, fvecc, scalef);

  IF (overfl) THEN
  begin
    Line1; Line0;
    Writeln(fp_out,'   *       OVERFLOW IN INITIAL OBJECTIVE FUNCTION EVALUATION.            *');
    Line0; Line1;
    instop:=true;
    goto RETURN
  end;

{ CHECK FOR SOLUTION USING SECOND STOPPING CRITERION }

  For i:=1 to n do
    IF ABS(fvecc^[i]) > ftol Then GOTO 130;

  instop:=true;
  Line1; Line0;
  Writeln(fp_out,'   *  WARNING: THIS IS ALREADY A SOLUTION BY THE CRITERIA OF THE SOLVER  *');
  Line0;

{    IF THE PROBLEM IS BADLY SCALED THE OBJECTIVE FUNCTION MAY MEET THE
     TOLERANCE ALTHOUGH THE INITIAL ESTIMATE IS NOT THE SOLUTION. }

  Writeln(fp_out,'   *  THIS MAY POSSIBLY BE ALLEVIATED BY RESCALING THE PROBLEM IF THE    *');
  Writeln(fp_out,'   *  INITIAL ESTIMATE IS KNOWN NOT TO BE A SOLUTION                     *');
  Line0; Line1;

{ CHECK FOR NEWTON'S METHOD REQUESTED BUT LINE SEARCH NOT BEING USED }

130:IF (newton) AND (NOT linesr) THEN
  begin
    linesr:=true;
    Line1; Line0;
    Writeln(fp_out,'   *  WARNING: INCOMPATIBLE OPTIONS: NEWTON=true AND LINESR=false.       *');
    Writeln(fp_out,'   *  LINESR SET TO true; EXECUTION OF NEWTON METHOD CONTINUING.         *');
    Line0; Line1
  end;

{ CHECK INTEGER CONSTANTS }

  IF (acptcr <> 1) AND (acptcr <> 12) THEN
  begin
    Line1; Line0;
    Writeln(fp_out,'   *  ACPTCR NOT AN ACCEPTABLE VALUE: ',acptcr:5,'                              *');
    instop:=true;
    Line0; Line1
  end;
  IF (jactyp < 0) OR (jactyp > 3) THEN
  begin
    Line1; Line0;
    Writeln(fp_out,'   *  JACTYP:',jactyp:5,' - NOT IN PROPER RANGE.                                *');
    instop:=true;
    Line0; Line1
  end;
  IF (stopcr <> 1) AND (stopcr <> 12) AND (stopcr <> 2) AND (stopcr <> 3) THEN
  begin
    Line1; Line0;
    Writeln(fp_out,'   *  STOPCR NOT AN ACCEPTABLE VALUE: ',stopcr:5,'                              *');
    instop:=true;
    Line0; Line1
  end;
  IF (qnupdm < 0) OR (qnupdm > 1) THEN
  begin
    Line1; Line0;
    Writeln(fp_out,'   *  QNUPDM:',qnupdm:5, ' - NOT IN PROPER RANGE.                               *');
    instop:=true;
    Line0; Line1
  end;
  IF (trupdm < 0) OR (trupdm > 1) THEN
  begin
    Line1; Line0;
    Writeln(fp_out,'   *  TRUPDM:',trupdm:5, ' - NOT IN PROPER RANGE.                               *');
    instop:=true;
    Line0; Line1
  end;
  IF (jupdm < 0) OR (jupdm > 2) THEN
  begin
    Line1; Line0;
    Writeln(fp_out,'   *  JUPDM:',jupdm:5, ' - NOT IN PROPER RANGE.                                 *');
    instop:=true;
    Line0; Line1
  end;

{  CHECK FOR INITIAL ESTIMATES NOT WITHIN SPECIFIED BOUNDS AND
   SET CONTYP TO 1 :=> AT LEAST ONE BOUND IS IN EFFECT.  }

  contyp:=0;
  For i:=1 to n do
    IF (boundl^[i] <> temp1) OR (boundu^[i] <> temp2) THEN
    begin
      contyp:=1;
      goto return
    end;

  frster:=true;
  IF contyp <> 0 THEN
  begin
    For i:=1 to n do
    begin 

{ CHECK FOR INITIAL ESTIMATES OUT OF RANGE AND LOWER
  BOUND GREATER THAN OR EQUAL TO THE UPPER BOUND.  }

      IF (xc^[i] < boundl^[i]) OR (xc^[i] > boundu^[i]) THEN
      begin
        IF (frster) THEN
        begin
          instop:=true;
          frster:=false;
          Line1; Line0;
          Writeln(fp_out,'   *       COMPONENTS MUST BE WITHIN BOUNDS:                             *');
          Line0;
          Writeln(fp_out,'   *      NO         XC               BOUNDL          BOUNDU             *');
          Line0
        end;
        Writeln(fp_out,'   *     ',i:3,'   ', xc^[i]:12:3,'      ',boundl^[i]:12:3,'    ', boundu^[i]:12:3,'            *')
      end
    end;

    IF (NOT frster) THEN
    begin
      Line0; Line1
    end;

    frster:=true;
    For i:=1 to n do
      IF boundl^[i] >= boundu^[i] THEN
      begin
        IF (frster) THEN
        begin
          frster:=false;
          Line1; Line0;
          Writeln(fp_out,'   *       LOWER BOUND MUST BE LESS THAN UPPER BOUND - VIOLATIONS LISTED.              *');
          Line0
        end;
        Writeln(fp_out,'   *       BOUNDL(',i:3,') = ', boundl^[i]:12:3,'    BOUNDU(',i:3,') = ', boundu^[i]:12:3,'      *')
      end;

    IF (NOT frster) THEN
    begin
      Line0; Line1
    end
  end;

return: End; {initch}


Procedure jaccd (n, nunit, output:integer; epsmch:double; fvecj1, fvecj2:pVec;
                 jacfdm: pMat; scalex, xc: pVec);
{----------------------------------------------------------------------
!    FEB. 11, 1991
!
!    THIS SUBROUTINE EVALUATES THE JACOBIAN USING CENTRAL DIFFERENCES.
!
!    FVECJ1 AND FVECJ2 ARE TEMPORARY VECTORS TO HOLD THE
!    RESIDUAL VECTORS FOR THE CENTRAL DIFFERENCE CALCULATION.
!---------------------------------------------------------------------}
Var
  overfl: boolean;
  j,k: integer;
  curtep, deltaj, tempj: double;
Begin
  overfl:=false;
  curtep:=Power1(epsmch, 0.33);

  For j:=1 to n do
  begin
    deltaj:=curtep*SIGN((MAX(ABS(xc^[j]),one/scalex^[j])),xc^[j]);
    tempj:=xc^[j];
    xc^[j]:=xc^[j] + deltaj;

{   NOTE: THIS STEP IS FOR FLOATING POINT ACCURACY ONLY }

    deltaj:=xc^[j] - tempj;

    fcn1(overfl, n, fvecj1, xc);
    xc^[j]:=tempj-deltaj;
    fcn1(overfl, n, fvecj2, xc);
    IF (overfl) AND (output > 2) AND (NOT wrnsup) THEN
    begin
      Line0;
      Writeln(fp_out,'   *    WARNING: OVERFLOW IN FUNCTION VECTOR IN "JACCD".                 *')
    end;
    For k:=1 to n do
      jacfdm^[k,j]:=(fvecj1^[k] - fvecj2^[k])/(two*deltaj);
    xc^[j]:=tempj
  end;

End; {jaccd}

Procedure jacfd (jactyp, n, nunit, output:integer; epsmch:double;
                  boundl, boundu, fvecc, fvecj1:pVec; jacfdm:pMat;
                  scalex, xc: pVec);
{-----------------------------------------------------------------------------
!   FEB. 15, 1991
!
!   THIS SUBROUTINE EVALUATES THE JACOBIAN USING ONE-SIDED FINITE DIFFERENCES.
!
!   JACTYP "1" SIGNIFIES FORWARD DIFFERENCES
!   JACTYP "2" SIGNIFIES BACKWARD DIFFERENCES
!
!   FVECJ1 IS A TEMPORARY VECTOR WHICH STORES THE RESIDUAL
!   VECTOR FOR THE FINITE DIFFERENCE CALCULATION.
!----------------------------------------------------------------------------}
Var
    deltaj, sqrtep, tempj: double;
    overfl: boolean;
    j,k: integer;
Begin

  sqrtep:=SQRT(epsmch);

{   FINITE-DIFFERENCE CALCULATION BY COLUMNS }

  For j:=1 to n do
  begin

{   DELTAJ IS THE STEP SIZE - IT IS ALWAYS POSITIVE }

    deltaj:=sqrtep*MAX(ABS(xc^[j]), one/scalex^[j]);
    tempj:=xc^[j];
{   temporary storage of xc(j)  }
    IF jactyp = 1 THEN
    begin
      IF xc^[j]+deltaj <= boundu^[j] THEN
      begin
{     STEP WITHIN BOUNDS - COMPLETE FORWARD DIFFERENCE }
        xc^[j]:=xc^[j] + deltaj;
        deltaj:=xc^[j] - tempj;
        fordif (overfl, j, n, deltaj, fvecc, fvecj1, jacfdm, xc)
      end
      ELSE
      begin
{     STEP WOULD VIOLATE BOUNDU - TRY BACKWARD DIFFERENCE }

        IF xc^[j]-deltaj >= boundl^[j] THEN
        begin
          xc^[j]:=xc^[j] - deltaj;
          bakdif (overfl, j, n, deltaj, tempj, fvecc, fvecj1, jacfdm, xc);
        end
        ELSE
        begin
{       STEP WOULD ALSO VIOLATE BOUNDL - IF THE DIFFERENCE IN THE
        BOUNDS, (BOUNDU-BOUNDL), IS GREATER THAN DELTAJ CALCULATE
        THE FUNCTION VECTOR AT EACH BOUND AND USE THIS DIFFERENCE -
        THIS REQUIRES ONE EXTRA FUNCTION EVALUATION.
        THE CURRENT FVECC IS STORED IN WV3, THEN REPLACED.  }
          IF boundu^[j]-boundl^[j] >= deltaj THEN
          begin
            bnddif (overfl, j, n, epsmch, boundl, boundu, fvecc, fvecj1, jacfdm, xc);
            IF (output > 2) AND (NOT wrnsup) AND (NOT overfl) THEN
            begin
              Line0;
              Writeln(fp_out,'   *    WARNING: BOUNDS TOO CLOSE FOR 1-SIDED FINITE-DIFFERENCES.        *');
              Writeln(fp_out,'   *           LOWER AND UPPER BOUNDS USED FOR JACOBIAN COLUMN: ',j:3,'      *');
              Writeln(fp_out,'   *           THIS REQUIRED ONE EXTRA FUNCTION EVALUATION.              *')
            end
          end
          ELSE
          begin
{         BOUNDS ARE EXTREMELY CLOSE ^[BUT NOT EQUAL OR
          THE PROGRAM WOULD HAVE STOPPED IN INITCH].     }

            IF (output > 2) AND (NOT wrnsup) AND (NOT overfl) THEN
            begin
              Line0;
              Writeln(fp_out,'   *    WARNING: BOUNDS TOO CLOSE FOR 1-SIDED FINITE-DIFFERENCES.        *');
              Writeln(fp_out,'   *           BOUNDS ARE EXTREMELY CLOSE FOR COMPONENT: ',j:3,'             *');
              Writeln(fp_out,'   *           FINITE DIFFERENCE JACOBIAN IS UNRELIABLE.                 *')
            end;
            bnddif (overfl, j, n, epsmch, boundl, boundu, fvecc, fvecj1, jacfdm, xc);
          end
        end
      end
    end
    ELSE  {jactyp<>1}
    begin
      IF xc^[j]-deltaj >= boundl^[j] THEN
      begin
        xc^[j]:=xc^[j] - deltaj;
        bakdif (overfl, j, n, deltaj, tempj, fvecc, fvecj1, jacfdm, xc)
      end
      ELSE
      begin
        IF xc^[j]+deltaj <= boundu^[j] THEN
        begin
          xc^[j]:=xc^[j] + deltaj;
          fordif (overfl, j, n, deltaj, fvecc, fvecj1, jacfdm, xc)
        end
        ELSE
        begin
          IF boundu^[j]-boundl^[j] >= deltaj THEN
          begin
            bnddif (overfl, j, n, epsmch, boundl, boundu, fvecc, fvecj1, jacfdm, xc);
            IF (output > 2) AND (NOT wrnsup) AND (NOT overfl) THEN
            begin
              Line0;
              Writeln(fp_out,'   *    WARNING: BOUNDS TOO CLOSE FOR 1-SIDED FINITE-DIFFERENCES.        *');
              Writeln(fp_out,'   *           LOWER AND UPPER BOUNDS USED FOR JACOBIAN COLUMN: ',j:3,'      *');
              Writeln(fp_out,'   *           THIS REQUIRED ONE EXTRA FUNCTION EVALUATION.              *')
            end
          end
          ELSE
          begin
            bnddif (overfl, j, n, epsmch, boundl, boundu, fvecc, fvecj1, jacfdm, xc);
            For k:=1 to n do jacfdm^[k,j]:=zero;
            IF (output > 2) AND (NOT wrnsup) AND (NOT overfl) THEN
            begin
              Line0;
              Writeln(fp_out,'   *    WARNING: BOUNDS TOO CLOSE FOR 1-SIDED FINITE-DIFFERENCES.        *');
              Writeln(fp_out,'   *           BOUNDS ARE EXTREMELY CLOSE FOR COMPONENT: ',j:3,'             *');
              Writeln(fp_out,'   *           FINITE DIFFERENCE JACOBIAN IS UNRELIABLE.                 *')
            end
          end
        end
      end
    end;
    IF (overfl) AND (output > 2) AND (NOT wrnsup) THEN
    begin
      Line0;
      Writeln(fp_out,'   *    WARNING: OVERFLOW IN FUNCTION VECTOR IN SUBROUTINE JACFD.        *');
      Writeln(fp_out,'   *           BOUNDS ARE EXTREMELY CLOSE FOR COMPONENT: ',j:3,'             *');
      For k:=1 to n do jacfdm^[k,j]:=zero
    end;
    xc^[j]:=tempj
  end
 
End; {JACFD}

Procedure jacobi (checkj:boolean; var jacerr, overfl:boolean; jactyp, n, nunit, output:integer;
                  epsmch:double; var fdtolj:double; boundl, boundu, fvecc, fvecj1, fvecj2:pVec;
                  jac, jacfdm: pMat; scalex, xc: pVec);
{-------------------------------------------------------------------------
!    APR. 13, 1991
!
!    THIS SUBROUTINE EVALUATES THE JACOBIAN.  IF CHECKJ IS TRUE
!    THEN THE ANALYTICAL JACOBIAN IS CHECKED NUMERICALLY.
!
!    jacob IS A USER-SUPPLIED ANALYTICAL JACOBIAN USED ONLY IF JACTYP=0.
!    THE JACOBIAN NAME MAY BE CHANGED BY USING THE EXTERNAL STATEMENT
!    IN THE MAIN DRIVER.
!
!    JACFD ESTIMATES THE JACOBIAN USING FINITE DIFFERENCES:
!    FORWARD IF JACTYP:=1 OR BACKWARD IF JACTYP:=2.
!
!    JACCD ESTIMATES THE JACOBIAN USING CENTRAL DIFFERENCES.
!
!    IF THE ANALYTICAL JACOBIAN IS CHECKED THE FINITE DIFFERENCE
!    JACOBIAN IS STORED IN "JACFDM" AND THEN COMPARED.
!
!    FRSTER  INDICATES FIRST ERROR - USED ONLY TO SET BORDERS FOR OUTPUT
!    JACERR  FLAG TO INDICATE TO THE CALLING PROGRAM AN ERROR
!            IN THE ANALYTICAL JACOBIAN
!------------------------------------------------------------------------}
Var
    frster: boolean;
    i, j: integer;
Begin
  frster:=true;
  jacerr:=false;
  overfl:=false;

  IF jactyp = 0 THEN
    jacob (overfl, n, jac, xc)
  ELSE IF (jactyp = 1) OR (jactyp = 2) THEN
    jacfd (jactyp, n, nunit, output, epsmch, boundl, boundu, fvecc, fvecj1, jac, scalex, xc)
  ELSE
    jaccd (n, nunit, output, epsmch, fvecj1, fvecj2, jac, scalex, xc);

  IF (jactyp = 0) AND (checkj) THEN
  begin

{   NOTE: JACTYP=0 SENT TO JACFD PRODUCES A FORWARD DIFFERENCE ESTIMATE OF THE JACOBIAN.  }

    jacfd (jactyp, n, nunit, output, epsmch, boundl, boundu, fvecc, fvecj1, jacfdm, scalex, xc);

    For j:=1 to n do
      For i:=1 to n do
      begin
        IF ABS((jacfdm^[i,j]-jac^[i,j])/MAX(MAX(ABS(jac^[i,j]),ABS(jacfdm^[i,j])),one)) > fdtolj THEN
        begin
          jacerr:=true;
          IF output >= 0 THEN
          begin
            IF (frster) THEN
            begin
              frster:=false;
              Line1
            end;
            Line0;
            Writeln(fp_out,'   *    CHECK JACOBIAN TERM (',i:3,',',j:3,'):                                 *');
            Line0;
            Writeln(fp_out,'   *    ANALYTICAL DERIVATIVE IS ',jac^[i,j]:12:3,'                            *');
            Writeln(fp_out,'   *     NUMERICAL DERIVATIVE IS ',jacfdm^[i,j]:12:3,'                            *');
            Line0; Line1
          end
        end
      end
  end

End; {jacobi}

Procedure line (var abort, absnew, deuflh, geoms, newton, overch, overfl,
                qnfail, qrsing, restrt, sclfch, sclxch:boolean; var acpcod,
                acptcr, contyp, isejac, itnum, jupdm, maxexp, maxlin, mgll, mnew,
                n, narmij, nfunc, nunit, output, qnupdm, stopcr, trmcod:integer;
                var alpha, confac, epsmch, fcnmax, fcnnew, fcnold, lam0,
                maxstp, newlen, sbrnrm, sigma:double;  a:pMat; boundl, boundu, delf,
                ftrack, fvec:pVec; h:pMat; hhpi:pVec; jac:pMat; rdiag, rhs, s, sbar,
                scalef, scalex, sn, strack, xc, xplus: pVec);
{--------------------------------------------------------------
!    SEPT. 9, 1991
!-------------------------------------------------------------}
Label 90, return;
Var
    {delfts, lambda, mu, newmax, norm, nrmpre: double; }
    convio: boolean;
    i,j,k : integer;
    wv2   : pVec;
Begin
  New(wv2);
  convio:=false;
  overfl:=false;
  abort := FALSE;
  qnfail := FALSE;

  IF (newton) OR (absnew) THEN
  begin
    IF (newton) THEN

{   FIND NEXT ITERATE FOR PURE NEWTON'S METHOD }
      For k:=1 to n do xplus^[k]:=xc^[k] + sn^[k]
    ELSE

{   FIND NEXT ITERATE FOR "ABSOLUTE" NEWTON'S METHOD.
    IF COMPONENT I WOULD BE OUTSIDE ITS BOUND THEN TAKE ABSOLUTE VALUE OF
    THE VIOLATION AND GO THIS DISTANCE INTO THE FEASIBLE REGION.
    ENSURE THAT THIS REFLECTION OFF ONE BOUND DOES NOT VIOLATE THE OTHER }

      For i:=1 to n do
      begin
        wv2^[i]:=zero;
        IF sn^[i] >= zero THEN
          IF xc^[i]+sn^[i] > boundu^[i] THEN
          begin
            convio:=true;
            wv2^[i]:=two;
            xplus^[i]:=MAX(two*boundu^[i]-xc^[i]-sn^[i], boundl^[i])
          end
          ELSE
            xplus^[i]:=xc^[i] + sn^[i]
        ELSE
          IF xc^[i]+sn^[i] < boundl^[i] THEN
          begin
            convio:=true;
            wv2^[i]:=-two;
            xplus^[i]:=MIN(two*boundl^[i]-xc^[i]-sn^[i], boundu^[i])
          end
          ELSE
            xplus^[i]:=xc^[i] + sn^[i]
      end;

    IF (convio) AND (output > 4) THEN
    begin
      Line0;
      Msg('       CONSTRAINT VIOLATORS IN ABSOLUTE NEWTON''S METHOD');
      Line0;
      Msg('       COMPONENT  PROPOSED POINT  VIOLATED BOUND  FEASIBLE VALUE');
      Line0;
      For i:=1 to n do
      begin
        IF wv2^[i] > zero THEN
          Writeln(fp_out,'   *     ',i:6,'     ',(xc^[i]+sn^[i]):12:3,'     ',boundu^[i]:12:3,'    ',
                  xplus^[i]:12:3,'        *')
        ELSE IF wv2^[i] < zero THEN
          Writeln(fp_out,'   *     ',i:6,'     ',(xc^[i]+sn^[i]):12:3,'     ',boundl^[i]:12:3,'    ',
                  xplus^[i]:12:3,'        *')
      end
    end;

    fcn1(overfl, n, fvec, xplus);

    nfunc:=nfunc + 1;
    IF (overfl) THEN
    begin
      overfl:=false;
      fcnnew:=Power(ten, maxexp);
      IF output > 2 THEN
      begin
        Line0;
        Msg('       POTENTIAL OVERFLOW DETECTED IN FUNCTION EVALUATION.')
      end;
      GOTO 90
    end;

    fcnevl (overfl, maxexp, n, nunit, output, epsmch, fcnnew, fvec, scalef);

{    RETURN FROM PURE NEWTON'S METHOD - OTHERWISE CONDUCT LINE SEARCH }

90: goto RETURN
  end;

  IF output > 3 THEN
  begin
    Line0; Line0;
    Msg('    SUMMARY OF LINE SEARCH');
    Line0
  end;

{    SHORTEN NEWTON STEP IF LENGTH IS GREATER THAN MAXSTP }

  IF newlen > maxstp THEN
  begin
    IF output > 3 THEN
    begin
      Line0;
      Msg('       LENGTH OF NEWTON STEP SHORTENED TO MAXSTP.')
    end;
    For k:=1 to n do sn^[k]:=sn^[k]*maxstp/newlen
  end;

{    CHECK DIRECTIONAL DERIVATIVE (MAGNITUDE AND SIGN). }

  innerp (overch, overfl, maxexp, n, n, n, nunit, output, delfts, delf, sn);

  IF (overfl) THEN
  begin
    overfl:=false;
    IF (output > 2) AND (NOT wrnsup) THEN
    begin
      Line0;
      Writeln(fp_out,'   *    WARNING: DIRECTIONAL DERIVATIVE, DELFTS, SET TO ',delfts:12:3,'     *')
    end
  end;

{   REVERSE SEARCH DIRECTION IF DIRECTIONAL DERIVATIVE IS POSITIVE }

  IF delfts > zero THEN
  begin
    For k:=1 to n do sn^[k]:=-sn^[k];
    delfts:=-delfts;
    IF (output > 2) AND (NOT wrnsup) THEN
    begin
      Line0;
      Msg('    WARNING: DIRECTIONAL DERIVATIVE IS POSITIVE: DIRECTION REVERSED.')
    end
  end;

{    OUTPUT INFORMATION }

  IF (output > 3) THEN
  begin
    if Abs(delfts) < 100000 then
      Writeln(fp_out,'   *       INNER PRODUCT OF DELF AND SN, DELFTS: ....... ',delfts:12:3,'      *')
    else
      Writeln(fp_out,'   *       INNER PRODUCT OF DELF AND SN, DELFTS: ....... ',delfts:12:-5,'      *');
    IF (NOT sclxch) THEN
      Writeln(fp_out,'   *       LENGTH OF NEWTON STEP, NEWLEN: .............. ',newlen:12:3,'      *')
    ELSE
      Writeln(fp_out,'   *       LENGTH OF SCALED NEWTON STEP, NEWLEN: ....... ',newlen:12:3,'      *');
    Writeln(fp_out,'   *       MAXIMUM STEP LENGTH ALLOWED, MAXSTP: ........ ',maxstp:12:3,'      *')
  end;

{    ESTABLISH INITIAL RELAXATION FACTOR }

  IF (deuflh) THEN

{    AT FIRST STEP IN DAMPED NEWTON OR AFTER EXPLICIT JACOBIAN EVALUATION
     IN QUASI-NEWTON OR IF THE STEP SIZE IS WITHIN STOPPING TOLERANCE BUT
     STOPCR=3, THE LINE SEARCH IS STARTED AT LAMBDA=1.  }

    IF ((isejac = 1) OR (trmcod = 1)) AND (stopcr = 3) THEN
      lambda:=lam0
    ELSE
    begin
      For k:=1 to n do wv2^[k]:=(sbar^[k] - sn^[k])*scalex^[k];

      twonrm (overfl, maxexp, n, epsmch, norm, wv2);

{     PREVENT DIVIDE BY ZERO IF NORM IS ZERO (UNDERFLOWS) }

      IF norm < epsmch THEN
{     START LINE SEARCH AT LAMBDA=LAM0, USE DUMMY MU }
        mu:=ten
      ELSE
        mu:=nrmpre*lambda/norm;

      IF output > 4 THEN
      begin
        Line0;
        Writeln(fp_out,'   *       DEUFLHARD TEST RATIO, MU: ',mu:11:3,'                         *')
      end;

{     SET INITIAL LAMBDA DEPENDING ON MU.   THIS IS A MODIFICATION OF
      DEUFLHARD'S METHOD WHERE THE CUTOFF VALUE WOULD BE 0.7 FOR LAM0=1.0 }

      IF mu > lam0/ten THEN
        lambda:=lam0
      ELSE
        lambda:=lam0/ten
    end
  ELSE
    lambda:=lam0;

{   STORE LENGTH OF NEWTON STEP.  IF NEWTON STEP LENGTH WAS
    GREATER THAN MAXSTP IT WAS SHORTENED TO MAXSTP. }

  nrmpre:=MIN(maxstp,newlen);

{   ESTABLISH FCNMAX AND NEWMAX FOR NONMONOTONIC LINE SEARCH }

  newmax:=newlen;
  fcnmax:=fcnold;
  IF isejac > narmij THEN
    IF isejac < narmij+mgll THEN
      For j:=1 to mnew do
      begin
        fcnmax:=MAX(fcnmax,ftrack^[j-1]);
        newmax:=MAX(newmax,strack^[j-1])
      end
    ELSE
      For j:=0 to mnew do
      begin
        fcnmax:=MAX(fcnmax, ftrack^[j]);
        newmax:=MAX(newmax, strack^[j])
      end;

  IF output > 3 THEN
  begin
    Line0; Line0;
    IF (NOT sclxch) THEN
      Msg('       LINE SEARCH')
    ELSE
      Msg('       LINE SEARCH (X''S GIVEN IN UNSCALED UNITS)')
  end;

{   CONDUCT LINE SEARCH }

  deufls (abort, deuflh, geoms, overch, overfl, qnfail, qrsing, restrt,
          sclfch, sclxch, acpcod, acptcr, contyp, itnum, jupdm, maxexp,
          maxlin, n, nfunc, nunit, output, qnupdm, stopcr, alpha, confac,
          delfts, epsmch, fcnmax, fcnnew, fcnold, lambda, newmax, sbrnrm,
          sigma, a, h, boundl, boundu, delf, fvec, hhpi, jac, rdiag, rhs,
          s, sbar, scalef, scalex, sn, xc, xplus);

return: Dispose(wv2)

End;  {line}

Procedure llfa (var overch, overfl, sclfch, sclxch:boolean; isejac, maxexp, n, nunit,
                output:integer; var epsmch, omega:double; a:pMat; delf, fvec, fvecc:
                pVec; jac, plee:pMat; rdiag, s, scalef, scalex, xc, xplus: pVec);
{---------------------------------------------------------
!    FEB. 23, 1991
!
!    THE LEE AND LEE QUASI-NEWTON METHOD IS APPLIED TO
!    THE FACTORED FORM OF THE JACOBIAN.
!
!    NOTE: T AND W ARE TEMPORARY WORKING VECTORS ONLY.
!--------------------------------------------------------}
Label return;
Var
  denom, sqrtep, SUM:double;
  t, w, wv3: pVec; temp1, temp2: pMat;
  skipup:boolean;
  i,j,k:integer;
Begin
  New(t); New(w); New(wv3); New(temp1); New(temp2);
  overfl:=false;
  sqrtep:=SQRT(epsmch);
  skipup:=true;

  For i:=1 to n do
  begin
    a^[i,i]:=rdiag^[i];
    s^[i]:=xplus^[i] - xc^[i]
  end;

{  R IS IN THE UPPER TRIANGLE OF A.

        T=RS   }

  uvmul(n, n, n, n, a, s, t);

{  FORM PART OF NUMERATOR AND CHECK TO SEE IF A SIGNIFICANT
   CHANGE WOULD BE MADE TO THE JACOBIAN.  }

  For i:=1 to n do
  begin
    For k:=1 to n do w^[k]:=jac^[1,k];
    innerp (overch, overfl, maxexp, n, n, n, nunit, output, sum, w, t);
    w^[i]:=scalef^[i]*(fvec^[i] - fvecc^[i]) - sum;

{  TEST TO ENSURE VECTOR W IS NONZERO.  IF W^[I]:=0 FOR
   ALL I THEN THE UPDATE IS SKIPPED - SKIPUP IS TRUE. }

    IF ABS(w^[i]) > sqrtep*scalef^[i]*(ABS(fvec^[i]) + ABS(fvecc^[i])) THEN
      skipup:=false
{  update to be performed }
    ELSE
      w^[i]:=zero
  end;

  IF (NOT skipup) THEN
  begin
{   T:=Q^W   Q^ IS STORED IN JAC }

    avmul (n, n, n, n, jac, w, t);

{   FIND DENOMINATOR; FORM W:=S^P (P IS SYMMETRIC SO PS IS FOUND) }

    avmul (n, n, n, n, plee, s, w);
    IF (sclxch) THEN
{   SCALE W TO FIND DENOMINATOR }
      For k:=1 to n do wv3^[k]:=w^[k]*scalex^[k]*scalex^[k]
    ELSE
    begin
      For k:= 1 to n do
      begin
        temp1^[k,1]:=w^[k];
        temp2^[k,1]:=wv3^[k]
      end;
      matcop (n, n, 1, 1, n, 1, temp1, temp2);
      For k:=1 to n do wv3^[k]:=temp2^[k,1]
    end;

    innerp (overch, overfl, maxexp, n, n, n, nunit, output, denom, wv3, s);

{   IF OVERFLOW WOULD OCCUR MAKE NO CHANGE TO JACOBIAN }

    IF (overfl) THEN
    begin
      IF output > 3 THEN
      begin
        Line0;
        Msg('    WARNING: JACOBIAN NOT UPDATED TO AVOID OVERFLOW IN DENOMINATOR OF');
        Msg('    LEE AND LEE UPDATE.')
      end;
      goto RETURN
    end;

{   IF DENOM IS ZERO THE SOLVER IS PROBABLY NEAR SOLUTION -
    AVOID OVERFLOW AND CONTINUE WITH SAME JACOBIAN. }

    IF denom = zero then goto RETURN;

{   THE SCALED VERSION OF S IS TAKEN TO THE UPDATE }

    For k:=1 to n do w^[k]:=w^[k]*scalex^[k]*scalex^[k]/denom;

{   UPDATE THE QR DECOMPOSITION USING A SERIES OF GIVENS (JACOBI) ROTATIONS }

    qrupda (overfl, maxexp, n, epsmch, a, jac, t, w);

{   RESET RDIAG AS DIAGONAL OF CURRENT R WHICH IS IN THE UPPER TRIANGLE OF A }

    For i:=1 to n do rdiag^[i]:=a^[i,i];

{   UPDATE P MATRIX }

    denom:=Power(omega,(isejac+2)) + denom;
    plee^[1,1]:=plee^[1,1] - wv3^[1]*wv3^[1]/denom;
    For j:=2 to n do
    begin
      For i:=1 to j-1 do
      begin
        plee^[i,j]:=plee^[i,j] - wv3^[i]*wv3^[j]/denom;
        plee^[j,i]:=plee^[i,j]
      end;
      plee^[j,j]:=plee^[j,j] - wv3^[j]*wv3^[j]/denom
    end
  end;

{   UPDATE THE GRADIENT VECTOR, DELF.
    DELF = (QR)^F = R^Q^F = R^JAC F }

  IF (sclfch) THEN
    For k:=1 to n do w^[k]:=fvec^[k]*scalef^[k]
  ELSE
  begin
    For k:= 1 to n do
    begin
      temp1^[k,1]:=fvec^[k];
      temp2^[k,1]:=w^[k]
    end;
    matcop (n, n, 1, 1, n, 1, temp1, temp2);
    For k:=1 to n do w^[k]:=temp2^[k,1]
  end;
  avmul (n, n, n, n, jac, w, t);
  mmulv (n, a, t, delf);

return:Dispose(t); Dispose(w); Dispose(wv3); Dispose(temp1); Dispose(temp2)
End; {llfa}

Procedure llun (overch:boolean; var overfl:boolean; isejac, maxexp, n, nunit,
                output:integer; epsmch,omega:double; fvec, fvecc:pVec; jac, plee:
                pMat; s, scalex, xc, xplus: pVec);
{--------------------------------------------------------
!
!    FEB. 13, 1991
!
!    UPDATE THE JACOBIAN USING THE LEE AND LEE METHOD.
!-------------------------------------------------------}
Label return;
Var
  denom, sqrtep, SUM, tempi:double;
  wv1,temp,temp1: pVec;
  i, j:integer;
Begin
  New(wv1); New(temp); New(temp1);
  sqrtep:=SQRT(epsmch);
  For i:=1 to n do s^[i]:=(xplus^[i] - xc^[i])*scalex^[i];
  For i:=1 to n do
  begin
    For j:=1 to n do temp^[j]:=plee^[j,i]; 
    wv1^[i] := DOT_PRODUCT(n, s, temp)
  end;

  innerp (overch, overfl, maxexp, n, n, n, nunit, output, denom, wv1, s);

{ IF OVERFLOW WOULD OCCUR MAKE NO CHANGE TO JACOBIAN }

  IF (overfl) THEN
  begin
    IF output > 3 THEN
    begin
      Line0;
      Writeln(fp_out,'    WARNING: JACOBIAN NOT UPDATED TO AVOID OVERFLOW IN DENOMINATOR OF');
      Writeln(fp_out,'    LEE AND LEE UPDATE.')
    end;
    goto RETURN
  end;

{  IF DENOM IS ZERO THE SOLVER MUST BE VERY NEAR SOLUTION -
   AVOID OVERFLOW AND CONTINUE WITH SAME JACOBIAN  }

  IF denom = zero then goto RETURN;
  For i:=1 to n do
  begin
    For j:= 1 to n do
    begin
      temp^[i]:=jac^[i,j];
      temp1^[i]:=xplus^[j] - xc^[j]
    end;
    sum := DOT_PRODUCT(n, temp, temp1);
    tempi:=fvec^[i] - fvecc^[i] - sum;
    IF ABS(tempi) >= sqrtep*(ABS(fvec^[i]) + ABS(fvecc^[i])) THEN
    begin
      tempi:=tempi/denom;
      For j:=1 to n do jac^[i,j]:=jac^[i,j] + tempi*wv1^[j]*scalex^[j]
    end
  end;

{ UPDATE P MATRIX }

  denom:=Power(omega,isejac+2) + denom;
  plee^[1,1]:=plee^[1,1] - wv1^[1]*wv1^[1]/denom;
  For j:=2 to n do
  begin
    For i:=1 to j-1 do
    begin
      plee^[i,j]:=plee^[i,j] - wv1^[i]*wv1^[j]/(denom*scalex^[i]*scalex^[j]);
      plee^[j,i]:=plee^[i,j]
    end;
    plee^[j,j]:=plee^[j,j] - wv1^[j]*wv1^[j]/denom
  end;

return: Dispose(wv1); Dispose(temp); Dispose(temp1)
End; {llun}

Procedure matprt (nrowpr, ncolpr: integer; a: pMat);
{-----------------------------------------------------------------------------
!    FEB. 6, 1991
!
!    THIS SUBROUTINE PRINTS RECTANGULAR BLOCKS STARTING WITH ELEMENT A^[1,1]
!    OF SIZE NROWPR BY NCOLPR FOR MATRIX A ^[WHICH HAS DECLARED SIZE NROWA BY
!    NCOLA].  THE MATRIX IS PRINTED AS A BLOCK FOR SIZES UP TO 5X5 OR BY
!    COLUMNS IF IT IS LARGER.
!
!    NROWPR IS THE NUMBER OF ROWS TO BE PRINTED
!    NCOLPR IS THE NUMBER OF COLUMNS TO BE PRINTED
!
!    IF MATRIX PRINTING IS TO BE SUPPRESSED THEN LOGICAL
!    VARIABLE MATSUP MUST BE SET TO TRUE BEFORE THE CALL TO NNES.
!----------------------------------------------------------------------------}
Label return;
Var
  i, j, k, limit: integer;
  s:array[0..79] of char; s1:array[0..20] of char; s2:array[0..3] of char;
Begin

  IF (matsup) THEN
  begin
    Line0;
    Msg('       MATRIX PRINTING SUPPRESSED.');
    goto return
  end;

{ FOR NCOLPR <= 5 WRITE MATRIX AS A WHOLE }

  strcopy(s,''); strcopy(s1,''); strcopy(s2,'');
  Line0;
  IF ncolpr <= 5 THEN
  begin
    Strcopy(s,'           ');
    For k:= 1 to ncolpr do
    begin
      Str(k:3,s2); StrCat(s,s2); StrCat(s,'        ');
    end;
    Msg(s); Line0;
    For i:=1 to nrowpr do
    begin
      strcopy(s,'');
      Str(i:3,s2); StrCat(s,s2); StrCat(s,'     ');
      For k:= 1 to ncolpr do
      begin
        Str(a^[i,k]:10:3,s1); StrCat(s,s1); StrCat(s,' ')
      end;
      Msg(s);
    end;
    Line0
  end
  ELSE
  begin

{   LIMIT IS THE NUMBER OF GROUPS OF 5 COLUMNS }
    limit:=ncolpr Div 5;

{   WRITE COMPLETE BLOCKS FIRST (LEFTOVERS LATER) }
    For j:=1 to limit do
    begin
      Strcopy(s,'           ');
      For k:= 1+(j-1)*5 to 5+(j-1)*5 do
      begin
        Str(k:3,s2); StrCat(s,s2); StrCat(s,'        ');
      end;
      Msg(s); Line0;
      For i:=1 to nrowpr do
      begin
        strcopy(s,'');
        Str(i:3,s2); StrCat(s,s2); StrCat(s,'     ');
        For k:= 1+(j-1)*5 to 5+(j-1)*5 do
        begin
          if abs(a^[i,k]) < 100000 then
          begin
            Str(a^[i,k]:10:3,s1); StrCat(s,s1); StrCat(s,' ')
          end
          else
          begin
            Str(a^[i,k]:10:-3,s1); StrCat(s,s1); StrCat(s,' ')
          end
        end;
        Msg(s)
      end;
      Line0
    end;

{   WRITE REMAINING ELEMENTS }
    Strcopy(s,'           ');
    For k:= 5*limit+1 to ncolpr do
    begin
      Str(k:3,s2); StrCat(s,s2); StrCat(s,'        ');
    end;
    Msg(s); Line0;
    For i:=1 to nrowpr do
    begin
      strcopy(s,'');
      Str(i:3,s1); StrCat(s,s1); StrCat(s,'     ');
      For k:= 5*limit+1 to ncolpr do
      begin
        if abs(a^[i,k]) < 100000.0 then
        begin
          Str(a^[i,k]:10:3,s1); StrCat(s,s1); StrCat(s,' ')
        end
        else
        begin
          Str(a^[i,k]:10:-3,s1); StrCat(s,s1); StrCat(s,' ')
        end
      end;
      Msg(s)
    end;
    Line0
  end;

return: End; {matprt}


Procedure maxst (var overfl:boolean; maxexp, n, nunit, output:integer;
                 epsmch:double; var maxstp:double; mstpf:double;
                 scalex, xc: pVec);
{--------------------------------------------------------------------------
!
!    FEB. 11, 1991
!
!    THIS SUBROUTINE ESTABLISHES A MAXIMUM STEP LENGTH BASED ON THE 2-NORMS
!    OF THE INITIAL ESTIMATES AND THE SCALING FACTORS MULTIPLIED BY A
!    FACTOR MSTPF.
!
!         MAXSTP=MSTPF*MAX(NORM1,NORM2)
!
!              WHERE   MSTPF  USER-CHOSEN FACTOR (DEFAULT: 1000)
!                      NORM1  2-NORM OF SCALED STARTING ESTIMATES
!                      NORM2  2-NORM OF COMPONENT SCALING FACTORS
!-------------------------------------------------------------------------}
Label return;
Var
    norm1, norm2: double;
    wv1:pVec;
    k:integer;
Begin
  New(wv1);
  overfl:=false;
  For k:=1 to n do wv1^[k]:=scalex^[k]*xc^[k];

  twonrm (overfl, maxexp, n, epsmch, norm1, wv1);

  IF (overfl) THEN
  begin
    maxstp:=Power(ten, maxexp);
    IF (output > 2) AND (NOT wrnsup) THEN
    begin
      Line0;
      Writeln(fp_out,'   *    WARNING: NORM OF SCALED INITIAL ESTIMATE SET TO ',norm1:12:3,'     *');
      Line0;
      Writeln(fp_out,'   *       MAXIMUM STEP SIZE, MAXSTP, SET TO ',maxstp:12:3,'                *')
    end;
    goto RETURN
  end;

  twonrm (overfl, maxexp, n, epsmch, norm2, scalex);

  IF (overfl) THEN
  begin
    maxstp:=Power(ten, maxexp);
    IF (output > 2) AND (NOT wrnsup) THEN
    begin
      Line0;
      Writeln(fp_out,'   *    WARNING: NORM OF SCALING FACTORS SET TO ',norm2:12:3,'             *');
      Line0;
      Writeln(fp_out,'   *       MAXIMUM STEP SIZE, MAXSTP, SET TO ',maxstp:12:3,'                *')
    end;
    goto RETURN
  end;

  maxstp:=mstpf*MAX(norm1,norm2);

return: Dispose(wv1)
End; {maxst}

Procedure nersl (newton, restrt, sclfch, sclxch:boolean; acpcod, jupdm,
                 n, nunit, output:integer; fcnnew:double; fvec, xplus:pVec);
{-------------------------------------------------------
!    SEPT. 2, 1991
!
!    THE RESULTS OF EACH ITERATION ARE PRINTED.
!------------------------------------------------------}
Var i:integer;
Begin
  Line0;
  Line0;
  IF (NOT sclxch) THEN
    Msg('    SUMMARY OF ITERATION RESULTS')
  ELSE
    Msg('    SUMMARY OF ITERATION RESULTS (X''S GIVEN IN UNSCALED UNITS)');
  Line0;
  Msg('        UPDATED ESTIMATES               UPDATED FUNCTION VALUES');
  Line0;
  For i:=1 to n do
    if (Abs(xplus^[i])<100000) and (Abs(fvec^[i])<100000) then
      Writeln(fp_out,'   *      X(',i:3,') = ',xplus^[i]:12:3,'               F(',i:3,') = ',fvec^[i]:12:3,'        *')
    else
      Writeln(fp_out,'   *      X(',i:3,') = ',xplus^[i]:12:-3,'               F(',i:3,') = ',fvec^[i]:12:-3,'        *');
  Line0;
  IF (NOT sclfch) THEN
    if Abs(fcnnew)<100000 then
      Writeln(fp_out,'   *      OBJECTIVE FUNCTION VALUE: ',fcnnew:12:3,'                           *')
    else
      Writeln(fp_out,'   *      OBJECTIVE FUNCTION VALUE: ',fcnnew:12:-3,'                           *')
  ELSE
    if Abs(fcnnew)<100000 then
      Writeln(fp_out,'   *      SCALED OBJECTIVE FUNCTION VALUE: ',fcnnew:12:3,'                    *')
    else
      Writeln(fp_out,'   *      SCALED OBJECTIVE FUNCTION VALUE: ',fcnnew:12:-3,'                    *');
  Line0;
  IF (output > 3) AND (NOT newton) THEN
    Writeln(fp_out,'   *      STEP ACCEPTANCE CODE, ACPCOD: ',acpcod:9,'                          *');
  IF (restrt) AND (output >= 3) AND (jupdm <> 0) THEN
  begin
    IF output > 3 Then Line0;
    Msg('      NOTE: JACOBIAN EVALUATED EXPLICITLY AT THIS STEP.');
    Line0
  end
End; {nersl}

Procedure nestop (var absnew, linesr, newton, sclfch, sclxch:boolean; var
                  acptcr, itnum, n, nac1, nac2, nac12, nfunc, njetot:integer;
                  nunit, output, stopcr:integer; var trmcod:integer; var
                  fcnnew:double; ftol, nsttol:double; var stpmax:double;
                  stptol:double; fvec, scalef, scalex, xc, xplus: pVec);
{----------------------------------------------------------------------------
!      FEB. 23, 1992
!
! *** THIS SUBROUTINE CHECKS TO SEE IF THE CONVERGENCE CRITERIA HAVE BEEN MET
!     AND PRINTS MAIN RESULTS TO OUTPUT FILE.  ***
!---------------------------------------------------------------------------}
Label 100, return;
Var
    MAX1, max2, ratio1, SUM: double;
    i: integer;
Begin

  IF output > 3 THEN
  begin
    Line0;
    Msg('    CONVERGENCE TESTING');
    Line0
  end;

{   IF THE NEWTON STEP WAS WITHIN TOLERANCE THEN TRMCOD IS 1 }

  IF trmcod = 1 THEN
  begin

    IF output > 3 THEN
    begin
      Line0;
      IF (NOT sclxch) THEN
        Writeln(fp_out,'   *    MAXIMUM NEWTON STEP LENGTH STPMAX: ',stpmax:12:3,'                  *')
      ELSE
        Writeln(fp_out,'   *    MAXIMUM SCALED NEWTON STEP LENGTH STPMAX: ',stpmax:12:3,'           *');
      Writeln(fp_out,'   *      FIRST CONVERGENCE CRITERION MET;  NSTTOL IS: ',nsttol:12:3,'     *');
      Line0;

{   SKIP CHECKING OTHER STEP SIZE CRITERION AS TRMCOD IS ALREADY 1 }
      GOTO 100

    end
  end;

{   IF THE NEWTON STEP WAS NOT WITHIN TOLERANCE THEN, IF
    STOPCR IS NOT EQUAL TO 2, THE SECOND STEP SIZE STOPPING
    CRITERION MUST BE CHECKED.  }

  IF (stopcr <> 2) AND (trmcod <> 1) THEN
  begin
    MAX1:=zero;
    For i:=1 to n do
    begin
      ratio1:=(ABS(xplus^[i]-xc^[i]))/MAX(ABS(xplus^[i]), one/scalex^[i]);
      MAX1:=MAX(MAX1, ratio1);
      IF output > 4 THEN
        Writeln(fp_out,'   *    RELATIVE STEP SIZE (',i:3,') = ',ratio1:12:3,'                            *')
    end;
    IF output > 4 then Line0;
    IF output > 3 THEN
    begin
      IF (NOT sclxch) THEN
        Writeln(fp_out,'   *    MAXIMUM STEP SIZE: ',MAX1:12:3,'   STPTOL: ',stptol:11:-3,'              *')
      ELSE
        Writeln(fp_out,'   *    MAXIMUM RELATIVE STEP SIZE: ',MAX1:12:3,'   STPTOL: ',stptol:11:-3,'     *');
      Line0
    end;
    IF MAX1 < stptol THEN trmcod:=1
  end;

{   NOTE: CONTINUATION AT LABEL 100 MEANS THAT TRMCOD WAS 1 ON ENTRY
          SO THE STEP SIZE CRITERION ABOVE DID NOT NEED TO BE CHECKED.

    THE SECOND STOPPING CRITERION IS CHECKED IF NEEDED.  }

100:IF (stopcr = 2) OR (stopcr = 12) OR ((stopcr = 3) AND (trmcod = 1)) THEN
  begin
    max2:=zero;

    For i:=1 to n do
    begin
      max2:=MAX(max2,scalef^[i]*ABS(fvec^[i]));
      IF output > 4 THEN
        IF (NOT sclfch) THEN
          Writeln(fp_out,'   *   ABSOLUTE FUNCTION VECTOR (',i:3,') = ',ABS(fvec^[i]):12:-3,'                       *')
        ELSE
          Writeln(fp_out,'   *   ABSOLUTE SCALED FUNCTION VECTOR (',i:3,') = ',(scalef^[i]*ABS(fvec^[i])):12:-3,
                         '                *');
    end;

    IF output > 4 then Line0;
    IF output > 3 THEN
    begin
      IF (NOT sclfch) THEN
        Writeln(fp_out,'   *    MAXIMUM ABSOLUTE FUNCTION: ',max2:12:-3,'     FTOL: ',ftol:11:-3,'      *')
      ELSE
        Writeln(fp_out,'   *    MAX ABSOLUTE SCALED FUNCTION: ',max2:12:-3,'     FTOL: ',ftol:11:-3,'   *');
      Line0
    end;

    IF max2 < ftol THEN
      IF (stopcr = 3) AND (trmcod = 1) THEN

{     BOTH NEEDED STOPPING CRITERIA HAVE BEEN MET }
        trmcod:=3

      ELSE IF (stopcr = 12) AND (trmcod = 1) THEN

{     BOTH STOPPING CRITERIA HAVE BEEN MET ALTHOUGH
      EITHER ONE WOULD BE SATISFACTORY.  }
        trmcod:=12

      ELSE IF (stopcr = 2) OR (stopcr = 12) THEN
        trmcod:=2
    ELSE IF stopcr = 3 THEN

{     ONLY THE FIRST STOPPING CRITERION WAS MET - TRMCOD
      MUST BE RESET FROM 1 BACK TO 0. }
      trmcod:=0

  end;

{    PRINT FINAL RESULTS IF CONVERGENCE REACHED }

  IF trmcod > 0 THEN
  begin
    IF output > 0 THEN
    begin
      IF output = 1 then Line1;
      Line0;
      Writeln(fp_out,'   *    CONVERGENCE REACHED; TERMINATION CODE: ...............',trmcod:6,'       *');

{   ITERATION RESULTS NOT PRINTED IN NERSL }

      Line0; Line0;
      IF (sclfch) OR (sclxch) THEN
      begin
        Msg('                           UNSCALED RESULTS');
        Line0
      end;
      Msg('          FINAL ESTIMATES                 FINAL FUNCTION VALUES');
      Line0;
      For i:=1 to n do
        Writeln(fp_out,'   *    X(',i:3,') = ',xplus^[i]:15:6,'            F(',i:3,') = ',fvec^[i]:15:6,'       *');
      Line0;

      IF (sclfch) THEN
      begin
{       NEED UNSCALED OBJECTIVE FUNCTION }
        sum := DOT_PRODUCT(n, fvec, fvec);
        Writeln(fp_out,'   *    FINAL OBJECTIVE FUNCTION VALUE: ',(sum/two):12:-3,'                       *')
      end
      ELSE
        Writeln(fp_out,'   *    FINAL OBJECTIVE FUNCTION VALUE: ',fcnnew:12:-3,'                       *');

      IF (sclfch) OR (sclxch) THEN
      begin
        Line0; Line0;
        Msg('                           SCALED RESULTS');
        Line0;
        Msg('          FINAL ESTIMATES                  FINAL FUNCTION VALUES');
        Line0;
        For i:=1 to n do
          Writeln(fp_out,'   *    X(',i:3,') = ',scalex^[i]*xplus^[i]:15:6,'            F(',i:3,') = ',
          scalef^[i]*fvec^[i]:15:6,'      *');
        Line0;
        Writeln(fp_out,'   *    FINAL OBJECTIVE FUNCTION VALUE: ',fcnnew:12:-3,'                       *')
      end
    end
  end
  ELSE
  begin
    IF output > 3 THEN
    begin
      Msg('    CONVERGENCE NOT REACHED.'); Line0
    end;
    goto RETURN
  end;

{ TERMINATION HAS BEEN REACHED }

  IF output > 0 THEN
  begin
    Line0; Line0;
    Writeln(fp_out,'   *    TOTAL NUMBER OF ITERATIONS: ..........................',itnum:6,'       *');
    IF (NOT newton) AND (NOT absnew) THEN
      IF (linesr) THEN
        Writeln(fp_out,'   *    TOTAL NUMBER OF LINE SEARCH FUNCTION EVALUATIONS: ....',nfunc:6,'       *')
      ELSE
        Writeln(fp_out,'   *    TOTAL NUMBER OF TRUST REGION FUNCTION EVALUATIONS: ...',nfunc:6,'       *');

    Writeln(fp_out,'   *    TOTAL NUMBER OF EXPLICIT JACOBIAN EVALUATIONS: .......',njetot:6,'       *');
    Writeln(fp_out,'   *    TOTAL NUMBER OF FUNCTION EVALUATIONS: ................',nfetot:6,'       *');
    Line0;

    IF (NOT newton) AND (NOT absnew) AND (acptcr <> 1) AND (output > 2) THEN
    begin
      Writeln(fp_out,'   *    NUMBER OF STEPS ACCEPTED BY FUNCTION VALUE ONLY: .....',nac1:6,'       *');
      Writeln(fp_out,'   *    NUMBER OF STEPS ACCEPTED BY STEP SIZE VALUE ONLY: ....',nac2:6,'       *');
      Writeln(fp_out,'   *    NUMBER OF STEPS ACCEPTED BY EITHER CRITERION: ........',nac12:6,'       *');
      Line0
    end;
    Line1
  end;

return:End; {nestop}


END.

{end of file unnes1.pas}