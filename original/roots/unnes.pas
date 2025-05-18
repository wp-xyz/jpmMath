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
*                                     Pascal Release By J-P Moreau, Paris.   *
*                                              (PART 1/3)                    *
*                                           (www.jpmoreau.fr)                *
*****************************************************************************}
Unit UNNES;

Interface

Uses WinCrt, Utils, Fcn;

Var
    bypass, matsup, wrnsup: Boolean;
    itnum, nfunc:integer;
    smallb, bigb, smalls, bigs, bigr: Double;


    Procedure ascalf (n:Integer; epsmch:Double; fvecc:pVec; jac:pMat; scalef:pVec);
    Procedure ascalx (n:integer; epsmch:double; jac:pMat; scalex:pVec);
    Procedure ataov (Var overfl:boolean; maxexp, n, nunit, output:integer;
                     a, b:pMat; scalef:pVec);
    Procedure atvov (Var overfl:boolean; maxexp, n, nunit, output:integer;
                     amat: pMat; bvec, cvec: pVec);
    Procedure atbmul (ncadec, ncaact, ncbdec, ncbact, ncdec, ncact:integer;
                      amat, bmat, cmat: pMat);
    Procedure avmul (nradec, nraact, ncdec, ncact:integer; amat:pMat; bvec, cvec:pVec);
    Procedure bakdif (Var overfl:boolean; j, n:integer; Var deltaj:double; tempj:double;
                      fvec, fvecj1:pVec; jacfdm:pMat; xc:pVec);
    Procedure bnddif (Var overfl:boolean; j, n:integer; epsmch:double; boundl, boundu,
                      fvecc, fvecj1:pVec; jacfdm:pMat; xc:pVec);
    Procedure broyfa (Var overch, overfl:boolean; sclfch, sclxch:boolean;  maxexp, n,
                      nunit, output:integer; epsmch:double; a:pMat;  delf, fvec,
                      fvecc:pVec; jac:pMat; rdiag, s, scalef, scalex, t, w, xc,
                      xplus:pVec);
    Procedure broyun (var overfl:boolean; maxexp, n, nunit, output:integer;
                      var epsmch:double; fvec, fvecc:pVec; jac:pMat;
                      scalex, xc, xplus:pVec);
    Procedure cholde (n:integer; var maxadd, maxffl:double; sqrtep:double; h, l: pMat);
    Procedure chsolv (overch:boolean; var overfl:boolean; maxexp, n, nunit,
                      output:integer; l:pMat; rhs, s:pVec);
    Procedure condno (overch:boolean; var overfl:boolean; maxexp, n, nunit,
                      output:integer; var connum:double; a:pMat; rdiag:pVec);
    Procedure delcau (cauchy, overch:boolean; var overfl:boolean; itnum,
                      maxexp, n, nunit, output:integer; var beta, caulen,
                      delta, epsmch, maxstp, newlen, sqrtz:double;
                      a: pMat; delf, scalex: pVec);
    Procedure deufls (var abort:boolean; deuflh, geoms:boolean; var overch,overfl,qnfail:boolean;
                      qrsing:boolean; var restrt:boolean; sclfch:boolean; var sclxch:boolean;
                      var acpcod:integer; acptcr, contyp, itnum,jupdm, maxexp, maxlin, n:integer;
                      var nfunc:integer; nunit, output, qnupdm:integer; var stopcr:integer;
                      alpha, confac, delfts, epsmch, fcnmax:double; var fcnnew:double;
                      fcnold:double; var lambda:double; newmax:double; var sbrnrm:double; sigma:double;
                      a, astore:pMat; boundl, boundu, delf, fvec, hhpi:pVec; jac:pMat; rdiag, rhs, s,
                      sbar, scalef, scalex, sn, xc, xplus: pVec);
    Procedure fcnevl (var overfl:boolean; maxexp, n, nunit, output:integer;
                      epsmch:double; var fcnnew:double; fvec, scalef:pVec);
    Procedure fordif (var overfl:boolean; j, n:integer; deltaj:double; fvec, fvecj1:pVec;
                      jacfdm:pMat; xc:pVec);
    Procedure innerp (overch:boolean; var overfl:boolean; maxexp, ldima, ldimb, n, nunit,
                      output:integer; var dtpro:double; a, b:pVec);
    Procedure matcop (nradec, nraact, ncadec, ncaact, nrbdec, ncbdec:integer; amat, bmat:pMat);
    Procedure qrsolv(overch:boolean; var overfl:boolean; maxexp, n, nunit, output:integer;
                     a:pMat; hhpi, rdiag, b: pVec);
    Procedure qrupda (var overfl:boolean; maxexp, n:integer; epsmch:double;
                      a, jac:pMat; u, v:pVec);
    Procedure rsolv (overch:boolean; var overfl:boolean; maxexp, n, nunit,
                     output:integer; a:pMat; rdiag, b:pVec);
    Procedure setup (var absnew, cauchy, deuflh, geoms, linesr, newton, overch:boolean;
                     var acptcr, itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns,
                     maxqns, minqns:integer; n:integer; var narmij, niejev, njacch,
                     output, qnupdm, stopcr, supprs, trupdm: integer; var alpha, confac,
                     delta, delfac, epsmch, etafac, fdtolj, ftol, lam0, mstpf, nsttol,
                     omega, ratiof, sigma, stptol: double; boundl, boundu, scalef,
                     scalex: pVec; var help:String);
    Procedure twonrm (var overfl:boolean; maxexp, n:integer; epsmch:double;
                      var eucnrm:double; v:pVec);
    Procedure update (mnew:integer; var mold:integer; n:integer; var trmcod:
                      integer; var fcnnew, fcnold:double;  fvec, fvecc, xc, xplus:pVec);
    Procedure utbmul (ncadec, ncaact, ncbdec, ncbact, ncdec, ncact:integer;
                      amat, bmat, cmat:pMat);
    Procedure uvmul (nradec, nraact, ncdec, ncact:integer; amat:pMat; bvec, cvec:pVec);


Implementation


Procedure abmul (nradec, nraact, ncbdec, ncbact, ncdec, ncact: Integer;
                 amat, bmat, cmat: pMat; arow:pVec);
{------------------------------------------------------------------------
!    FEB. 8, 1991
!
!    MATRIX MULTIPLICATION AB=C
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4
!    EACH ROW OF MATRIX A IS SAVED AS A COLUMN, AROW, BEFORE USE.
!
!    NRADEC IS 1ST DIM. OF AMAT; NRAACT IS ACTUAL LIMIT FOR 1ST INDEX

!    NCBDEC IS 2ND DIM. OF BMAT; NCBACT IS ACTUAL LIMIT FOR 2ND INDEX
!    NCDEC IS COMMON DIMENSION OF AMAT & BMAT; NCACT IS ACTUAL LIMIT
!
!    I.E. NRADEC IS THE NUMBER OF ROWS OF A DECLARED
!         NCBDEC IS THE NUMBER OF COLUMNS OF B DECLARED
!         NCDEC IS THE COMMON DECLARED DIMENSION
!
!    MODIFICATION OF MATRIX MULTIPLIER DONATED BY PROF. JAMES MACKINNON,
!    QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!-----------------------------------------------------------------------}
Var

    i, j, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, ncc32r: Integer;
    SUM: Double;

Begin
{       FIND NUMBER OF GROUPS OF SIZE 32, 16... }

  ncc32:=ncact Div 32;
  ncc32r:=ncact - 32*ncc32;
  ncc16:=ncc32r Div 16;
  ncc16r:=ncc32r - 16*ncc16;
  ncc8:=ncc16r Div 8;
  ncc8r:=ncc16r - 8*ncc8;
  ncc4:=ncc8r Div 4;
  ncc4r:=ncc8r - 4*ncc4;

{ REASSIGN ROWS TO VECTOR AROW }

  For  i:=1 to nraact do
  begin
    k:=0;
    IF ncc32 > 0 THEN
      For  kk:=1 to ncc32 do
      begin
      k:=k+32;
      arow^[k-31]:=amat^[i,k-31];
      arow^[k-30]:=amat^[i,k-30];
      arow^[k-29]:=amat^[i,k-29];
      arow^[k-28]:=amat^[i,k-28];
      arow^[k-27]:=amat^[i,k-27];
      arow^[k-26]:=amat^[i,k-26];
      arow^[k-25]:=amat^[i,k-25];
      arow^[k-24]:=amat^[i,k-24];
      arow^[k-23]:=amat^[i,k-23];
      arow^[k-22]:=amat^[i,k-22];
      arow^[k-21]:=amat^[i,k-21];
      arow^[k-20]:=amat^[i,k-20];
      arow^[k-19]:=amat^[i,k-19];
      arow^[k-18]:=amat^[i,k-18];
      arow^[k-17]:=amat^[i,k-17];
      arow^[k-16]:=amat^[i,k-16];
      arow^[k-15]:=amat^[i,k-15];
      arow^[k-14]:=amat^[i,k-14];
      arow^[k-13]:=amat^[i,k-13];
      arow^[k-12]:=amat^[i,k-12];
      arow^[k-11]:=amat^[i,k-11];
      arow^[k-10]:=amat^[i,k-10];
      arow^[k-9]:=amat^[i,k-9];
      arow^[k-8]:=amat^[i,k-8];
      arow^[k-7]:=amat^[i,k-7];
      arow^[k-6]:=amat^[i,k-6];
      arow^[k-5]:=amat^[i,k-5];
      arow^[k-4]:=amat^[i,k-4];
      arow^[k-3]:=amat^[i,k-3];
      arow^[k-2]:=amat^[i,k-2];
      arow^[k-1]:=amat^[i,k-1];
      arow^[k]:=amat^[i,k]
      end;
    IF ncc16 > 0 THEN
      For  kk:=1 to ncc16 do
      begin
      k:=k+16;
      arow^[k-15]:=amat^[i,k-15];
      arow^[k-14]:=amat^[i,k-14];
      arow^[k-13]:=amat^[i,k-13];
      arow^[k-12]:=amat^[i,k-12];
      arow^[k-11]:=amat^[i,k-11];
      arow^[k-10]:=amat^[i,k-10];
      arow^[k-9]:=amat^[i,k-9];
      arow^[k-8]:=amat^[i,k-8];
      arow^[k-7]:=amat^[i,k-7];
      arow^[k-6]:=amat^[i,k-6];
      arow^[k-5]:=amat^[i,k-5];
      arow^[k-4]:=amat^[i,k-4];
      arow^[k-3]:=amat^[i,k-3];
      arow^[k-2]:=amat^[i,k-2];
      arow^[k-1]:=amat^[i,k-1];
      arow^[k]:=amat^[i,k]
      end;
    IF ncc8 > 0 THEN
      For  kk:=1 to ncc8 do
      begin
      k:=k+8;
      arow^[k-7]:=amat^[i,k-7];
      arow^[k-6]:=amat^[i,k-6];
      arow^[k-5]:=amat^[i,k-5];
      arow^[k-4]:=amat^[i,k-4];
      arow^[k-3]:=amat^[i,k-3];
      arow^[k-2]:=amat^[i,k-2];
      arow^[k-1]:=amat^[i,k-1];
      arow^[k]:=amat^[i,k]
      end;
    IF ncc4 > 0 THEN
    For  kk:=1 to ncc4 do
    begin
      k:=k+4;
      arow^[k-3]:=amat^[i,k-3];
      arow^[k-2]:=amat^[i,k-2];
      arow^[k-1]:=amat^[i,k-1];
      arow^[k]:=amat^[i,k]
    end;
    IF ncc4r > 0 THEN
    For  kk:=1 to ncc4r do
    begin
      k:=k+1;
      arow^[k]:=amat^[i,k]
    end;

{   FIND ENTRY FOR MATRIX C USING COLUMN VECTOR AROW }

    For  j:=1 to ncbact do
    begin
      sum:=zero;
      k:=0;
      IF ncc32 > 0 THEN
        For kk:=1 to ncc32 do
        begin
          k:=k+32;
          sum:=sum + arow^[k-31]*bmat^[k-31,j] + arow^[k-30]*bmat^[k-30,j] +
            arow^[k-29]*bmat^[k-29,j] + arow^[k-28]*bmat^[k-28,j] +
            arow^[k-27]*bmat^[k-27,j] + arow^[k-26]*bmat^[k-26,j] +
            arow^[k-25]*bmat^[k-25,j] + arow^[k-24]*bmat^[k-24,j];
          sum:=sum + arow^[k-23]*bmat^[k-23,j] + arow^[k-22]*bmat^[k-22,j] +
            arow^[k-21]*bmat^[k-21,j] + arow^[k-20]*bmat^[k-20,j] +
            arow^[k-19]*bmat^[k-19,j] + arow^[k-18]*bmat^[k-18,j] +
            arow^[k-17]*bmat^[k-17,j] + arow^[k-16]*bmat^[k-16,j];
          sum:=sum + arow^[k-15]*bmat^[k-15,j] + arow^[k-14]*bmat^[k-14,j] +
            arow^[k-13]*bmat^[k-13,j] + arow^[k-12]*bmat^[k-12,j] +
            arow^[k-11]*bmat^[k-11,j] + arow^[k-10]*bmat^[k-10,j] +
            arow^[k-9]*bmat^[k-9,j] + arow^[k-8]*bmat^[k-8,j];
          sum:=sum  +  arow^[k-7]*bmat^[k-7,j] + arow^[k-6]*bmat^[k-6,j] +
            arow^[k-5]*bmat^[k-5,j] + arow^[k-4]*bmat^[k-4,j] +
            arow^[k-3]*bmat^[k-3,j] + arow^[k-2]*bmat^[k-2,j] +
            arow^[k-1]*bmat^[k-1,j] + arow^[k]*bmat^[k,j]
        end;

      IF ncc16 > 0 THEN
        For  kk:=1 to ncc16 do
        begin
          k:=k+16;
          sum:=sum + arow^[k-15]*bmat^[k-15,j] + arow^[k-14]*bmat^[k-14,j] +
            arow^[k-13]*bmat^[k-13,j] + arow^[k-12]*bmat^[k-12,j] +
            arow^[k-11]*bmat^[k-11,j] + arow^[k-10]*bmat^[k-10,j] +
            arow^[k-9]*bmat^[k-9,j] + arow^[k-8]*bmat^[k-8,j];
          sum:=sum + arow^[k-7]*bmat^[k-7,j] + arow^[k-6]*bmat^[k-6,j] +
            arow^[k-5]*bmat^[k-5,j] + arow^[k-4]*bmat^[k-4,j] +
            arow^[k-3]*bmat^[k-3,j] + arow^[k-2]*bmat^[k-2,j] +
            arow^[k-1]*bmat^[k-1,j] + arow^[k]*bmat^[k,j]
        end;

      IF (ncc8 > 0) THEN
        For kk:=1 to ncc8 do
        begin
          k:=k+8;
          sum:=sum + arow^[k-7]*bmat^[k-7,j] + arow^[k-6]*bmat^[k-6,j] +
            arow^[k-5]*bmat^[k-5,j] + arow^[k-4]*bmat^[k-4,j] +
            arow^[k-3]*bmat^[k-3,j] + arow^[k-2]*bmat^[k-2,j] +
            arow^[k-1]*bmat^[k-1,j] + arow^[k]*bmat^[k,j]
        end;

      IF (ncc4 > 0) THEN
        For  kk:=1 to ncc4 do
        begin
          k:=k+4;
          sum:=sum + arow^[k-3]*bmat^[k-3,j] + arow^[k-2]*bmat^[k-2,j] +
            arow^[k-1]*bmat^[k-1,j] + arow^[k]*bmat^[k,j]
        end;

      IF ncc4r > 0 THEN
        For  kk:=1 to ncc4r do
        begin
          k:=k+1;
          sum:=sum + arow^[k]*bmat^[k,j]
        end;

      cmat^[i,j]:=sum

    end  {j loop} 
  end  {i loop}

End; {abmul}



Procedure ascalf (n:Integer; epsmch:Double; fvecc:pVec; jac:pMat; scalef:pVec);
{------------------------------------------------------------------------------
!    FEB. 13, 1991
!
!    THIS SUBROUTINE ESTABLISHES SCALING FACTORS FOR THE RESIDUAL VECTOR
!    IF FUNCTION ADAPTIVE SCALING IS CHOSEN USING INTEGER VARIABLE ITSCLF.
!
!    NOTE: IN QUASI-NEWTON METHODS THE SCALING FACTORS ARE
!          UPDATED ONLY WHEN THE JACOBIAN IS EVALUATED EXPLICITLY.
!
!    SCALING FACTORS ARE DETERMINED FROM THE INFINITY NORMS OF THE ROWS
!    OF THE JACOBIAN AND THE VALUES OF THE CURRENT FUNCTION VECTOR, FVECC.
!
!    A MINIMUM TOLERANCE ON THE SCALING FACTOR IS THE SQUARE
!    ROOT OF THE MACHINE PRECISION, SQRTEP.
!------------------------------------------------------------------------------}
Var
    amax, sqrtep: double;
    i, j: integer;
Begin

  sqrtep:=SQRT(epsmch);

{ I COUNTS THE ROWS }

  For  i:=1 to n do
  begin
    amax:=zero;

{ FIND MAXIMUM ENTRY IN ROW I }

    For j:=1 to n do amax:=MAX(amax, ABS(jac^[i,j]));

    amax:=MAX(amax, fvecc^[i]);

{ SET SCALING FACTOR TO A DEFAULT OF ONE IF ITH ROW IS ZEROS }

    IF amax = zero then amax:=one;
    amax:=MAX(amax, sqrtep);
    scalef^[i]:=one/amax

  end

End; {ascalf}


Procedure ascalx (n:integer; epsmch:double; jac:pMat; scalex:pVec);
{-------------------------------------------------------------------------
!    FEB. 13, 1991
!
!    THIS SUBROUTINE ESTABLISHES SCALING FACTORS FOR THE COMPONENET VECTOR
!    IF ADAPTIVE SCALING IS CHOSEN USING INTEGER ITSCLX.
!
!    NOTE: IN QUASI-NEWTON METHODS THE SCALING FACTORS ARE
!          UPDATED ONLY WHEN THE JACOBIAN IS EVALUATED EXPLICITLY.
!
!    SCALING FACTORS ARE DETERMINED FROM THE INFINITY NORMS
!    OF THE COLUMNS OF THE JACOBIAN.
!
!    A MINIMUM TOLERANCE ON THE SCALING FACTOR IS THE SQUARE
!    ROOT OF THE MACHINE PRECISION, SQRTEP.
!------------------------------------------------------------------------}
Var
    amax, sqrtep: double;
    i, j: integer;

Begin
  sqrtep:=SQRT(epsmch);

{ J COUNTS COLUMNS }

  For j:=1 to n do
  begin

    amax:=zero;

{   FIND MAXIMUM ENTRY IN JTH COLUMN }

    For  i:=1 to n do  amax:=MAX(amax, ABS(jac^[i,j]));

{   IF A COLUMN IS ALL ZEROS SET AMAX TO ONE }

    IF amax = zero then amax:=one;
    scalex^[j]:=MAX(amax, sqrtep)
  end

End; {ascalx}



Procedure atamul (nradec, ncadec, nraact, ncaact, nrbdec, ncbdec: integer;
                  amat, bmat: pMat);
{-----------------------------------------------------------------------
!    FEB. 8, 1991
!
!    MATRIX MULTIPLICATION:   A^A:=B
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4.
!
!    NRADEC IS NUMBER OF ROWS IN A DECLARED
!    NCADEC IS NUMBER OF COLUMNS IN A DECLARED
!    NRAACT IS THE LIMIT FOR THE 1ST INDEX IN A
!    NCAACT IS THE LIMIT FOR THE 2ND INDEX IN A
!    NRBDEC IS NUMBER OF ROWS IN B DECLARED
!    NCBDEC IS NUMBER OF COLUMNS IN B DECLARED
!
!    MODIFIED VERSION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES
!    MACKINNON, QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!---------------------------------------------------------------------}
Var
     SUM: Double;
     i,j,k,kk,ncc4,ncc4r,ncc8,ncc8r,ncc16, ncc16r, ncc32, ncc32r: integer;
Begin

{ FIND NUMBER OF GROUPS OF SIZE 32, 16...}
  ncc32:=nraact Div 32;
  ncc32r:=nraact - 32*ncc32;
  ncc16:=ncc32r Div 16;
  ncc16r:=ncc32r - 16*ncc16;
  ncc8:=ncc16r Div 8;
  ncc8r:=ncc16r - 8*ncc8;
  ncc4:=ncc8r Div 4;
  ncc4r:=ncc8r - 4*ncc4;

{ FIND ENTRY IN MATRIX B }

  For  i:=1 to ncaact do
  begin
    For  j:=i to ncaact do
    begin
      sum:=zero;
      k:=0;
      IF ncc32 > 0 THEN
        For kk:=1 to ncc32 do
        begin
          k:=k+32;
          sum:=sum + amat^[k-31,i]*amat^[k-31,j] + amat^[k-30,i]*amat^[k-30,j] +
            amat^[k-29,i]*amat^[k-29,j] + amat^[k-28,i]*amat^[k-28,j] +
            amat^[k-27,i]*amat^[k-27,j] + amat^[k-26,i]*amat^[k-26,j] +
            amat^[k-25,i]*amat^[k-25,j] + amat^[k-24,i]*amat^[k-24,j];
          sum:=sum + amat^[k-23,i]*amat^[k-23,j] + amat^[k-22,i]*amat^[k-22,j] +
            amat^[k-21,i]*amat^[k-21,j] + amat^[k-20,i]*amat^[k-20,j] +
            amat^[k-19,i]*amat^[k-19,j] + amat^[k-18,i]*amat^[k-18,j] +
            amat^[k-17,i]*amat^[k-17,j] + amat^[k-16,i]*amat^[k-16,j];
          sum:=sum + amat^[k-15,i]*amat^[k-15,j] + amat^[k-14,i]*amat^[k-14,j] +
            amat^[k-13,i]*amat^[k-13,j] + amat^[k-12,i]*amat^[k-12,j] +
            amat^[k-11,i]*amat^[k-11,j] + amat^[k-10,i]*amat^[k-10,j] +
            amat^[k-9,i]*amat^[k-9,j] + amat^[k-8,i]*amat^[k-8,j];
          sum:=sum + amat^[k-7,i]*amat^[k-7,j] + amat^[k-6,i]*amat^[k-6,j] +
            amat^[k-5,i]*amat^[k-5,j] + amat^[k-4,i]*amat^[k-4,j] +
            amat^[k-3,i]*amat^[k-3,j] + amat^[k-2,i]*amat^[k-2,j] +
            amat^[k-1,i]*amat^[k-1,j] + amat^[k,i]*amat^[k,j]
        end;

      IF ncc16 > 0 THEN
        For kk:=1 to ncc16 do
        begin
          k:=k+16;
          sum:=sum + amat^[k-15,i]*amat^[k-15,j] + amat^[k-14,i]*amat^[k-14,j] +
            amat^[k-13,i]*amat^[k-13,j] + amat^[k-12,i]*amat^[k-12,j] +
            amat^[k-11,i]*amat^[k-11,j] + amat^[k-10,i]*amat^[k-10,j] +
            amat^[k-9,i]*amat^[k-9,j] + amat^[k-8,i]*amat^[k-8,j];
          sum:=sum + amat^[k-7,i]*amat^[k-7,j] + amat^[k-6,i]*amat^[k-6,j] +
            amat^[k-5,i]*amat^[k-5,j] + amat^[k-4,i]*amat^[k-4,j] +
            amat^[k-3,i]*amat^[k-3,j] + amat^[k-2,i]*amat^[k-2,j] +
            amat^[k-1,i]*amat^[k-1,j] + amat^[k,i]*amat^[k,j]
        end;

      IF ncc8 > 0 THEN
        For kk:=1 to ncc8 do
        begin
          k:=k+8;
          sum:=sum + amat^[k-7,i]*amat^[k-7,j] + amat^[k-6,i]*amat^[k-6,j] +
            amat^[k-5,i]*amat^[k-5,j] + amat^[k-4,i]*amat^[k-4,j] +
            amat^[k-3,i]*amat^[k-3,j] + amat^[k-2,i]*amat^[k-2,j] +
            amat^[k-1,i]*amat^[k-1,j] + amat^[k,i]*amat^[k,j]
        end;

      IF ncc4 > 0 THEN
        For kk:=1 to ncc4 do
        begin
          k:=k+4;
          sum:=sum + amat^[k-3,i]*amat^[k-3,j] + amat^[k-2,i]*amat^[k-2,j] +
            amat^[k-1,i]*amat^[k-1,j] + amat^[k,i]*amat^[k,j]
        end;

      IF ncc4r > 0 THEN
        For kk:=1 to ncc4r do
        begin
          k:=k+1;
          sum:=sum + amat^[k,i]*amat^[k,j]
        end;

      bmat^[i,j]:=sum;
      IF i <> j then bmat^[j,i]:=bmat^[i,j]
    end {j loop}
  End; {i loop}

End; {atamul}


Procedure ataov (Var overfl:boolean; maxexp, n, nunit, output:integer;
                 a, b:pMat; scalef:pVec);
{----------------------------------------------------------------------
!    SEPT. 8, 1991
!
!    THIS SUBROUTINE FINDS THE PRODUCT OF THE TRANSPOSE OF THE MATRIX A
!    AND MATRIX A.  EACH ENTRY IS CHECKED BEFORE BEING ACCEPTED.
!    IF IT WOULD CAUSE AN OVERFLOW 10^MAXEXP IS INSERTED IN ITS PLACE.
!---------------------------------------------------------------------}
Label 40, 70;
Var
    eps, SUM: double;
    i, j, k: integer;
Begin
  eps:=one/Power(ten, maxexp);
  overfl:=false;

  For i:=1 to n do
  begin
    For j:=i+1 to n do
    begin
      sum:=zero;
      For k:=1 to n do
      begin
        IF Ln(ABS(a^[k,i])+eps) + LOG(ABS(a^[k,j])+eps) +
                             two*LOG(scalef^[k]) > maxexp THEN
        begin
          overfl:=true;
          b^[i,j]:=SIGN(Power(ten,maxexp),a^[k,i])*SIGN(one,a^[k,j]);
          IF (output > 2) AND (NOT wrnsup) THEN
          begin
            Line0;
            Writeln(fp_out,'   *    WARNING: COMPONENT MATRIX-MATRIX PRODUCT SET TO ', b^[i,j]:12:3,'     *')
          end;
          GOTO 40;
        end;
        sum:=sum + a^[k,i]*a^[k,j]*scalef^[k]*scalef^[k]
      end; {k loop}
      b^[i,j]:=sum;
      b^[j,i]:=sum;
40: end; {j loop}

    sum:=zero;
    For k:=1 to n do
    begin
      IF two*(LOG(ABS(a^[k,i])+eps) + LOG(scalef^[k])) > maxexp THEN
      begin
        overfl:=true;
        b^[i,i]:=Power(ten,maxexp);
        IF (output > 2) AND (NOT wrnsup) THEN
        begin
          Line0;
          Writeln(fp_out,'   *    WARNING: COMPONENT MATRIX-MATRIX PRODUCT SET TO ', b^[i,i]:12:3,'     *')
        end;
        GOTO 70
      end;
      sum:=sum + a^[k,i]*a^[k,i]*scalef^[k]*scalef^[k]
    end;
    b^[i,i]:=sum;
70: End {i loop}

End; {ataov}


Procedure atbmul (ncadec, ncaact, ncbdec, ncbact, ncdec, ncact:integer;
                  amat, bmat, cmat: pMat);
{----------------------------------------------------------------------
!    FEB. 8, 1991
!
!    MATRIX MULTIPLICATION:   A^B:=C
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4.
!
!    NCADEC IS 2ND DIM. OF AMAT; NCAACT IS ACTUAL LIMIT FOR 2ND INDEX
!    NCBDEC IS 2ND DIM. OF BMAT; NCBACT IS ACTUAL LIMIT FOR 2ND INDEX
!    NCDEC IS COMMON DIMENSION OF AMAT & BMAT; NCACT IS ACTUAL LIMIT
!
!    I.E.   NCADEC IS NUMBER OF COLUMNS OF A DECLARED
!           NCBDEC IS NUMBER OF COLUMNS OF B DECLARED
!           NCDEC  IS THE NUMBER OF ROWS IN BOTH A AND B DECLARED
!
!    MODIFIED VERSION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES
!    MACKINNON, QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!----------------------------------------------------------------------}
Var
    i, j, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32,
    ncc32r: Integer;
    SUM: Double;
Begin
{   FIND NUMBER OF GROUPS OF SIZE 32, 16... }

  ncc32:=ncact Div 32;
  ncc32r:=ncact - 32*ncc32;
  ncc16:=ncc32r Div 16;
  ncc16r:=ncc32r - 16*ncc16;
  ncc8:=ncc16r Div 8;
  ncc8r:=ncc16r - 8*ncc8;
  ncc4:=ncc8r Div 4;
  ncc4r:=ncc8r - 4*ncc4;

{   FIND ENTRY IN MATRIX C }

  For i:=1 to ncaact do
  begin
    For j:=1 to ncbact do
    begin
      sum:=zero;
      k:=0;
      IF ncc32 > 0 THEN
        For kk:=1 to ncc32 do
        begin
          k:=k+32;
          sum:=sum + amat^[k-31,i]*bmat^[k-31,j] + amat^[k-30,i]*bmat^[k-30,j] +
            amat^[k-29,i]*bmat^[k-29,j] + amat^[k-28,i]*bmat^[k-28,j] +
            amat^[k-27,i]*bmat^[k-27,j] + amat^[k-26,i]*bmat^[k-26,j] +
            amat^[k-25,i]*bmat^[k-25,j] + amat^[k-24,i]*bmat^[k-24,j];
          sum:=sum + amat^[k-23,i]*bmat^[k-23,j] + amat^[k-22,i]*bmat^[k-22,j] +
            amat^[k-21,i]*bmat^[k-21,j] + amat^[k-20,i]*bmat^[k-20,j] +
            amat^[k-19,i]*bmat^[k-19,j] + amat^[k-18,i]*bmat^[k-18,j] +
            amat^[k-17,i]*bmat^[k-17,j] + amat^[k-16,i]*bmat^[k-16,j];
          sum:=sum + amat^[k-15,i]*bmat^[k-15,j] + amat^[k-14,i]*bmat^[k-14,j] +
            amat^[k-13,i]*bmat^[k-13,j] + amat^[k-12,i]*bmat^[k-12,j] +
            amat^[k-11,i]*bmat^[k-11,j] + amat^[k-10,i]*bmat^[k-10,j] +
            amat^[k-9,i]*bmat^[k-9,j] + amat^[k-8,i]*bmat^[k-8,j];
          sum:=sum + amat^[k-7,i]*bmat^[k-7,j] + amat^[k-6,i]*bmat^[k-6,j] +
            amat^[k-5,i]*bmat^[k-5,j] + amat^[k-4,i]*bmat^[k-4,j] +
            amat^[k-3,i]*bmat^[k-3,j] + amat^[k-2,i]*bmat^[k-2,j] +
            amat^[k-1,i]*bmat^[k-1,j] + amat^[k,i]*bmat^[k,j]
        end;

      IF (ncc16 > 0) THEN
        For kk:=1 to ncc16 do
        begin
          k:=k+16;
          sum:=sum + amat^[k-15,i]*bmat^[k-15,j] + amat^[k-14,i]*bmat^[k-14,j] +
            amat^[k-13,i]*bmat^[k-13,j] + amat^[k-12,i]*bmat^[k-12,j] +
            amat^[k-11,i]*bmat^[k-11,j] + amat^[k-10,i]*bmat^[k-10,j] +
            amat^[k-9,i]*bmat^[k-9,j] + amat^[k-8,i]*bmat^[k-8,j];
          sum:=sum + amat^[k-7,i]*bmat^[k-7,j] + amat^[k-6,i]*bmat^[k-6,j] +
            amat^[k-5,i]*bmat^[k-5,j] + amat^[k-4,i]*bmat^[k-4,j] +
            amat^[k-3,i]*bmat^[k-3,j] + amat^[k-2,i]*bmat^[k-2,j] +
            amat^[k-1,i]*bmat^[k-1,j] + amat^[k,i]*bmat^[k,j]
        end;

      IF ncc8 > 0 THEN
        For kk:=1 to ncc8 do
        begin
          k:=k+8;
          sum:=sum + amat^[k-7,i]*bmat^[k-7,j] + amat^[k-6,i]*bmat^[k-6,j] +
            amat^[k-5,i]*bmat^[k-5,j] + amat^[k-4,i]*bmat^[k-4,j] +
            amat^[k-3,i]*bmat^[k-3,j] + amat^[k-2,i]*bmat^[k-2,j] +
            amat^[k-1,i]*bmat^[k-1,j] + amat^[k,i]*bmat^[k,j]
        end;

      IF ncc4 > 0 THEN
        For kk:=1 to ncc4 do
        begin
          k:=k+4;
          sum:=sum + amat^[k-3,i]*bmat^[k-3,j] + amat^[k-2,i]*bmat^[k-2,j] +
            amat^[k-1,i]*bmat^[k-1,j] + amat^[k,i]*bmat^[k,j]
        end;

      IF ncc4r > 0 THEN
        For kk:=1 to ncc4r do
        begin
          k:=k+1;
          sum:=sum + amat^[k,i]*bmat^[k,j]
        end;
      cmat^[i,j]:=sum
    end {j loop}
  End; {i loop}

End; {atbmul}


Procedure atvov (Var overfl:boolean; maxexp, n, nunit, output:integer;
                 amat: pMat; bvec, cvec: pVec);
{----------------------------------------------------------------------
!    FEB. 8 ,1991
!
!    THIS SUBROUTINE FINDS THE PRODUCT OF THE TRANSPOSE OF THE MATRIX A
!    AND THE VECTOR B WHERE EACH ENTRY IS CHECKED TO PREVENT OVERFLOWS.
!---------------------------------------------------------------------}
Label 40;
Var
    eps, SUM: double;
    i, j, k: integer;
Begin
  eps:=one/Power(ten, maxexp);
  overfl:=false;

  For i:=1 to n do
  begin
    sum:=zero;
    For j:=1 to n do
    begin
      IF LOG(ABS(amat^[j,i])+eps) + LOG(ABS(bvec^[j])+eps) > maxexp THEN
      begin
        overfl:=true;
        cvec^[i]:=SIGN(Power(ten,maxexp), amat^[j,i])*SIGN(one,bvec^[j]);
        IF (output > 2)  AND (NOT wrnsup) THEN
        begin
          Line0;
          Writeln(fp_out,'   *    WARNING: COMPONENT MATRIX-VECTOR PRODUCT SET TO ', cvec^[i]:12:3,'     *')
        end;
        GOTO 40
      end
    end;
    sum:=sum + amat^[j,i]*bvec^[j]
  end;  
  cvec^[i]:=sum;

40: End; {atvov}

Procedure avmul (nradec, nraact, ncdec, ncact:integer; amat:pMat; bvec, cvec:pVec);
{----------------------------------------------------------------------------------
!    FEB. 8, 1991
!
!    MATRIX-VECTOR MULTIPLICATION AB:=C
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4
!    EACH ROW OF MATRIX A IS SAVED AS A COLUMN BEFORE USE.
!
!    NRADEC IS 1ST DIM. OF AMAT; NRAACT IS ACTUAL LIMIT FOR 1ST INDEX
!    NCDEC IS COMMON DIMENSION OF AMAT & BVEC; NCACT IS ACTUAL LIMIT
!
!    I.E. NRADEC IS THE NUMBER OF ROWS OF A DECLARED
!         NCDEC IS THE COMMON DECLARED DIMENSION (COLUMNS OF A AND ROWS OF B)
!
!    MODIFICATION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES MACKINNON,
!    QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!---------------------------------------------------------------------------------}
Var
    i, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, ncc32r: integer;
    SUM: double;
Begin

{   FIND NUMBER OF GROUPS OF SIZE 32, 16... }

  ncc32:=ncact Div 32;
  ncc32r:=ncact - 32*ncc32;
  ncc16:=ncc32r Div 16;
  ncc16r:=ncc32r - 16*ncc16;
  ncc8:=ncc16r Div 8;
  ncc8r:=ncc16r - 8*ncc8;
  ncc4:=ncc8r Div 4;
  ncc4r:=ncc8r - 4*ncc4;

  For i:=1 to nraact do
  begin

{   FIND ENTRY FOR VECTOR C }

    sum:=zero;
    k:=0;
    IF ncc32 > 0 THEN
      For  kk:=1 to ncc32 do
      begin
        k:=k+32;
        sum:=sum + amat^[i,k-31]*bvec^[k-31] + amat^[i,k-30]*bvec^[k-30] +
          amat^[i,k-29]*bvec^[k-29] + amat^[i,k-28]*bvec^[k-28] +
          amat^[i,k-27]*bvec^[k-27] + amat^[i,k-26]*bvec^[k-26] +
          amat^[i,k-25]*bvec^[k-25] + amat^[i,k-24]*bvec^[k-24];
        sum:=sum + amat^[i,k-23]*bvec^[k-23] + amat^[i,k-22]*bvec^[k-22] +
          amat^[i,k-21]*bvec^[k-21] + amat^[i,k-20]*bvec^[k-20] +
          amat^[i,k-19]*bvec^[k-19] + amat^[i,k-18]*bvec^[k-18] +
          amat^[i,k-17]*bvec^[k-17] + amat^[i,k-16]*bvec^[k-16];
        sum:=sum + amat^[i,k-15]*bvec^[k-15] + amat^[i,k-14]*bvec^[k-14] +
          amat^[i,k-13]*bvec^[k-13] + amat^[i,k-12]*bvec^[k-12] +
          amat^[i,k-11]*bvec^[k-11] + amat^[i,k-10]*bvec^[k-10] +
          amat^[i,k-9]*bvec^[k-9] + amat^[i,k-8]*bvec^[k-8];
        sum:=sum + amat^[i,k-7]*bvec^[k-7] + amat^[i,k-6]*bvec^[k-6] +
          amat^[i,k-5]*bvec^[k-5] + amat^[i,k-4]*bvec^[k-4] + amat^[i,k-3]*bvec^[k-3] +
          amat^[i,k-2]*bvec^[k-2] + amat^[i,k-1]*bvec^[k-1] + amat^[i,k]*bvec^[k]
      end;
    IF ncc16 > 0 THEN
      For  kk:=1 to ncc16 do
      begin
        k:=k+16;
        sum:=sum + amat^[i,k-15]*bvec^[k-15] + amat^[i,k-14]*bvec^[k-14] +
          amat^[i,k-13]*bvec^[k-13] + amat^[i,k-12]*bvec^[k-12] +
          amat^[i,k-11]*bvec^[k-11] + amat^[i,k-10]*bvec^[k-10] +
          amat^[i,k-9]*bvec^[k-9] + amat^[i,k-8]*bvec^[k-8];
        sum:=sum + amat^[i,k-7]*bvec^[k-7] + amat^[i,k-6]*bvec^[k-6] +
          amat^[i,k-5]*bvec^[k-5] + amat^[i,k-4]*bvec^[k-4] + amat^[i,k-3]*bvec^[k-3] +
          amat^[i,k-2]*bvec^[k-2] + amat^[i,k-1]*bvec^[k-1] + amat^[i,k]*bvec^[k]
      end;
    IF ncc8 > 0 THEN
      For  kk:=1 to ncc8 do
      begin
        k:=k+8;
        sum:=sum + amat^[i,k-7]*bvec^[k-7] + amat^[i,k-6]*bvec^[k-6] +
          amat^[i,k-5]*bvec^[k-5] + amat^[i,k-4]*bvec^[k-4] + amat^[i,k-3]*bvec^[k-3] +
          amat^[i,k-2]*bvec^[k-2] + amat^[i,k-1]*bvec^[k-1] + amat^[i,k]*bvec^[k]
      end;
    IF ncc4 > 0 THEN
      For  kk:=1 to ncc4 do
      begin
        k:=k+4;
        sum:=sum + amat^[i,k-3]*bvec^[k-3] + amat^[i,k-2]*bvec^[k-2] +
          amat^[i,k-1]*bvec^[k-1] + amat^[i,k]*bvec^[k]
      end;
    IF ncc4r > 0 THEN
      For  kk:=1 to ncc4r do
      begin
        k:=k+1;
        sum:=sum + amat^[i,k]*bvec^[k]
      end;
    cvec^[i]:=sum
  end; {i loop}

End; {avmul}


Procedure bakdif (Var overfl:boolean; j, n:integer; Var deltaj:double; tempj:double;
                  fvec, fvecj1:pVec; jacfdm:pMat; xc:pVec);
{   FEB. 6, 1991 }
Var i:integer;
Begin
  deltaj:=tempj - xc^[j];
  fcn1(overfl, n, fvecj1, xc);
  IF (NOT overfl) THEN
    For i:=1 to n do
     jacfdm^[i,j]:=(fvec^[i] - fvecj1^[i])/deltaj
End; {bakdif}


Procedure bnddif (Var overfl:boolean; j, n:integer; epsmch:double; boundl, boundu,
                  fvecc, fvecj1:pVec; jacfdm:pMat; xc:pVec);
{-------------------------------------------------------------------------------
!   FEB. 15, 1991
!
!   FINITE DIFFERENCE CALCULATION WHEN THE BOUNDS FOR COMPONENT J ARE SO CLOSE
!   THAT NEITHER A FORWARD NOR BACKWARD DIFFERENCE CAN BE PERFORMED.
!------------------------------------------------------------------------------}
Var
    i: integer;
    eps3q: double;
    wv3:pVec;
Begin
  eps3q:=Power1(epsmch,0.75);
  New(wv3);
{ STORE CURRENT }
  For i:=1 to n do wv3^[i]:=fvecc^[i];
  xc^[j]:=boundu^[j];
  fcn1 (overfl, n, fvecj1, xc);
  IF (NOT overfl) THEN
  begin
    xc^[j]:=boundl^[j];
    fcn1 (overfl, n, fvecc, xc);
    IF (NOT overfl) THEN
      For i:=1 to n do
{     ENSURE THAT THE JACOBIAN CALCULATION ISN'T JUST NOISE }
        IF fvecj1^[i]-fvecc^[i] > eps3q THEN
          jacfdm^[i,j]:=(fvecj1^[i]-fvecc^[i])/(boundu^[j]-boundl^[j])
        ELSE
          jacfdm^[i,j]:=zero
  end;
  For i:=1 to n do fvecc^[i]:=wv3^[i];
  Dispose(wv3)
End; {bnddif}

Procedure broyfa (Var overch, overfl:boolean; sclfch, sclxch:boolean;  maxexp, n,
                  nunit, output:integer; epsmch:double; a:pMat;  delf, fvec,
                  fvecc:pVec; jac:pMat; rdiag, s, scalef, scalex, t, w, xc,
                  xplus:pVec);
{--------------------------------------------------------------------------
!    FEB. 23, 1992
!
!    THE BROYDEN QUASI-NEWTON METHOD IS APPLIED TO THE FACTORED
!    FORM OF THE JACOBIAN.
!
!    NOTE: T AND W ARE TEMPORARY WORKING VECTORS ONLY.
!
!    THE UPDATE OCCURS ONLY IF A SIGNIFICANT CHANGE IN THE JACOBIAN
!    WOULD RESULT,  I.E., NOT ALL THE VALUES IN VECTOR W ARE LESS THAN
!    THE THRESHOLD IN MAGNITUDE.  IF THE VECTOR W IS ESSENTIALLY ZERO THEN
!    THE LOGICAL VARIABLE SKIPUP REMAINS SET AT TRUE.
!-------------------------------------------------------------------------}
Label return;
Var
    denom, eps, SUM, sqrtep: double;
    skipup: boolean;
    i,k: integer; tmp:pVec;
    tmp1,tmp2: pMat;
Begin
  New(tmp); New(tmp1); New(tmp2);
  overfl:=false;
  eps:=one/Power(ten,maxexp);
  sqrtep:=SQRT(epsmch);

  For  i:=1 to n do
  begin
    a^[i,i]:=rdiag^[i];
    s^[i]:=xplus^[i] - xc^[i]
  end;

{   R IS NOW IN THE UPPER TRIANGLE OF A }

  skipup:=true;

{       THE BROYDEN UPDATE IS CONDENSED INTO THE FORM

        A(NEW) := A(OLD) + T S^
 
        THE PRODUCT A*S IS FORMED IN TWO STAGES AS R IS IN THE UPPER
        TRIANGLE OF MATRIX A AND Q^ IS IN JAC.

        FIRST MULTIPLY R*S (A IS CONSIDERED UPPER TRIANGULAR)
}
  uvmul (n, n, n, n, a, s, t);

{       NOTE: THIS T IS TEMPORARY - IT IS THE T FROM BELOW WHICH
              IS SENT TO SUBROUTINE QRUPDA.  }

  For  i:=1 to n do
  begin
    For k:=1 to n do tmp^[k]:=jac^[k,i];
    innerp (overch, overfl, maxexp, n, n, n, nunit, output, sum, tmp, t);
    w^[i]:=scalef^[i]*(fvec^[i]-fvecc^[i]) - sum;

{   TEST TO ENSURE VECTOR W IS NONZERO.  ANY VALUE GREATER
    THAN THE THRESHOLD WILL SET SKIPUP TO FALSE.  }

    IF ABS(w^[i]) > sqrtep*scalef^[i]*(ABS(fvec^[i]) + ABS(fvecc^[i])) THEN
      skipup:= false
    ELSE
      w^[i]:=zero
  end;

{   IF W(I):=0 FOR ALL I THEN THE UPDATE IS SKIPPED }

  IF (NOT skipup) THEN
  begin
{   T:=Q^W; Q^ IS IN JAC }

    avmul (n, n, n, n, jac, w, t);
    IF (sclxch) THEN
      For k:=1 to n do w^[k]:=s^[k]*scalex^[k]
    ELSE
    begin
      For k:=1 to n do
      begin
        tmp1^[k,1]:=s^[k]; tmp2^[k,1]:=w^[k]
      end;
      matcop (n, n, 1, 1, n, 1, tmp1, tmp2);
      For k:=1 to n do
      begin
        s^[k]:=tmp1^[k,1]; w^[k]:=tmp2^[k,1]
      end;
    end;

    twonrm (overfl, maxexp, n, epsmch, denom, w);

{   IF OVERFLOW WOULD OCCUR MAKE NO CHANGE TO JACOBIAN }

    IF (overfl) OR (LOG(denom+eps) > 0.5*maxexp) THEN
    begin
      IF (output > 3) THEN
      begin
        Line0;
        Msg('    WARNING: JACOBIAN NOT UPDATED TO AVOID OVERFLOW IN DENOMINATOR OF');
        Msg('    BROYDEN UPDATE.')
       end;
      goto return
    end
    ELSE
      denom:=denom*denom;

{   IF DENOM IS ZERO AVOID DIVIDE BY ZERO AND CONTINUE WITH SAME JACOBIAN }

    IF denom = zero  then goto return;

{   THE SCALED VERSION OF S REPLACES THE ORIGINAL BEFORE
    BEING SENT TO QRUPDA.   }

    For k:=1 to n do
      s^[k]:=s^[k]*scalex^[k]*scalex^[k]/denom;

{   UPDATE THE QR DECOMPOSITION USING A SERIES OF GIVENS ROTATIONS. }

    qrupda (overfl, maxexp, n, epsmch, a, jac, t, s);

{   RESET RDIAG AS DIAGONAL OF CURRENT R WHICH IS IN
    THE UPPER TRIANGE OF A.  }

    For i:=1 to n do rdiag^[i]:=a^[i,i]
  end;

{   UPDATE THE GRADIENT VECTOR, DELF.  THE NEW Q^ IS IN JAC.

    DELF := (QR)^F := R^Q^F := R^JAC F        }

  IF (sclfch) THEN
    For k:=1 to n do w^[k]:=fvec^[k]*scalef^[k]
  ELSE
  begin
    For k:=1 to n do
    begin
      tmp1^[k,1]:=fvec^[k]; tmp2^[k,1]:=w^[k]
    end;
    matcop (n, n, 1, 1, n, 1, tmp1, tmp2);
    For k:=1 to n do
    begin
      fvec^[k]:=tmp1^[k,1]; w^[k]:=tmp2^[k,1]
    end
  end;

  avmul (n, n, n, n, jac, w, t);
  For k:=1 to n do
  begin
    tmp1^[k,1]:=t^[k]; tmp2^[k,1]:=delf^[k]
  end;
  utbmul (n, n, 1, 1, n, n, a, tmp1, tmp2);
  For k:=1 to n do
  begin
    t^[k]:=tmp1^[k,1]; delf^[k]:=tmp2^[k,1]
  end;

  Dispose(tmp);  Dispose(tmp1);  Dispose(tmp2);
  
return: End; {broyfa}



Procedure broyun (var overfl:boolean; maxexp, n, nunit, output:integer;
                  var epsmch:double; fvec, fvecc:pVec; jac:pMat;
                  scalex, xc, xplus:pVec);
{---------------------------------------------------------------------
!    FEB. 23, 1992
!
!    UPDATE THE JACOBIAN USING BROYDEN'S METHOD.
!--------------------------------------------------------------------}
Label return;
Var
    i, j, k:integer;
    denom, eps, sqrtep, SUM, tempi:double;
    wv1,tmp1,tmp2:pVec;
Begin
  New(wv1); New(tmp1); New(tmp2);
  eps:=one/Power(ten,maxexp);
  sqrtep:=SQRT(epsmch);

  For k:=1 to n do
    wv1^[k]:=(xplus^[k] - xc^[k])*scalex^[k];

  twonrm (overfl, maxexp, n, epsmch, denom, wv1);

{ IF OVERFLOW WOULD OCCUR MAKE NO CHANGE TO JACOBIAN }

  IF (overfl) OR (LOG(denom+eps) > 0.5*maxexp) THEN
  begin
    IF (output > 3) THEN
    begin
      Line0;
      Msg('    WARNING: JACOBIAN NOT UPDATED TO AVOID OVERFLOW IN DENOMINATOR OF');
      Msg('    BROYDEN UPDATE.')
    end;
    goto return
  end
  ELSE
    denom:=denom*denom;

{ IF DENOM IS ZERO, AVOID OVERFLOW, CONTINUE WITH SAME JACOBIAN }

  IF denom = zero then goto return;

{ UPDATE JACOBIAN BY ROWS }

  For i:=1 to n do
  begin
    For k:=1 to n do
    begin
      tmp1^[k]:=jac^[i,k]; tmp2^[k]:=xplus^[k] - xc^[k]
    end;
    sum := DOT_PRODUCT(n, tmp1, tmp2);
    tempi:=fvec^[i] - fvecc^[i] - sum;

{ CHECK TO ENSURE THAT SOME MEANINGFUL CHANGE IS BEING MADE
  TO THE APPROXIMATE JACOBIAN; IF NOT, SKIP UPDATING ROW I. }

    IF ABS(tempi) >= sqrtep*(ABS(fvec^[i]) + ABS(fvecc^[i])) THEN
    begin
      tempi:=tempi/denom;
      For j:=1 to n do
        jac^[i,j]:=jac^[i,j] + tempi*(xplus^[j]-xc^[j])*scalex^[j]*scalex^[j]
    end
  end;
Return: Dispose(wv1)
End; {broyun}



Procedure cholde (n:integer; var maxadd, maxffl:double; sqrtep:double;
                  h, l: pMat);
{-------------------------------------------------------------------------
!    FEB. 23, 1992
!
!    THIS SUBROUTINE FINDS THE CHOLESKY DECOMPOSITION OF THE MATRIX, H,
!    AND RETURNS IT IN THE LOWER TRIANGLE OF MATRIX, L.
!------------------------------------------------------------------------}
Var
    minl, minl2, minljj, SUM:double;
    i,j,k: integer; tmp1,tmp2: pVec;
Begin
    New(tmp1); New(tmp2);

{   minl:=SQRT(sqrtep)*maxffl

    MAXFFL EQUALS 0 WHEN THE MATRIX IS KNOWN TO BE POSITIVE DEFINITE }

  IF maxffl = zero THEN
  begin
{   FIND SQUARE ROOT OF LARGEST MAGNITUDE DIAGONAL ELEMENT
    AND SET MINL2.    }
    For i:=1 to n do maxffl:=MAX(maxffl, ABS(h^[i,i]));
    maxffl:=SQRT(maxffl);
    minl2:=sqrtep*maxffl
  end;

{  MAXADD CONTAINS THE MAXIMUM AMOUNT WHICH IS IMPLICITLY ADDED
   TO ANY DIAGONAL ELEMENT OF MATRIX H. }

  maxadd:=zero;
  For j:=1 to n do
  begin
    For k:=1 to j-1 do tmp1^[k]:=l^[j,k];
    sum := DOT_PRODUCT(j-1, tmp1, tmp1);
    l^[j,j]:=h^[j,j]-sum;
    minljj:=zero;
    For i:=j+1 to n do
    begin
      For k:=1 to j-1 do
      begin
        tmp1^[k]:=l^[i,k];
        tmp2^[k]:=l^[j,k]
      end;
      sum := DOT_PRODUCT(j-1, tmp1, tmp2);
      l^[i,j]:=h^[j,i] - sum;
      minljj:=MAX(minljj, ABS(l^[i,j]))
    end;
    minljj:=MAX(minljj/maxffl, minl);
    IF l^[j,j] > minljj*minljj THEN

{   NORMAL CHOLESKY DECOMPOSITION }
      l^[j,j]:=SQRT(l^[j,j])
    ELSE
{   IMPLICITLY PERTURB DIAGONAL OF H }
    begin
      IF (minljj < minl2) THEN minljj:=minl2;
      maxadd:=MAX(maxadd, minljj*minljj-l^[j,j]);
      l^[j,j]:=minljj
    end;
    For k:=j+1 to n do l^[j+k,j]:=l^[j+k,j]/l^[j,j]
  end;
  Dispose(tmp1); Dispose(tmp2)
End; {cholde}

Procedure lsolv (overch:boolean; var overfl:boolean; maxexp, n, nunit,
                 output:integer; l:pMat; b, rhs:pVec); Forward;
Procedure ltsolv (overch:boolean; var overfl:boolean; maxexp, n, nunit,
                  output:integer; l:pMat; y, b:pVec); Forward;

Procedure chsolv (overch:boolean; var overfl:boolean; maxexp, n, nunit,
                  output:integer; l:pMat; rhs, s:pVec);
{---------------------------------------------------------------------
!
!    FEB. 14, 1991
!
!    THIS SUBROUTINE USES FORWARD/BACKWARD SUBSTITUTION TO SOLVE THE
!    SYSTEM OF LINEAR EQUATIONS:
!
!         (LL^)S:=RHS
!---------------------------------------------------------------------}
Var
    wv2:pVec;
Begin
  New(wv2);
  lsolv (overch, overfl, maxexp, n, nunit, output, l, wv2, rhs);
  ltsolv (overch, overfl, maxexp, n, nunit, output, l, s, wv2);
  Dispose(wv2)
End; {chsolv}

Procedure condno (overch:boolean; var overfl:boolean; maxexp, n, nunit,
                  output:integer; var connum:double; a:pMat; rdiag:pVec);
{-----------------------------------------------------------------------------
! N.B. Arguments P, PM & Q have been removed.
!
!    FEB. 14, 1991
!
!    THIS SUBROUTINE ESTIMATES THE CONDITION NUMBER OF A QR-DECOMPOSED
!    MATRIX USING THE METHOD OF CLINE, MOLER, STEWART AND WILKINSON
!    (SIAM J. N.A. 16 P368 (1979) ).
!
!    IF A POTENTIAL OVERFLOW IS DETECTED AT ANY POINT THEN A CONDITION NUMBER
!    EQUIVALENT TO THAT OF A SINGULAR MATRIX IS ASSIGNED BY THE CALLING
!    SUBROUTINE.
!----------------------------------------------------------------------------}
Label return;
Var 
    i,j,k: integer;
    eps, qm, qnorm, qp, temp, tempm: double;
    p, pm, q: pVec;
Begin
  New(p); New(pm); New(q);
  overfl:=false;
  eps:=one/Power(ten, maxexp);

  connum:=ABS(rdiag^[1]);
  For j:=2 to n do
  begin
    temp:=zero;
    For i:=1 to j-1 do
    begin
      IF (overch) THEN
        IF ABS(a^[i,j]) > Power(ten,maxexp-1) THEN
        begin
          overfl:=true;
          goto RETURN
        end;
      temp:=temp + ABS(a^[i,j])
    end;
    temp:=temp + ABS(rdiag^[j]);
    connum:=MAX(connum,temp)
  end;
  q^[1]:=one/rdiag^[1];
  For i:=2 to n do
  begin
    IF (overch) THEN
      IF LOG(ABS(q^[1])+eps) + LOG(ABS(a^[1,i]) + eps) > maxexp THEN
      begin
        overfl:=true;
        goto RETURN
      end;
    p^[i]:=a^[1,i]*q^[1]
  end;
  For j:=2 to n do
  begin
    IF (overch) THEN
      IF LOG(ABS(p^[j])+eps) - LOG(ABS(rdiag^[j]) + eps) > maxexp THEN
      begin
        overfl:=true;
        goto RETURN
      end;
    qp:=(one-p^[j])/rdiag^[j];
    qm:=(-one-p^[j])/rdiag^[j];
    temp:=ABS(qp);
    tempm:=ABS(qm);
    For i:=j+1 to n do
    begin
      IF (overch) THEN
        IF LOG(ABS(a^[j,i])+eps) + LOG(ABS(qm)+eps) > maxexp THEN
        begin
          overfl:=true;
          goto RETURN
        end;
      pm^[i]:=p^[i] + a^[j,i]*qm;
      IF (overch) THEN
        IF LOG(ABS(pm^[i])+eps) - LOG(ABS(rdiag^[i])+eps) > maxexp THEN
        begin
          overfl:=true;
          goto RETURN
        end;
      tempm:=tempm + (ABS(pm^[i])/ABS(rdiag^[i]));
      IF (overch) THEN
        IF tempm > Power(ten,maxexp-1) THEN
        begin
          overfl:=true;
          goto RETURN
        end;
      IF (overch) THEN
        IF LOG(ABS(a^[j,i])+eps) + LOG(ABS(qp)+eps) > maxexp THEN
        begin
          overfl:=true;
          goto RETURN
        end;
      p^[i]:=p^[i] + a^[j,i]*qp;
      IF (overch) THEN
        IF LOG(ABS(p^[i])+eps) - LOG(ABS(rdiag^[i])+eps) > maxexp THEN
        begin
          overfl:=true;
          goto RETURN
        end;
      temp:=temp + (ABS(p^[i])/ABS(rdiag^[i]));
      IF (overch) THEN
        IF temp > Power(ten, maxexp-1) THEN
        begin
          overfl:=true;
          goto RETURN
        end
    end; {i loop}
    IF temp >= tempm THEN
      q^[j]:=qp
    ELSE
    begin
      q^[j]:=qm;
      For k:=1 to n do p^[j+k]:=pm^[j+k]
    end
  end; {j loop}
  qnorm:=zero;
  For j:=1 to n do
  begin
    qnorm:=qnorm + ABS(q^[j]);
    IF (overch) THEN
      IF qnorm > Power(ten, maxexp-1) THEN
      begin
        overfl:=true;
        goto RETURN
      end
  end;
  IF LOG(connum) - LOG(qnorm) > maxexp THEN
  begin
    overfl:=true;
    goto RETURN
  end;
  connum:=connum/qnorm;
  rsolv (overch, overfl, maxexp, n, nunit, output, a, rdiag, q);
  IF (overfl) then goto RETURN;
  qnorm:=zero;
  For j:=1 to n do
  begin
    qnorm:=qnorm + ABS(q^[j]);
    IF (overch) THEN
      IF qnorm > Power(ten, maxexp-1) THEN
      begin
        overfl:=true;
        goto RETURN
      end
  end;
  connum:=connum*qnorm;

Return: Dispose(p); Dispose(pm); Dispose(q)
End; {condno}

Procedure delcau (cauchy, overch:boolean; var overfl:boolean; itnum,
                  maxexp, n, nunit, output:integer; var beta, caulen,
                  delta, epsmch, maxstp, newlen, sqrtz:double;
                  a: pMat; delf, scalex: pVec);
{--------------------------------------------------------------------
!
!    FEB. 23, 1992
!
!    THIS SUBROUTINE ESTABLISHES AN INITIAL TRUST REGION, DELTA,
!    IF ONE IS NOT SPECIFIED BY THE USER AND FINDS THE LENGTH OF
!    THE SCALED CAUCHY STEP, CAULEN, AT EACH STEP IF THE DOUBLE
!    DOGLEG OPTION IS BEING USED.
!
!    THE USER HAS TWO CHOICES FOR THE INITIAL TRUST REGION:
!       1)  BASED ON THE SCALED CAUCHY STEP
!       2)  BASED ON THE SCALED NEWTON STEP
!-------------------------------------------------------------------}
Label 120;
Const three = 3.0;
Var
    eps, temp:double;
    wv1:pVec;
    i,j,k:integer;
Begin
  New(wv1);
  overfl:=false;
  eps:=one/Power(ten, maxexp);

{ IF DELTA IS NOT GIVEN EVALUATE IT USING EITHER THE CAUCHY
  STEP OR THE NEWTON STEP AS SPECIFIED BY THE USER.

  THE SCALED CAUCHY LENGTH, CAULEN, IS REQUIRED IN TWO CASES.
   1) WHEN SELECTED AS THE CRITERION FOR THE INITIAL DELTA
   2) IN THE DOUBLE DOGLEG STEP REGARDLESS OF (1)          }

  IF output > 3 THEN
  begin
    Line0;
    Msg('    DETERMINATION OF SCALED CAUCHY STEP LENGTH, CAULEN');
  end;

{ FIND FACTOR WHICH GIVES CAUCHY POINT WHEN MULTPLYING
  STEEPEST DESCENT DIRECTION, DELF.

     CAULEN = ZETA^1.5/BETA
            =  SQRTZ^3/BETA        }
  For k:=1 to n do wv1^[k]:=delf^[k]/scalex^[k];
  twonrm (overfl, maxexp, n, epsmch, sqrtz, wv1);
  IF (overfl) THEN
  begin
    caulen:=Power(ten,maxexp);
    IF output > 4 THEN
    begin
      Line0;
      Writeln(fp_out,'   *       ZETA SET TO ', sqrtz:11:3, ' TO AVOID OVERFLOW                       *');
      Writeln(fp_out,'   *       SCALED CAUCHY LENGTH, CAULEN SET TO ',caulen:9:2,' TO AVOID OVERFLOW *');
      IF itnum = 1 THEN
      begin
        Line0;
        Msg('       THE PROBLEM SHOULD BE RESCALED OR A NEW STARTING POINT CHOSEN');
        Msg('       EXECUTION CONTINUES WITH SUBSTITUTIONS AS LISTED ABOVE.')
      end
    end
  end
  ELSE
    IF output > 4 THEN
    begin
      Line0;
      Writeln(fp_out,'   *       SQUARE ROOT OF ZETA, SQRTZ: ', sqrtz:12:3, '            *')
    end;

{ NOTE: THE LOWER TRIANGLE OF MATRIX A NOW CONTAINS THE
  TRANSPOSE OF R WHERE A:=QR.  }

  beta:=zero;
  For  i:=1 to n do
  begin
    temp:=zero;
    For j:=i to n do
    begin
      IF (overch) THEN
        IF LOG(ABS(a^[j,i])+eps) + LOG(ABS(delf^[j])+eps) > maxexp THEN
        begin
          caulen:=SQRT(epsmch);
          GOTO 120
        end;
      temp:=temp + a^[j,i]*delf^[j]/(scalex^[j]*scalex^[j])
    end;
    beta:=beta + temp*temp
  end;
  IF output > 4 THEN
  begin
    Line0;
    Writeln(fp_out,'   *       BETA:',BETA:11:3,'      NOTE: CAULEN = ZETA^1.5/BETA             *');
    Line0
  end;

{   AVOID OVERFLOWS IN FINDING CAULEN }

  IF (three*LOG(sqrtz+eps) - LOG(beta+eps) < maxexp) AND (NOT overfl) AND (beta <> zero) THEN
  begin
{   NORMAL DETERMINATION }

    caulen:=sqrtz*sqrtz*sqrtz/beta;

{   THIS STEP AVOIDS DIVIDE BY ZERO IN DOGLEG IN THE (RARE) CASE
    WHERE DELF(I):=0 FOR ALL I BUT THE POINT IS NOT YET A SOLUTION -
    MOST LIKELY A BAD STARTING ESTIMATE.  }

    caulen:=MAX(caulen,one/Power(ten, maxexp))
  end
  ELSE
{   SUBSTITUTION TO AVOID OVERFLOW }
    caulen:=Power(ten, maxexp);

120:IF output > 3 THEN
  begin
    Line0;
    Writeln(fp_out,'   *       SCALED CAUCHY LENGTH, CAULEN: ', caulen:12:3,'        *')
  end;

{   ESTABLISH INITIAL TRUST REGION IF NEEDED }

  IF (delta < zero) THEN
  begin
{   USE DISTANCE TO CAUCHY POINT OR LENGTH OF NEWTON STEP }

    IF (cauchy) THEN
      delta:=MIN(caulen,maxstp)
    ELSE
      delta:=MIN(newlen,maxstp);

    IF output > 3 THEN
    begin
      Line0;
      Writeln(fp_out,'   *       INITIAL TRUST REGION SIZE, DELTA: ', delta:12:3,'       *')
    end
  end;
  Dispose(wv1)
End; {delcau}


Procedure deufls (var abort:boolean; deuflh, geoms:boolean; var overch,overfl,qnfail:boolean;
                   qrsing:boolean; var restrt:boolean; sclfch:boolean; var sclxch:boolean;
                   var acpcod:integer; acptcr, contyp, itnum,jupdm, maxexp, maxlin, n:integer;
                   var nfunc:integer; nunit, output, qnupdm:integer; var stopcr:integer;
                   alpha, confac, delfts, epsmch, fcnmax:double; var fcnnew:double;
                   fcnold:double; var lambda:double; newmax:double; var sbrnrm:double; sigma:double;
                   a, astore:pMat; boundl, boundu, delf, fvec, hhpi:pVec; jac:pMat; rdiag, rhs, s,
                   sbar, scalef, scalex, sn, xc, xplus: pVec);
{----------------------------------------------------------------------------
!    FEB. 23, 1992
!
!    THIS SUBROUTINE CONDUCTS A LINE SEARCH IN THE NEWTON DIRECTION
!    IF NO CONSTRAINTS ARE VIOLATED.  IF THE FIRST TRIAL IS A FAILURE
!    THERE ARE TWO TYPES OF LINE SEARCH AVAILABLE:
!      1)  REDUCE THE RELAXATION FACTOR, LAMBDA, TO SIGMA*LAMBDA
!          WHERE SIGMA IS USER-SPECIFIED (GEOMETRIC Line0 SEARCH)
!      2)  AT THE FIRST STEP MINIMIZE A QUADRATIC THROUGH THE OBJECTIVE
!          FUNCTION AT THE CURRENT POINT AND TRIAL ESTIMATE (WHICH MUST BE A
!          FAILURE) WITH INITIAL SLOPE DELFTS.  AT SUBSEQUENT STEPS MINIMIZE
!          A CUBIC THROUGH THE OBJECTIVE FUNCTION AT THE TWO MOST RECENT
!          FAILURES AND THE CURRENT POINT, AGAIN USING THE INITIAL SLOPE,
!          DELFTS.
!
!    CONVIO  INDICATES A CONSTRAINT VIOLATION BY ONE OR MORE COMPONENTS
!    FRSTST  INDICATES THE FIRST STEP IN THE LINE SEARCH.
!
!    RATIO   RATIO OF PROPOSED STEP LENGTH IN (I)TH DIRECTION
!            TO DISTANCE FROM (I)TH COMPONENT TO BOUNDARY VIOLATED
!    RATIOM  MINIMUM OF RATIOS FOR ALL CONSTAINTS VIOLATED
!---------------------------------------------------------------------------}
Const  point1 = 0.1; point5 = 0.5; three = 3.0;
Var
    acubic, bcubic, disc, dlftsm, eps, factor, fplpre, lampre, lamtmp,
    ratio, ratiom: double;
    wv2:pVec;
    i, ii, j, k:integer;
    convio, frstst:boolean;
Begin
  New(wv2);
  frstst:=true;
  overfl:=false;
  abort := FALSE;
  eps:=one/Power(ten, maxexp);
  {lambda:=one;  added by JPM}

  For k:=1 to maxlin do
  begin
    ratiom:=one;
    convio:=false;

{   FIND TRIAL POINT AND CHECK IF CONSTRAINTS VIOLATED (IF
    CONTYP IS NOT EQUAL TO ZERO).  }

    For i:=1 to n do
    begin
{   NOTE: WV2 MARKS VIOLATIONS.  WV2(I) CHANGES TO 1 FOR LOWER BOUND
          VIOLATIONS AND TO 2 FOR UPPER BOUND VIOLATIONS.
          CONSTRAINT VIOLATIONS CAN ONLY OCCUR AT THE FIRST STEP.  }
      wv2^[i]:=-one;
      xplus^[i]:=xc^[i] + lambda*sn^[i];
      IF (contyp > 0) AND (frstst) THEN
        IF xplus^[i] < boundl^[i] THEN
        begin
          convio:=true;
          wv2^[i]:=one
        end
        ELSE IF xplus^[i] > boundu^[i] THEN
        begin
          convio:=true;
          wv2^[i]:=two
        end
        ELSE
          wv2^[i]:=-one
    end; {i loop}

{     IF CONSTRAINTS ARE VIOLATED FIRST REDUCE THE STEP SIZES FOR THE
      VIOLATING COMPONENTS TO OBTAIN A FEASIBLE POINT.
      IF THE DIRECTION TO THIS MODIFIED POINT IS NOT A DESCENT DIRECTION
      OR IF THE MODIFIED STEP DOES NOT LEAD TO AN ACCEPTABLE POINT THEN
      RETURN TO THE NEWTON DIRECTION AND START A LINE SEARCH AT A FEASIBLE
      POINT WHERE THE COMPONENT WHICH HAS THE SMALLEST VALUE OF RATIO
      (DEFINED BELOW) IS TAKEN TO "CONFAC" OF THE DISTANCE TO THE BOUNDARY.
      DEFAULT VALUE OF CONFAC IS 0.95.      }

    IF (convio) THEN
    begin
      IF output > 3 THEN
      begin
        Line0;
        Writeln(fp_out,'   *          LINE SEARCH STEP:', k:3, '                                         *');
        Line0;
        Writeln(fp_out,'   *          LAMBDA FOR ATTEMPTED STEP: ',lambda:12:3,'                       *');
        Msg('          CONSTRAINT VIOLATED');
        Line0;
        Msg('          TRIAL ESTIMATES (VIOLATIONS MARKED).');
        Line0;
        For i:=1 to n do
          IF wv2^[i] > zero THEN
            Writeln(fp_out,'   *             XPLUS(',i:3,') = ',xplus^[i]:12:3,'*                                *')
          ELSE
            Writeln(fp_out,'   *             XPLUS(',i:3,') = ',xplus^[i]:12:3,'                                 *')
      end;

      For  i:=1 to n do
        IF wv2^[i] > zero THEN
        begin
{         FIND RATIO FOR THIS VIOLATING COMPONENT }

          IF wv2^[i] = one THEN
            ratio:=-(xc^[i]-boundl^[i])/(xplus^[i]-xc^[i])
          ELSE IF wv2^[i] = two THEN
            ratio:=(boundu^[i]-xc^[i])/(xplus^[i]-xc^[i]);

{         NOTE: THIS Line0 IS FOR OUTPUT PURPOSES ONLY }

          wv2^[i]:=ratio;

          ratiom:=MIN(ratiom,ratio);
          IF ratio > point5 THEN
            s^[i]:=confac*ratio*lambda*sn^[i]
          ELSE
{         WITHIN BUFFER ZONE - ONLY TAKE 1/2
          OF THE STEP YOU WOULD TAKE OTHERWISE.  }
            s^[i]:=point5*confac*ratio*lambda*sn^[i];

{        ESTABLISH MODIFIED TRIAL POINT }

          xplus^[i]:=xc^[i] + s^[i]
        end
        ELSE
        begin
{         FOR NONVIOLATORS XPLUS REMAINS UNCHANGED BUT THE COMPONENT
          OF S IS LOADED TO CHECK THE DIRECTIONAL DERIVATIVE.  }

          s^[i]:=lambda*sn^[i]
        end;

      IF output > 3 THEN
      begin
        Line0;
        Msg('       NEW S AND XPLUS VECTORS (WITH RATIOS FOR VIOLATIONS)');
        Line0;
        Msg('       NOTE: RATIOS ARE RATIO OF LENGTH TO BOUNDARY FROM CURRENT');
        Msg('             X VECTOR TO MAGNITUDE OF CORRESPONDING PROPOSED STEP.');
        Line0;
        For i:=1 to n do
          IF wv2^[i] < zero THEN
            Writeln(fp_out,'   *       S(',i:3,') = ',s^[i]:12:3,'    XPLUS(',i:3,') = ',xplus^[i]:12:3,'              *')
          ELSE
            Writeln(fp_out,'   *       S(',i:3,') = ',s^[i]:12:3,'    XPLUS(',i:3,') = ',xplus^[i]:12:3,' ',wv2^[i]:11:3,' *');
        Line0;
        Writeln(fp_out,'   *       MINIMUM OF RATIOS, RATIOM: ',ratiom:12:3,'                        *')
      end;

{   CHECK DIRECTIONAL DERIVATIVE FOR MODIFIED POINT, DLFTSM }

      innerp (overch, overfl, maxexp, n, n, n, nunit, output, dlftsm, delf, s);
      overfl:=false;
      IF output > 3 THEN
      begin
        Line0;
        Writeln(fp_out,'   *       INNER PRODUCT OF DELF AND S FOR MODIFIED S: ', dlftsm:12:3,'       *')
      end;

{   IF INNER PRODUCT IS POSITIVE RETURN TO NEWTON DIRECTION
    AND CONDUCT A LINE SEARCH WITHIN THE FEASIBLE REGION.   }

      IF dlftsm > zero THEN
      begin
        IF output > 3 THEN
        begin
          Msg('     DELFTS > 0  START LINE SEARCH AT LAMBDA = CONFAC*LAMBDA*RATIOM');
          Msg('     NOTE: NO TRIAL POINT WAS EVALUATED AT THIS STEP OF LINE SEARCH.')
        end;

{   THE STARTING POINT IS SET AT JUST INSIDE THE MOST
    VIOLATED BOUNDARY.  }

        lambda:=confac*ratiom*lambda

{   LAMBDA IS ALREADY SET - SKIP NORMAL PROCEDURE }

        {CYCLE}
      end

    end; {if convio}

{   NO CONSTRAINTS VIOLATED - EVALUATE RESIDUAL VECTOR AT NEW POINT }

    fcn1 (overfl, n, fvec, xplus);

    Inc(nfunc);

{   CHECK FOR OVERFLOW IN FUNCTION VECTOR EVALUATION.
    IF SO, REDUCE STEP LENGTH AND CONTINUE LINE SEARCH. }

    IF (overfl) THEN
    begin 
      IF output > 3 THEN
        
        Msg('       OVERFLOW IN FUNCTION VECTOR - STEP LENGTH REDUCED.');

{       FORCE STEP TO BE WITHIN CONSTRAINTS - DON'T CALL
        THIS THE FIRST STEP, I.E. FRSTST STAYS AT TRUE.  }

      IF (convio) THEN
        lambda:=ratiom*confac*lambda
      ELSE
        lambda:=sigma*lambda;

{     LAMBDA IS ALREADY SET - SKIP NORMAL PROCEDURE }

      {CYCLE}
    end;

{     EVALUATE (POSSIBLY SCALED) OBJECTIVE FUNCTION AT NEW POINT }

    fcnevl (overfl, maxexp, n, nunit, output, epsmch, fcnnew, fvec, scalef);
    
    IF (overfl) THEN
    begin
      IF output > 3 THEN
        Msg('       OVERFLOW IN OBJECTIVE FUNCTION - STEP LENGTH REDUCED.');

{       FORCE STEP TO BE WITHIN CONSTRAINTS - DON'T CALL
        THIS THE FIRST STEP, I.E. FRSTST STAYS AT TRUE.  }

      IF (convio) THEN
        lambda:=ratiom*confac*lambda
      ELSE
        lambda:=sigma*lambda;
      {CYCLE}
    end;

{       IF DEUFLHARD'S METHOD IS BEING USED FOR EITHER RELAXATION FACTOR
        INITIALIZATION OR THE SECOND ACCEPTANCE CRITERION THEN EVALUATE SBAR.
        EVALUATION METHOD DEPENDS UPON WHETHER THE JACOBIAN WAS PERTURBED
        IN THE SOLUTION OF THE LINEAR SYSTEM.
        LOGICAL VARIABLE QRSING IS TRUE IF PERTURBATION TOOK PLACE.   }

    IF (deuflh) OR (acptcr = 12) THEN
    begin
      IF (qrsing) THEN
      begin
{       FORM -J^F AS RIGHT HAND SIDE - METHOD DEPENDS ON WHETHER
        QNUPDM EQUALS 0 OR 1 IF A QUASI-NEWTON UPDATE IS BEING USED.
        IF JUPDM IS 0 THEN THE NEWTON STEP HAS BEEN FOUND IN SUBROUTINE
        NSTPUN.    }

        IF (jupdm = 0) OR (qnupdm = 0) THEN
        begin
{         UNSCALED JACOBIAN IN MATRIX JAC }
          For i:=1 to n do
            IF (sclfch) THEN
              wv2^[i]:=-fvec^[i]*scalef^[i]*scalef^[i]
            ELSE
              wv2^[i]:=-fvec^[i];
          mmulv (n, jac, wv2, rhs)
        end
        ELSE
        begin
{         R IN UPPER TRIANGLE OF A PLUS RDIAG AND Q^ IN JAC
          - FROM QR DECOMPOSITION OF SCALED JACOBIAN.  }
          For i:=1 to n do
          begin
            wv2^[i]:=zero;
            For j:=1 to n do wv2^[i]:=wv2^[i] - jac^[i,j]*fvec^[j]*scalef^[j]
          end;
          rhs^[1]:=rdiag^[1]*wv2^[1];
          For j:=2 to n do
          begin
            rhs^[j]:=zero;
            For i:=1 to j-1 do rhs^[j]:=rhs^[j] + a^[i,j]*wv2^[i];
            rhs^[j]:=rhs^[j] + rdiag^[j]*wv2^[j]
          end
        end;
        chsolv(overch, overfl, maxexp, n, nunit, output, a, rhs, sbar)
      end
      ELSE
      begin
{       RIGHT HAND SIDE IS -FVEC }
        IF (qnupdm = 0) OR (jupdm = 0) THEN
        begin
{       QR DECOMPOSITION OF SCALED JACOBIAN STORED IN ASTORE }
          For ii:=1 to n do sbar^[ii]:=-fvec^[ii]*scalef^[ii];
          qrsolv (overch, overfl, maxexp, n, nunit, output, astore, hhpi, rdiag, sbar)
        end
        ELSE
        begin
{       SET UP RIGHT HAND SIDE - MULTIPLY -FVEC BY Q (STORED IN JAC). }
          For ii:=1 to n do wv2^[ii]:=-fvec^[ii]*scalef^[ii];
          avmul (n, n, n, n, jac, wv2, sbar);
          rsolv (overch, overfl, maxexp, n, nunit, output, a, rdiag, sbar)
        end
      end;

{       NORM OF SCALED SBAR IS NEEDED FOR SECOND ACCEPTANCE TEST }

      IF acptcr = 12 THEN
      begin
        For ii:=1 to n do wv2^[ii]:=scalex^[ii]*sbar^[ii];
        twonrm (overfl, maxexp, n, epsmch, sbrnrm, wv2)
      end
    end;

    IF output > 3 THEN
    begin
      Line0;
      Writeln(fp_out,'   *          LINE SEARCH STEP:',k:3, '                                         *');
      Line0;
      Writeln(fp_out,'   *          LAMBDA FOR ATTEMPTED STEP: ',lambda:12:3,'                      *');
      Line0;
      IF (NOT convio) THEN
        Msg('          NEW COMPONENT/FCN VECTORS (XPLUS(I) = XC(I) + LAMBDA*SN(I))')
      ELSE
        Msg('          NEW FUNCTION VECTORS AT MODIFIED POINT');
      Line0;
      For i:=1 to n do
        if (Abs(xplus^[i])<100000) and (Abs(fvec^[i])<100000) then
          Writeln(fp_out,'   *          XPLUS(',i:3,') = ',xplus^[i]:12:3,'     FVEC(',i:3,') = ',fvec^[i]:12:3,'       *')
        else
          Writeln(fp_out,'   *          XPLUS(',i:3,') = ',xplus^[i]:12:-3,'     FVEC(',i:3,') = ',fvec^[i]:12:-3,'       *');
      Line0;
      IF (NOT sclfch) THEN
        if Abs(fcnnew)<100000 then
          Writeln(fp_out,'   *          OBJECTIVE FUNCTION VALUE AT XPLUS: ',fcnnew:12:3,'              *')
        else
          Writeln(fp_out,'   *          OBJECTIVE FUNCTION VALUE AT XPLUS: ',fcnnew:12:-3,'              *')
      ELSE
        if Abs(fcnnew)<100000 then
          Writeln(fp_out,'   *          SCALED OBJECTIVE FUNCTION VALUE AT XPLUS: ',fcnnew:12:3,'       *')
        else
          Writeln(fp_out,'   *          SCALED OBJECTIVE FUNCTION VALUE AT XPLUS: ',fcnnew:12:-3,'       *');
      if Abs(fcnmax+alpha*lambda*delfts)<100000 then
        Writeln(fp_out,'   *          FCNMAX + ALPHA*LAMBDA*DELFTS: ',(fcnmax+alpha*lambda*delfts):12:3,'                   *')
      else
        Writeln(fp_out,'   *          FCNMAX + ALPHA*LAMBDA*DELFTS: ',(fcnmax+alpha*lambda*delfts):12:-3,
        '                   *');

      IF (deuflh) OR (acptcr = 12) THEN
        IF itnum > 0 THEN
          IF (NOT sclxch) THEN
          begin
            Line0;
            Msg('          DEUFLHARD SBAR VECTOR:');
            Line0;
            For i:=1 to n do
              Writeln(fp_out,'   *          SBAR(',i:3,') = ', sbar^[i]:12:3,'                                     *')
          end
          ELSE
          begin
            Line0;
            Msg('          DEUFLHARD SBAR VECTOR IN SCALE DX UNITS:');
            Line0;
            For i:=1 to n do
              Writeln(fp_out,'   *          SBAR(',i:3,') = ', scalex^[i]*sbar^[i]:12:3,'                          *')
          end;

      IF acptcr = 12 THEN
      begin
        Line0;
        IF (NOT sclxch) THEN
          Writeln(fp_out,'   *          VALUE OF SBRNRM AT XPLUS: ..................', sbrnrm:12:3,'     *')
        ELSE
          Writeln(fp_out,'   *          VALUE OF SCALED SBRNRM AT XPLUS: ...........', sbrnrm:12:3,'     *');

        Writeln(fp_out,'   *          NEWMAX:: ...................................', newmax:12:3,'     *');
      end
    end;

{     *** CHECK FOR ACCEPTABLE STEP *** }

    IF fcnnew < fcnmax + alpha*lambda*delfts THEN
    begin
      acpcod:=1;

{     NOTE: STEP WILL BE ACCEPTED REGARDLESS OF NEXT TEST.
            THIS SECTION IS FOR BOOKKEEPING ONLY.   }

      IF acptcr = 12 THEN
        IF sbrnrm < newmax THEN acpcod:=12;

      exit
    end;

    IF (acptcr = 12) AND (sbrnrm < newmax) THEN
    begin
      acpcod:=2;
      exit
    end;

{   FAILURE OF STEP ACCEPTANCE TEST }

    IF (convio) THEN
      lambda:=confac*ratiom*lambda;

{     LAMBDA IS ALREADY SET - SKIP NORMAL PROCEDURE }

      {CYCLE}

    IF lambda = zero THEN
    begin
      IF output > 0 THEN
      begin
        Line0;
        Msg('       LAMBDA IS 0.0:  NO PROGRESS POSSIBLE - CHECK BOUNDS OR START.')
      end;
      abort:=true;
      exit
    end;

    IF (geoms) THEN

{   GEOMETRIC LINE SEARCH }

      lambda:=sigma*lambda

    ELSE
    begin
      IF (frstst) THEN
      begin
        frstst:=false;

{       FIND MINIMUM OF QUADRATIC AT FIRST STEP }

        lamtmp:=-lambda*lambda*delfts/(two*(fcnnew-fcnold-lambda* delfts));
        IF output > 4 THEN
        begin
          Line0;
          Writeln(fp_out,'   *          TEMPORARY LAMBDA FROM QUADRATIC MODEL: ', lamtmp:11:3, '          *')
        end
      end
      ELSE
      begin
{       FIND MINIMUM OF CUBIC AT SUBSEQUENT STEPS }

        factor:=one/(lambda-lampre);
        IF lambda*lambda = zero THEN lambda:=sigma*lambda;

{       NOTE: IF THIS LAMBDA^2 WAS ZERO ANY SUBSEQUENT
              LAMBDA^2 WILL ALSO BE ZERO.

        CYCLE }

        acubic:=factor*((one/lambda*lambda)*(fcnnew-fcnold-lambda*delfts)-((one/lampre*lambda)*(fplpre-fcnold-lampre*delfts)));
        bcubic:=factor*((-lampre/lambda*lambda)*(fcnnew-fcnold-lambda*delfts)+((lambda/lampre*lambda)
        *(fplpre-fcnold-lampre*delfts)));
        IF two*LOG(ABS(bcubic)+eps) > maxexp THEN
        begin
          lamtmp:=sigma*lambda;
          IF output > 4 THEN
          begin
            Line0;
            Msg('           POTENTIAL OVERFLOW IN CALCULATING TRIAL LAMBDA FROM');
            Msg('           CUBIC MODEL - LAMBDA SET TO SIGMA*LAMBDA.')
          end
        end
        ELSE
        begin
          disc:=bcubic*bcubic - three*acubic*delfts;
          IF ABS(acubic) <= epsmch THEN
            lamtmp:=-delfts/(two*bcubic)
          ELSE
            IF disc < zero THEN
              lamtmp:=sigma*lambda
            ELSE
              lamtmp:=(-bcubic + SQRT(disc))/(three*acubic);
          IF output > 4 THEN
          begin
            Line0;
            Writeln(fp_out,'   *            TEMPORARY LAMBDA FROM CUBIC MODEL : .....', lamtmp:11:3,'      *')
          end;
          IF lamtmp > sigma*lambda THEN
          begin
            lamtmp:=sigma*lambda;
            IF output > 4 THEN
            begin
              Line0;
              Msg('            LAMTMP TOO LARGE - REDUCED TO SIGMA*LAMBDA.')
            end
          end
        end {3rd else}
      end {2nd else}
    end; {1st else}
    lampre:=lambda;
    fplpre:=fcnnew;
    IF lamtmp < point1*lambda THEN
    begin
      IF output > 4 THEN
      begin
        Line0;
        Msg('           LAMTMP TOO SMALL - INCREASED TO 0.1*LAMBDA.')
      end;
      lambda:=point1*lambda
    end
    ELSE
    begin
      IF (output > 4) AND (lamtmp <> sigma*lambda) THEN
      begin
        Line0;
        Msg('           LAMTMP WITHIN LIMITS - LAMBDA SET TO LAMTMP.')
      end;
      lambda:=lamtmp
    end
  end; {k loop!}

{   FAILURE IN LINE SEARCH }

  acpcod:=0;

{   IF A QUASI-NEWTON STEP HAS FAILED IN THE LINE SEARCH THEN SET QNFAIL
    TO TRUE ANS RETURN TO SUBROUTINE NNES.  THIS WILL CAUSE THE JACOBIAN
    TO BE RE-EVALUATED EXPLICITLY AND A LINE SEARCH IN THE NEW DIRECTION
    CONDUCTED.   }

  IF (NOT restrt) THEN
  begin
    qnfail:=true;
    exit
  end;

{   FALL THROUGH MAIN LOOP - WARNING GIVEN }

  IF (output > 2) AND (NOT wrnsup) THEN
  begin
    Line0;
    Writeln(fp_out,'   *       WARNING: ', maxlin:3, ' CYCLES COMPLETED IN LINE SEARCH WITHOUT SUCCESS   *')
  end;
  IF (stopcr = 2) OR (stopcr = 3) THEN
  begin
    stopcr:=12;
    Line0;
    Msg('       STOPPING CRITERION RESET FROM 2 TO 12 TO AVOID HANGING.')
  end;

  Dispose(wv2)

End; {deufls}

Procedure fcnevl (var overfl:boolean; maxexp, n, nunit, output:integer;
                  epsmch:double; var fcnnew:double; fvec, scalef:pVec);
{----------------------------------------------------------------------
!    FEB. 23, 1992
!
!    THE OBJECTIVE FUNCTION, FCNNEW, DEFINED BY:
!
!       FCNNEW = 1/2(SCALEF*FVEC^SCALEF*FVEC)
!
!    IS EVALUATED BY THIS SUBROUTINE.
!---------------------------------------------------------------------}
Label return;
Var
  eps: double; k:integer;
  wv1: pVec;
Begin
  New(wv1);
  overfl:=false;
  eps:=one/Power(ten, maxexp);
  For k:=1 to n do wv1^[k]:=fvec^[k]*scalef^[k];
  twonrm (overfl, maxexp, n, epsmch, fcnnew, wv1);

{ IF AN OVERFLOW WOULD OCCUR SUBSTITUTE A LARGE VALUE FOR FCNNEW }

  IF (overfl) OR (two*LOG(fcnnew+eps) > maxexp) THEN
  begin
    overfl:=true;
    fcnnew:=Power(ten, maxexp);
    IF (output > 2) AND (NOT wrnsup) THEN
    begin
      Line0;
      Writeln('   *    WARNING: TO AVOID OVERFLOW, OBJECTIVE FUNCTION SET TO: ',fcnnew:11:3,'*')
    end;
    goto RETURN
  end;
  fcnnew:=fcnnew*fcnnew/two;
return: Dispose(wv1)
End; {fcnevl}

Procedure fordif (var overfl:boolean; j, n:integer; deltaj:double; fvec, fvecj1:pVec;
                 jacfdm:pMat; xc:pVec);
{----------------------------------------------
!       FEB. 6, 1991
!---------------------------------------------}
var k:integer;
Begin
  fcn1 (overfl, n, fvecj1, xc);
  IF (NOT overfl) THEN
    For k:=1 to n do
      jacfdm^[k,j]:=(fvecj1^[k] - fvec^[k])/deltaj;
End;

Procedure innerp (overch:boolean; var overfl:boolean; maxexp, ldima, ldimb, n, nunit,
                  output:integer; var dtpro:double; a, b:pVec);
{--------------------------------------------------------------------
!    FEB. 14, 1991
!
!    THIS SUBROUTINE FINDS THE INNER PRODUCT OF TWO VECTORS, A AND B.
!    IF OVERCH IS FALSE, UNROLLED LOOPS ARE USED.
!
!    LDIMA IS THE DIMENSION OF A
!    LDIMB IS THE DIMENSION OF B
!    N IS THE DEPTH INTO A AND B THE INNER PRODUCT IS DESIRED.
!    ^[USUALLY LDIMA:=LDIMB:=N)
!-------------------------------------------------------------------}
Label return;
Var
    eps:double;
    i, k, kk, ng4, ng4r, ng8, ng8r, ng16, ng16r, ng32, ng32r:integer;
Begin
  eps:=one/Power(ten, maxexp);
  overfl:=false;

  dtpro:=zero;
  IF (overch) THEN
    For i:=1 to n do
    begin
      IF LOG(ABS(a^[i])+eps) + LOG(ABS(b^[i])+eps) > maxexp THEN
      begin
        overfl:=true;
        dtpro:=SIGN(Power(ten,maxexp),a^[i])*SIGN(one,b^[i]);
        IF (output > 2) AND (NOT wrnsup) THEN
        begin
          Line0;
          Writeln(fp_out,'   *    WARNING: TO AVOID OVERFLOW, INNER PRODUCT SET TO ', dtpro:12:3, '     *')
        end;
        goto RETURN
      end;
      dtpro:=dtpro + a^[i]*b^[i]
    end
  ELSE
  begin

{   SET NUMBER OF GROUPS OF EACH SIZE }

    ng32:=n Div 32;
    ng32r:=n - 32*ng32;
    ng16:=ng32r Div 16;
    ng16r:=ng32r - 16*ng16;
    ng8:=ng16r Div 8;
    ng8r:=ng16r - 8*ng8;
    ng4:=ng8r Div 4;
    ng4r:=ng8r - 4*ng4;

{   FIND INNER PRODUCT }

    k:=0;
    IF ng32 > 0 THEN
      For kk:=1 to ng32 do
      begin
        k:=k+32;
        dtpro:=dtpro + a^[k-31]*b^[k-31] + a^[k-30]*b^[k-30] + a^[k-29]*b^[k-29] +
            a^[k-28]*b^[k-28] + a^[k-27]*b^[k-27] + a^[k-26]*b^[k-26] +
            a^[k-25]*b^[k-25] + a^[k-24]*b^[k-24];
        dtpro:=dtpro + a^[k-23]*b^[k-23] + a^[k-22]*b^[k-22] + a^[k-21]*b^[k-21] +
            a^[k-20]*b^[k-20] + a^[k-19]*b^[k-19] + a^[k-18]*b^[k-18] +
            a^[k-17]*b^[k-17] + a^[k-16]*b^[k-16];
        dtpro:=dtpro + a^[k-15]*b^[k-15] + a^[k-14]*b^[k-14] + a^[k-13]*b^[k-13] +
            a^[k-12]*b^[k-12] + a^[k-11]*b^[k-11] + a^[k-10]*b^[k-10] +
            a^[k-9]*b^[k-9] + a^[k-8]*b^[k-8];
        dtpro:=dtpro + a^[k-7]*b^[k-7] + a^[k-6]*b^[k-6] + a^[k-5]*b^[k-5] + a^[k-4]
          *b^[k-4] + a^[k-3]*b^[k-3] + a^[k-2]*b^[k-2] + a^[k-1]*b^[k-1] + a^[k]*b^[k]
      end;

    IF ng16 > 0 THEN
      For kk:=1 to ng16 do
      begin
        k:=k+16;
        dtpro:=dtpro + a^[k-15]*b^[k-15] + a^[k-14]*b^[k-14] + a^[k-13]*b^[k-13] +
            a^[k-12]*b^[k-12] + a^[k-11]*b^[k-11] + a^[k-10]*b^[k-10] +
            a^[k-9]*b^[k-9] + a^[k-8]*b^[k-8];
        dtpro:=dtpro + a^[k-7]*b^[k-7] + a^[k-6]*b^[k-6] + a^[k-5]*b^[k-5] +
            a^[k-4]*b^[k-4] + a^[k-3]*b^[k-3] + a^[k-2]*b^[k-2] + a^[k-1]*b^[k-1] +
            a^[k]*b^[k]
      end;

    IF ng8 > 0 THEN
      For kk:=1 to ng8 do
      begin
        k:=k+8;
        dtpro:=dtpro + a^[k-7]*b^[k-7] + a^[k-6]*b^[k-6] + a^[k-5]*b^[k-5] +
            a^[k-4]*b^[k-4] + a^[k-3]*b^[k-3] + a^[k-2]*b^[k-2] + a^[k-1]*b^[k-1] +
            a^[k]*b^[k]
      end;

    IF ng4 > 0 THEN
      For kk:=1 to ng4 do
      begin
        k:=k+4;
        dtpro:=dtpro + a^[k-3]*b^[k-3] + a^[k-2]*b^[k-2] + a^[k-1]*b^[k-1] + a^[k]*b^[k]
      end;

    IF ng4r > 0 THEN
      For kk:=1 to ng4r do
      begin
        k:=k+1;
        dtpro:=dtpro + a^[k]*b^[k]
      end

  end;

Return: End; {innerp}


Procedure jacrot (var overfl:boolean; i, maxexp, n:integer; arot, brot, epsmch:double;
                   a, jac:pMat);
{----------------------------------------------------------
!    FEB. 11, 1991
!
!    JACOBI (OR GIVENS) ROTATION.
!---------------------------------------------------------}
Var
    c,denom,s,w,y: double;
    hold: array[1..2] of double;
    j, ldhold:integer;
Begin
  IF arot = zero THEN
  begin
    c:=zero;
    s:=-SIGN(one,brot)
  end
  ELSE
  begin
    hold[1]:=arot;
    hold[2]:=brot;
    ldhold:=2;
    twonrm (overfl, maxexp, ldhold, epsmch, denom, @hold);
    c:=arot/denom;
    s:=-brot/denom
  end;
  For j:=i to n do
  begin
    y:=a^[i,j];
    w:=a^[i+1,j];
    a^[i,j]:=c*y - s*w;
    a^[i+1,j]:=s*y + c*w
  end;
  For j:=1 to n do
  begin
    y:=jac^[i,j];
    w:=jac^[i+1,j];
    jac^[i,j]:=c*y - s*w;
    jac^[i+1,j]:=s*y + c*w
  end
End; {jacrot}

Procedure lsolv (overch:boolean; var overfl:boolean; maxexp, n, nunit,
                 output:integer; l:pMat; b, rhs:pVec);
{-----------------------------------------------------------------------------
!    FEB. 14, 1991
!
!    THIS SUBROUTINE SOLVES:
!
!           LB=RHS
!
!           WHERE    L     IS TAKEN FROM THE CHOLESKY DECOMPOSITION
!                    RHS   IS A GIVEN RIGHT HAND SIDE WHICH IS NOT OVERWRITTEN
!                    B     IS THE SOLUTION VECTOR
!
!    FRSTER IS USED FOR OUTPUT PURPOSES ONLY.
!----------------------------------------------------------------------------}
Label 30;
Var
    i, j, jstar:integer;
    eps, maxlog, SUM, tmplog:double;
    frster:boolean;
Begin

  frster:=true;
  overfl:=false;
  eps:=one/Power(ten,maxexp);

  IF (overch) THEN
    IF LOG(ABS(rhs^[1])+eps) - LOG(ABS(l^[1,1])+eps) > maxexp THEN
    begin
      overfl:=true;
      b^[1]:=SIGN(Power(ten,maxexp), rhs^[1])*SIGN(one,l^[1,1]);
      IF output > 3 THEN
      begin
        Line0;
        Writeln(fp_out,'   *    WARNING: COMPONENT 1 SET TO ',b^[1]:12:3,'                          *')
      end;
      GOTO 30
    end;

  b^[1]:=rhs^[1]/l^[1,1];

30: For i:=2 to n do
  begin
    IF (overch) THEN
    begin
{     CHECK TO FIND IF ANY TERMS IN THE EVALUATION WOULD OVERFLOW }
      maxlog:=LOG(ABS(rhs^[i])+eps) - LOG(ABS(l^[i,i])+eps);
      jstar:=0;
      For j:=1 to i-1 do
      begin
        tmplog:=LOG(ABS(l^[i,j])+eps) + LOG(ABS(b^[j])+eps) - LOG(ABS(l^[i,i])+eps);
        IF tmplog > maxlog THEN
        begin
          jstar:=j;
          maxlog:=tmplog
        end
      end;

{     IF AN OVERFLOW WOULD OCCUR ASSIGN A VALUE FOR THE
      TERM WITH CORRECT SIGN. }

      IF maxlog > maxexp THEN
      begin
        overfl:=true;
        IF jstar = 0 THEN
          b^[i]:=SIGN(Power(ten,maxexp), rhs^[i])*SIGN(one,l^[i,i])
        ELSE
          b^[i]:=-SIGN(Power(ten,maxexp), l^[i,jstar])*SIGN(one,b^[jstar])*SIGN(one,l^[i,i]);
        IF (frster) THEN
        begin
          frster:=false;
          Line0
        end;
        IF output > 3 then
          Writeln(fp_out,'   *    WARNING: COMPONENT ',i:3,' SET TO ',b^[1]:12:3,'                        *')
      end
    end;

{   SUM FOR EACH TERM, ORDERING OPERATIONS TO MINIMIZE
    POSSIBILITY OF OVERFLOW.  }

    sum:=zero;
    For j:=1 to i-1 do
      sum:=sum + (MIN(ABS(l^[i,j]),ABS(b^[j]))/l^[i,i])*(MAX(ABS(l^[i,j]),ABS(b^[j])))*SIGN(one,l^[i,j])*SIGN(one,b^[j]);
    b^[i]:=rhs^[i]/l^[i,i] - sum
  end {i loop}

End; {lsolv}

Procedure ltsolv (overch:boolean; var overfl:boolean; maxexp, n, nunit,
                  output:integer; l:pMat; y, b:pVec);
{--------------------------------------------------------------------------
!    FEB. 14, 1991
!
!    THIS SUBROUTINE SOLVES:
!
!           LY=B
!
!           WHERE    L  IS TAKEN FROM THE CHOLESKY DECOMPOSITION
!                    B  IS A GIVEN RIGHT HAND SIDE WHICH IS NOT OVERWRITTEN
!                    Y  IS THE SOLUTION VECTOR
!
!    FRSTER IS USED FOR OUTPUT PURPOSES ONLY.
!-------------------------------------------------------------------------}
Label 30;
Var
    eps, maxlog, SUM, tmplog:double;
    frster:boolean;
    i, j, jstar:integer;
Begin
  frster:=true;
  overfl:=false;
  eps:=one/Power(ten,maxexp);

  IF (overch) THEN
    IF LOG(ABS(b^[n])+eps) - LOG(ABS(l^[n,n])+eps) > maxexp THEN
    begin
      overfl:=true;
      y^[n]:=SIGN(Power(ten,maxexp),b^[n])*SIGN(one,l^[n,n]);
      IF output > 3 THEN
      begin
        frster:=false;
        Line0;
        Writeln(fp_out,'   *    WARNING: COMPONENT ',n:3,' SET TO ',y^[n]:11:3,'                         *')
      end;
      GOTO 30
    end;

  y^[n]:=b^[n]/l^[n,n];

30: For i:=n-1 Downto 1 do
  begin
    IF (overch) THEN
    begin

{     CHECK TO FIND IF ANY TERMS IN THE EVALUATION WOULD OVERFLOW }

      maxlog:=LOG(ABS(b^[i])+eps) - LOG(ABS(l^[i,i])+eps);
      jstar:=0;
      For j:=i+1 to n do
      begin
        tmplog:=LOG(ABS(l^[j,i])+eps) + LOG(ABS(y^[j])+eps) - LOG(ABS(l^[i,i])+eps);
        IF tmplog > maxlog THEN
        begin
          jstar:=j;
          maxlog:=tmplog
        end
      end;

{     IF AN OVERFLOW WOULD OCCUR ASSIGN A VALUE FOR THE
      TERM WITH CORRECT SIGN. }

      IF maxlog > maxexp THEN
      begin
        overfl:=true;
        IF jstar = 0 THEN
          y^[i]:=SIGN(Power(ten,maxexp),b^[i])*SIGN(one,l^[i,i])
        ELSE
          y^[i]:=-SIGN(Power(ten,maxexp),l^[jstar,i])*SIGN(one,y^[jstar])*SIGN(one,l^[i,i]);
        IF (frster) THEN
        begin
          frster:=false;
          Line0
        end;
        IF output > 3 then
          Writeln(fp_out,'   *    WARNING: COMPONENT ',i:3,' SET TO ',y^[i]:11:3,'                         *')
      end
    end;

{   SUM FOR EACH TERM ORDERING OPERATIONS TO MINIMIZE
    POSSIBILITY OF OVERFLOW. }

    sum:=zero;
    For j:=i+1 to n do
      sum:=sum + (MIN(ABS(l^[j,i]),ABS(y^[j]))/l^[i,i])*(MAX(ABS(l^[j,i]),ABS(y^[j])))*SIGN(one,l^[j,i])*SIGN(one,y^[j]);
    y^[i]:=b^[i]/l^[i,i] - sum
  end

End; {ltsolv}

Procedure matcop (nradec, nraact, ncadec, ncaact, nrbdec, ncbdec:integer; amat, bmat:pMat);
{------------------------------------------------------------------
!    SEPT. 15, 1991
!
!    COPY A CONTINGUOUS RECTANGULAR PORTION OF ONE MATRIX
!    INTO ANOTHER ^[ELEMENT ^[1,1] MUST BE INCLUDED].
!
!    NRADEC IS 1ST DIMENSION OF AMAT, NRAACT IS LIMIT OF 1ST INDEX
!    NCADEC IS 2ND DIMENSION OF AMAT, NCAACT IS LIMIT OF 2ND INDEX
!    NRBDEC IS 1ST DIMENSION OF BMAT
!    NCBDEC IS 2ND DIMENSION OF BMAT
!-----------------------------------------------------------------}
Var
    j, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, ncc32r:integer;
Begin

{   FIND NUMBER OF GROUPS OF SIZE 32, 16...}

  ncc32:=nraact Div 32;
  ncc32r:=nraact - 32*ncc32;
  ncc16:=ncc32r Div 16;
  ncc16r:=ncc32r - 16*ncc16;
  ncc8:=ncc16r Div 8;
  ncc8r:=ncc16r - 8*ncc8;
  ncc4:=ncc8r Div 4;
  ncc4r:=ncc8r - 4*ncc4;

  For j:=1 to ncaact do
  begin

{   COPY ENTRIES INTO MATRIX B BY COLUMN }
    k:=0;
    IF ncc32 > 0 THEN
      For kk:=1 to ncc32 do
      begin
        k:=k+32;
        bmat^[k-31,j]:=amat^[k-31,j];
        bmat^[k-30,j]:=amat^[k-30,j];
        bmat^[k-29,j]:=amat^[k-29,j];
        bmat^[k-28,j]:=amat^[k-28,j];
        bmat^[k-27,j]:=amat^[k-27,j];
        bmat^[k-26,j]:=amat^[k-26,j];
        bmat^[k-25,j]:=amat^[k-25,j];
        bmat^[k-24,j]:=amat^[k-24,j];
        bmat^[k-23,j]:=amat^[k-23,j];
        bmat^[k-22,j]:=amat^[k-22,j];
        bmat^[k-21,j]:=amat^[k-21,j];
        bmat^[k-20,j]:=amat^[k-20,j];
        bmat^[k-19,j]:=amat^[k-19,j];
        bmat^[k-18,j]:=amat^[k-18,j];
        bmat^[k-17,j]:=amat^[k-17,j];
        bmat^[k-16,j]:=amat^[k-16,j];
        bmat^[k-15,j]:=amat^[k-15,j];
        bmat^[k-14,j]:=amat^[k-14,j];
        bmat^[k-13,j]:=amat^[k-13,j];
        bmat^[k-12,j]:=amat^[k-12,j];
        bmat^[k-11,j]:=amat^[k-11,j];
        bmat^[k-10,j]:=amat^[k-10,j];
        bmat^[k-9,j]:=amat^[k-9,j];
        bmat^[k-8,j]:=amat^[k-8,j];
        bmat^[k-7,j]:=amat^[k-7,j];
        bmat^[k-6,j]:=amat^[k-6,j];
        bmat^[k-5,j]:=amat^[k-5,j];
        bmat^[k-4,j]:=amat^[k-4,j];
        bmat^[k-3,j]:=amat^[k-3,j];
        bmat^[k-2,j]:=amat^[k-2,j];
        bmat^[k-1,j]:=amat^[k-1,j];
        bmat^[k,j]:=amat^[k,j]
      end;

    IF ncc16 > 0 THEN
      For kk:=1 to ncc16 do
      begin
        k:=k+16;
        bmat^[k-15,j]:=amat^[k-15,j];
        bmat^[k-14,j]:=amat^[k-14,j];
        bmat^[k-13,j]:=amat^[k-13,j];
        bmat^[k-12,j]:=amat^[k-12,j];
        bmat^[k-11,j]:=amat^[k-11,j];
        bmat^[k-10,j]:=amat^[k-10,j];
        bmat^[k-9,j]:=amat^[k-9,j];
        bmat^[k-8,j]:=amat^[k-8,j];
        bmat^[k-7,j]:=amat^[k-7,j];
        bmat^[k-6,j]:=amat^[k-6,j];
        bmat^[k-5,j]:=amat^[k-5,j];
        bmat^[k-4,j]:=amat^[k-4,j];
        bmat^[k-3,j]:=amat^[k-3,j];
        bmat^[k-2,j]:=amat^[k-2,j];
        bmat^[k-1,j]:=amat^[k-1,j];
        bmat^[k,j]:=amat^[k,j]
      end;

    IF ncc8 > 0 THEN
      For kk:=1 to ncc8 do
      begin
        k:=k+8;
        bmat^[k-7,j]:=amat^[k-7,j];
        bmat^[k-6,j]:=amat^[k-6,j];
        bmat^[k-5,j]:=amat^[k-5,j];
        bmat^[k-4,j]:=amat^[k-4,j];
        bmat^[k-3,j]:=amat^[k-3,j];
        bmat^[k-2,j]:=amat^[k-2,j];
        bmat^[k-1,j]:=amat^[k-1,j];
        bmat^[k,j]:=amat^[k,j]
      end;

    IF ncc4 > 0 THEN
      For kk:=1 to ncc4 do
      begin
        k:=k+4;
        bmat^[k-3,j]:=amat^[k-3,j];
        bmat^[k-2,j]:=amat^[k-2,j];
        bmat^[k-1,j]:=amat^[k-1,j];
        bmat^[k,j]:=amat^[k,j]
      end;

    IF ncc4r > 0 THEN
      For kk:=1 to ncc4r do
      begin
        k:=k+1;
        bmat^[k,j]:=amat^[k,j]
      end
  end
End; {matcop}

Procedure qrsolv(overch:boolean; var overfl:boolean; maxexp, n, nunit, output:integer;
                 a:pMat; hhpi, rdiag, b: pVec);
{---------------------------------------------------------------------------
!    FEB. 2, 1991
!
!    THIS SUBROUTINE SOLVES
!
!         ^[QR]X:=B
!
!         WHERE  Q AND R ARE OBTAINED FROM THE QR DECOMPOSITION
!                B IS A GIVEN RIGHT HAND SIDE WHICH IS OVERWRITTEN
!
!                R IS CONTAINED IN THE STRICT UPPER TRIANGLE OF
!                  MATRIX A AND THE VECTOR RDIAG
!                Q IS "CONTAINED" IN THE LOWER TRIANGLE OF MATRIX A
!
!    FRSTOV  INDICATES FIRST OVERFLOW - USED ONLY TO SET BORDER FOR OUTPUT
!--------------------------------------------------------------------------}
Var
    frstov:boolean;
    eps, tau:double;
    i, j:integer;
Begin
  eps:=one/Power(ten, maxexp);
  frstov:=true;
  overfl:=false;

{ MULTIPLY RIGHT HAND SIDE BY Q THEN SOLVE USING R STORED IN MATRIX A }

  For j:=1 to n-1 do
  begin
    tau:=zero;
    For i:=j to n do
    begin
      IF (overch) THEN
        IF LOG(ABS(a^[i,j])+eps) + LOG(ABS(b^[i])+eps) - LOG(hhpi^[j]+eps) > maxexp THEN
        begin
          overfl:=true;
          tau:=SIGN(Power(ten,maxexp), a^[i,j])*SIGN(one,b^[j])
        end;
      tau:=tau + a^[i,j]*b^[i]/hhpi^[j]
    end;

    For i:=j to n do
    begin
      IF (overch) THEN
        IF LOG(ABS(tau)+eps) + LOG(ABS(a^[i,j])+eps) > maxexp THEN
        begin
          overfl:=true;
          b^[i]:=-SIGN(Power(ten,maxexp),tau)*SIGN(one,a^[i,j]);
          IF (output > 2) AND (NOT wrnsup) THEN
          begin
            IF (frstov) THEN Line0;
            Writeln(fp_out,'   *    WARNING: COMPONENT ',i:3,' SET TO ',b^[i]:11:3,' IN QRSOLV BEFORE RSOLV  *')
          end
        end;
      b^[i]:=b^[i] - tau*a^[i,j]
    end
  end;

  rsolv (overch, overfl, maxexp, n, nunit, output, a, rdiag, b);

End; {qrsolv}


Procedure qrupda (var overfl:boolean; maxexp, n:integer; epsmch:double;
                  a, jac:pMat; u, v:pVec);
{-----------------------------------------------------------------
!    FEB. 12, 1991
!
!    UPDATE QR DECOMPOSITION USING A SERIES OF GIVENS ROTATIONS.
!----------------------------------------------------------------}
Label return;
Var
    eucnrm:double;
    hold: array[1..2] of double; 
    i, k, l, ldhold:integer;
Begin

{   REPLACE SUBDIAGONAL WITH ZEROS SO THAT WHEN R IS MULTIPLIED BY GIVENS
    (JACOBI) ROTATIONS THE SUBDIAGONAL ELEMENTS DO NOT AFFECT THE OUTCOME }

  For i:=2 to n do a^[i,i-1]:=zero;

{  FIND LARGEST K FOR WHICH U^[K] DOES NOT EQUAL ZERO }

  k:=n;
  For l:=1 to n do
    IF u^[k] = zero THEN
      IF k > 1 THEN
        k:=k-1
      ELSE
        goto return
    ELSE
      goto return;

{    MULTIPLY UV^ BY A SERIES OF ROTATIONS SO THAT ALL BUT THE
     TOP ROW IS MADE ZERO ^[THEORETICALLY THIS IS WHAT HAPPENS
     ALTHOUGH THIS MATRIX ISN'T ACTUALLY FORMED]. }

  For i:=k-1 Downto 1 do
  begin
    jacrot (overfl, i, maxexp, n, u^[i], u^[i+1], epsmch, a, jac);
    IF u^[i] = zero THEN
{    THIS STEP JUST AVOIDS ADDING ZERO }
      u^[i]:=ABS(u^[i+1])
    ELSE
    begin
      hold[1]:=u^[i];
      hold[2]:=u^[i+1];
      ldhold:=2;
      twonrm (overfl, maxexp, ldhold, epsmch, eucnrm, @hold);
      u^[i]:=eucnrm
    end
  end;

{   ADD THE TOP ROW TO THE TOP ROW OF A - THIS FORMS THE
    UPPER HESSENBERG MATRIX. }
  For i:= 1 to n do  a^[1,i]:=a^[1,i] + u^[1]*v^[i];

{   FORM THE UPPER TRIANGULAR R MATRIX BY A SERIES OF ROTATIONS
    TO ZERO OUT THE SUBDIAGONALS. }

  For  i:=1 to k-1 do
    jacrot (overfl, i, maxexp, n, a^[i,i], a^[i+1,i], epsmch, a, jac);

return: End; {qrupda}


Procedure rsolv (overch:boolean; var overfl:boolean; maxexp, n, nunit,
                 output:integer; a:pMat; rdiag, b:pVec);
{--------------------------------------------------------------------------
!    FEB. 14, 1991
!
!    THIS SUBROUTINE SOLVES, BY BACKWARDS SUBSTITUTION,
!
!           RX:=B
!
!           WHERE    R IS TAKEN FROM THE QR DECOMPOSITION AND IS STORED IN
!                      THE STRICT UPPER TRIANGLE OF MATRIX A AND THE VECTOR,
!                      RDIAG
!                    B IS A GIVEN RIGHT HAND SIDE WHICH IS OVERWRITTEN
!
!    FRSTOV  INDICATES FIRST OVERFLOW - USED ONLY TO SET BORDERS FOR OUTPUT
!-------------------------------------------------------------------------}
Label 30;
Var
    eps, maxlog, SUM, tmplog:double;
    frstov:boolean;
    i, j, jstar:integer;
Begin
  frstov:=true;
  overfl:=false;
  eps:=one/Power(ten, maxexp);

  IF (overch) THEN
    IF LOG(ABS(b^[n])+eps) - LOG(ABS(rdiag^[n])+eps) > maxexp THEN
    begin
      overfl:=true;
      b^[n]:=SIGN(Power(ten,maxexp),b^[n])*SIGN(one,rdiag^[n]);
      IF (output > 2) AND (NOT wrnsup) THEN
      begin
        frstov:=false;
        Line0;
        Writeln(fp_out,'   *    WARNING: COMPONENT ',n:3,' SET TO ',b^[n]:12:3,'                         *')
      end;
      GOTO 30
    end;

  b^[n]:=b^[n]/rdiag^[n];

30: For i:=n-1 Downto 1 do
  begin
    IF (overch) THEN
    begin
{   CHECK TO FIND IF ANY TERMS IN THE EVALUATION WOULD OVERFLOW }

      maxlog:=LOG(ABS(b^[i])+eps) - LOG(ABS(rdiag^[i])+eps);
      jstar:=0;
      For j:=i+1 to n do
      begin
        tmplog:=LOG(ABS(a^[i,j])+eps) + LOG(ABS(b^[j])+eps) - LOG(ABS(rdiag^[i])+eps);
        IF tmplog > maxlog THEN
        begin
          jstar:=j;
          maxlog:=tmplog
        end
      end;

{     IF AN OVERFLOW WOULD OCCUR ASSIGN A VALUE FOR THE
      TERM WITH CORRECT SIGN.  }

      IF maxlog > maxexp THEN
      begin
        overfl:=true;
        IF jstar = 0 THEN
          b^[i]:=SIGN(Power(ten,maxexp),b^[i])*SIGN(one,rdiag^[i])
        ELSE
          b^[i]:=-SIGN(Power(ten,maxexp),a^[i,jstar])*SIGN(one,b^[jstar])*SIGN(one,rdiag^[i]);
        IF (frstov) THEN
        begin
          frstov:=false;
          Line0;
        end;
        IF (output > 2) AND (NOT wrnsup) Then
          Writeln(fp_out,'   *    WARNING: COMPONENT ',i:3,' SET TO ',b^[i]:12:3,'                         *')
      end
    end;

{   SUM FOR EACH TERM ORDERING OPERATIONS TO MINIMIZE
    POSSIBILITY OF OVERFLOW.  }

    sum:=zero;
    For j:=i+1 to n do
      sum:=sum + (MIN(ABS(a^[i,j]),ABS(b^[j]))/rdiag^[i])*(MAX(ABS(a^[i,j]),ABS(b^[j])))*SIGN(one,a^[i,j])*SIGN(one,b^[j]);
    b^[i]:=b^[i]/rdiag^[i]-sum
  end

End; {rsolv}

Procedure setup (var absnew, cauchy, deuflh, geoms, linesr, newton, overch:boolean;
                 var acptcr, itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns,
                 maxqns, minqns:integer; n:integer; var narmij, niejev, njacch,
                 output, qnupdm, stopcr, supprs, trupdm: integer; var alpha, confac,
                 delta, delfac, epsmch, etafac, fdtolj, ftol, lam0, mstpf, nsttol,
                 omega, ratiof, sigma, stptol: double; boundl, boundu, scalef,
                 scalex: pVec; var help:String);
{------------------------------------------------------------------------
!    DEC. 7, 1991
!
!    SUBROUTINE SETUP ASSIGNS DEFAULT VALUES TO ALL REQUISITE PARAMETERS.
!-----------------------------------------------------------------------}
Var
    temp, xmax: double;
    i, ibeta, it, maxebb, minebb: integer;
Begin

{ LOGICAL VALUES }

  absnew:=false;
  bypass:=false;
  cauchy:=false;
  deuflh:=true;
  geoms:= true;
  linesr:=true;
  matsup:=false;
  newton:=false;
  overch:=false;
  wrnsup:=false;

{ INTEGER VALUES }

  acptcr:=12;
  itsclf:=0;
  itsclx:=0;
  jactyp:=1;
  jupdm:=0;
  maxit:=100;  {max. number of iterations}
  maxns:=50;
  maxqns:=10;
  minqns:=7;
  narmij:=1;
  nfetot:=0;
  niejev:=1;
  njacch:=1;
  output:=2;   {set to 4 to have more outputs}
  qnupdm:=1;
  stopcr:=12;
  supprs:=0;
  trupdm:=0;

{ REAL VALUES }

  alpha:=1E-04;
  confac:=0.95;
  delta:=-1.0;
  delfac:=2.0;
  etafac:=0.2;
  lam0:=1.0;
  mstpf:=1E3;
  omega:=0.1;
  ratiof:=0.7;
  sigma:=0.5;

{ CHARACTER VARIABLE }

  help:='NONE';

{ NOTE: NOTATIONAL CHANGES IN CALLING PROGRAM FROM MACHAR
        1)  EPSMCH DENOTES MACHINE EPSILON
        2)  MINEBB DENOTES MINIMUM EXPONENT BASE BETA
        3)  MAXEBB DENOTES MAXIMUM EXPONENT BASE BETA  }

{Values extracted from Fortan 90 version}
  it     := 53;
  ibeta  := 2;
  minebb := -1021;
  maxebb := 1024;
  epsmch := 2.220446049250313E-16;
  xmax   := 1.7976931348E+300;
  maxexp := 308;

{ VALUES FOR TWO-NORM CALCULATIONS }

  smallb:=Power(one*ibeta,(minebb+1) div 2);
  bigb:=Power(one*ibeta,(maxebb-it+1) div 2);
  smalls:=Power(one*ibeta,(minebb-1) div 2);
  bigs:=Power(one*ibeta,(maxebb+it-1) div 2);
  bigr:=xmax;

{ SET STOPPING CRITERIA PARAMETERS }

  fdtolj:=1E-06;
  ftol:=Power1(epsmch,0.333);
  nsttol:=ftol*ftol;
  stptol:=nsttol;

{ VECTOR VALUES }

  temp:=-Power(ten, maxexp);
  For i:=1 to n do
  begin
    boundl^[i]:=temp;
    boundu^[i]:=-temp;
    scalef^[i]:=one;
    scalex^[i]:=one
  end

End; {setup}


Procedure twonrm (var overfl:boolean; maxexp, n:integer; epsmch:double;
                  var eucnrm:double; v:pVec);
{----------------------------------------------------------------------
!    FEB. 23 ,1992
!
!    THIS SUBROUTINE EVALUATES THE EUCLIDEAN NORM OF A VECTOR, V.
!    IT FOLLOWS THE ALGORITHM OF J.L. BLUE IN ACM TOMS V4 15 (1978)
!    BUT USES SLIGHTLY DIFFERENT CUTS.   THE CONSTANTS IN COMMON BLOCK
!    NNES_5 ARE CALCULATED IN THE SUBROUTINE MACHAR OR ARE PROVIDED
!    BY THE USER IN THE DRIVER.
!---------------------------------------------------------------------}
Label return;
Var
    abig, absvi, amed, asmall, sqrtep, ymax, ymin:double;
    i:integer;
Begin
  overfl:=false;
  sqrtep:=SQRT(epsmch);

  asmall:=zero;
  amed:=zero;
  abig:=zero;

  For i:=1 to n do
  begin
{   ACCUMULATE SUMS OF SQUARES IN ONE OF THREE ACCULULATORS,
    ABIG, AMED AND ASMALL, DEPENDING ON THEIR SIZES. }

    absvi:=ABS(v^[i]);

{   THIS COMPARISON RESTRICTS THE MAXIMUM VALUE OF AMED TO BE
    B/N => CANNOT SUM SO THAT AMED OVERFLOWS. }

    IF absvi > bigb/(one*n) THEN
{   THIS DIVISOR OF 10N RESTRICTS ABIG FROM (PATHALOGICALLY)
    OVERFLOWING FROM SUMMATION. }
      abig:=abig + ((v^[i]/Sqr(ten*n*bigs)))
    ELSE IF absvi < smalls THEN
      asmall:=asmall + Sqr(v^[i]/smalls)
    ELSE
      amed:=amed + v^[i]*v^[i]
  end;

  IF abig > zero THEN
  begin

{   IF OVERFLOW WOULD OCCUR ASSIGN BIGR AS NORM AND SIGNAL TO
    CALLING SUBROUTINE VIA OVERFL. }

    IF LOG(abig)/two + LOG(bigs) + one + LOG(one*n) > maxexp THEN
    begin
      eucnrm:=bigr;
      overfl:=true;
      goto RETURN
    end;

{   IF AMED IS POSITIVE IT COULD CONTRIBUTE TO THE NORM -
    DETERMINATION IS DELAYED UNTIL LATER TO SAVE REPEATING CODE. }

    IF amed > zero THEN
    begin
      ymin:=MIN(SQRT(amed),ten*n*bigs*SQRT(abig));
      ymax:=MAX(SQRT(amed),ten*n*bigs*SQRT(abig))
    end
    ELSE
    begin
{   AMED DOESN'T CONTRIBUTE AND ASMALL WON'T MATTER IF
    ABIG IS NONZERO - FIND NORM USING ABIG AND RETURN. }
      eucnrm:=ten*n*bigs*SQRT(abig);
      goto RETURN
    end
  end
  ELSE IF (asmall > zero) THEN
    IF (amed > zero) THEN
    begin
      ymin:=MIN(SQRT(amed),smalls*SQRT(asmall));
      ymax:=MAX(SQRT(amed),smalls*SQRT(asmall))
    end
    ELSE
    begin
{   ABIG AND AMED ARE ZERO SO USE ASMALL ONLY. }
      eucnrm:=smalls*SQRT(asmall);
      goto RETURN
    end
  ELSE
  begin
    eucnrm:=SQRT(amed);
    goto RETURN
  end;

  IF ymin < sqrtep*ymax THEN
{ SMALLER PORTION DOES NOT CONTRIBUTE TO NORM. }
    eucnrm:=ymax
  ELSE
{ SMALLER PORTION CONTRIBUTES TO NORM. }
    eucnrm:=ymax*SQRT((one+ymin/ymax)*(one+ymin/ymax));

Return: End; {twonrm}

Procedure update (mnew:integer; var mold:integer; n:integer; var trmcod:
                  integer; var fcnnew, fcnold:double;  fvec, fvecc, xc, xplus:pVec);
{------------------------------------------------------------------------
!    FEB. 9, 1991
!
!    THIS SUBROUTINE RESETS CURRENT ESTIMATES OF SOLUTION AND UPDATES THE
!    OBJECTIVE FUNCTION VALUE, M (USED TO SET HOW MANY PREVIOUS VALUES TO
!    LOOK AT IN THE NON- MONOTONIC COMPARISONS) AND THE TERMINATION CODE,
!    TRMCOD.
!-----------------------------------------------------------------------}
Var
    i: integer;
Begin
  fcnold:=fcnnew;
  mold:=mnew;
  trmcod:=0;
  For  i:=1 to n do
  begin
    fvecc^[i]:=fvec^[i];
    xc^[i]:=xplus^[i]
  end
End; { update}

Procedure utbmul (ncadec, ncaact, ncbdec, ncbact, ncdec, ncact:integer;
                  amat, bmat, cmat:pMat);
{---------------------------------------------------------------------------
!    FEB. 8, 1991
!
!    MATRIX MULTIPLICATION:   A^B:=C    WHERE A IS UPPER TRIANGULAR
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4.
!
!    NCADEC IS 2ND DIM. OF AMAT; NCAACT IS ACTUAL LIMIT FOR 2ND INDEX
!    NCBDEC IS 2ND DIM. OF BMAT; NCBACT IS ACTUAL LIMIT FOR 2ND INDEX
!    NCDEC IS COMMON DIMENSION OF AMAT & BMAT; NCACT IS ACTUAL LIMIT
!
!    I.E.   NCADEC IS NUMBER OF COLUMNS OF A DECLARED
!           NCBDEC IS NUMBER OF COLUMNS OF B DECLARED
!           NCDEC IS THE NUMBER OF ROWS IN BOTH A AND B DECLARED
!
!    MODIFICATION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES MACKINNON,
!    QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!--------------------------------------------------------------------------}
Var 
    i, j, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, ncc32r,
    nend: integer;  SUM:double;
Begin

  For  i:=1 to ncaact do
  begin
{   FIND NUMBER OF GROUPS OF SIZE 32, 16...}
    nend:=IMIN(i,ncact);

{   THIS ADJUSTMENT IS REQUIRED WHEN NCACT IS LESS THAN NCAACT }
    ncc32:=nend Div 32;
    ncc32r:=nend - 32*ncc32;
    ncc16:=ncc32r Div 16;
    ncc16r:=ncc32r - 16*ncc16;
    ncc8:=ncc16r Div 8;
    ncc8r:=ncc16r - 8*ncc8;
    ncc4:=ncc8r Div 4;
    ncc4r:=ncc8r - 4*ncc4;

{   FIND ENTRY IN MATRIX C }
    For j:=1 to ncbact do
    begin
      sum:=zero;
      k:=0;
      IF ncc32 > 0 THEN
        For kk:=1 to ncc32 do
        begin
          k:=k+32;
          sum:=sum + amat^[k-31,i]*bmat^[k-31,j]+amat^[k-30,i]*bmat^[k-30,
            j]+amat^[k-29,i]*bmat^[k-29,j]+amat^[k-28,i]*bmat^[k-28,j]+
            amat^[k-27,i]*bmat^[k-27,j]+amat^[k-26,i]*bmat^[k-26,j]+
            amat^[k-25,i]*bmat^[k-25,j]+amat^[k-24,i]*bmat^[k-24,j];
          sum:=sum + amat^[k-23,i]*bmat^[k-23,j]+amat^[k-22,i]*bmat^[k-22,
            j]+amat^[k-21,i]*bmat^[k-21,j]+amat^[k-20,i]*bmat^[k-20,j]+
            amat^[k-19,i]*bmat^[k-19,j]+amat^[k-18,i]*bmat^[k-18,j]+
            amat^[k-17,i]*bmat^[k-17,j]+amat^[k-16,i]*bmat^[k-16,j];
          sum:=sum + amat^[k-15,i]*bmat^[k-15,j]+amat^[k-14,i]*bmat^[k-14,
            j]+amat^[k-13,i]*bmat^[k-13,j]+amat^[k-12,i]*bmat^[k-12,j]+
            amat^[k-11,i]*bmat^[k-11,j]+amat^[k-10,i]*bmat^[k-10,j]+
            amat^[k-9,i]*bmat^[k-9,j]+amat^[k-8,i]*bmat^[k-8,j];
          sum:=sum + amat^[k-7,i]*bmat^[k-7,j]+amat^[k-6,i]*bmat^[k-6,j]+
            amat^[k-5,i]*bmat^[k-5,j]+amat^[k-4,i]*bmat^[k-4,j]+amat^[k-3,
            i]*bmat^[k-3,j]+amat^[k-2,i]*bmat^[k-2,j]+amat^[k-1,i]*
            bmat^[k-1,j]+amat^[k,i]*bmat^[k,j]
        end;

      IF ncc16 > 0 THEN
        For kk:=1 to ncc16 do
        begin
          k:=k+16;
          sum:=sum + amat^[k-15,i]*bmat^[k-15,j]+amat^[k-14,i]*bmat^[k-14,
            j]+amat^[k-13,i]*bmat^[k-13,j]+amat^[k-12,i]*bmat^[k-12,j]+
            amat^[k-11,i]*bmat^[k-11,j]+amat^[k-10,i]*bmat^[k-10,j]+
            amat^[k-9,i]*bmat^[k-9,j]+amat^[k-8,i]*bmat^[k-8,j];
          sum:=sum + amat^[k-7,i]*bmat^[k-7,j]+amat^[k-6,i]*bmat^[k-6,j]+
            amat^[k-5,i]*bmat^[k-5,j]+amat^[k-4,i]*bmat^[k-4,j]+amat^[k-3,
            i]*bmat^[k-3,j]+amat^[k-2,i]*bmat^[k-2,j]+amat^[k-1,i]*
            bmat^[k-1,j]+amat^[k,i]*bmat^[k,j]
        end;

      IF ncc8 > 0 THEN
        For kk:=1 to ncc8 do
        begin
          k:=k+8;
          sum:=sum + amat^[k-7,i]*bmat^[k-7,j]+amat^[k-6,i]*bmat^[k-6,j]+
            amat^[k-5,i]*bmat^[k-5,j]+amat^[k-4,i]*bmat^[k-4,j]+amat^[k-3,
            i]*bmat^[k-3,j]+amat^[k-2,i]*bmat^[k-2,j]+amat^[k-1,i]*
            bmat^[k-1,j]+amat^[k,i]*bmat^[k,j]
        end;

      IF ncc4 > 0 THEN
        For kk:=1 to ncc4 do
        begin
          k:=k+4;
          sum:=sum + amat^[k-3,i]*bmat^[k-3,j]+amat^[k-2,i]*bmat^[k-2,j]+
              amat^[k-1,i]*bmat^[k-1,j]+amat^[k,i]*bmat^[k,j]
        end;

      IF ncc4r > 0 THEN
        For kk:=1 to ncc4r do
        begin
          k:=k+1;
          sum:=sum + amat^[k,i]*bmat^[k,j]
        end;
      cmat^[i,j]:=sum
    end
  end;
End; {utbmul}

Procedure uvmul (nradec, nraact, ncdec, ncact:integer; amat:pMat; bvec, cvec:pVec);
{-------------------------------------------------------------------------------
!    FEB. 8, 1991
!
!    MATRIX-VECTOR MULTIPLICATION:  AB=C  WHERE A IS UPPER TRIANGULAR
!
!    VERSION WITH INNER LOOP UNROLLED TO DEPTHS 32, 16, 8 AND 4
!    EACH ROW OF MATRIX A IS SAVED AS A COLUMN BEFORE USE.
!
!    NRADEC IS 1ST DIM. OF AMAT; NRAACT IS ACTUAL LIMIT FOR 1ST INDEX
!    NCDEC IS COMMON DIMENSION OF AMAT & BVEC; NCACT IS ACTUAL LIMIT
!
!    I.E. NRADEC IS THE NUMBER OF ROWS OF A DECLARED
!         NCDEC IS THE COMMON DECLARED DIMENSION (COLUMNS OF A AND ROWS OF B)
!
!    MODIFICATION OF THE MATRIX MULTIPLIER DONATED BY PROF. JAMES MACKINNON,
!    QUEEN'S UNIVERSITY, KINGSTON, ONTARIO, CANADA
!------------------------------------------------------------------------------} 
Var
    i, k, kk, ncc4, ncc4r, ncc8, ncc8r, ncc16, ncc16r, ncc32, ncc32r:integer;
    SUM:double;
Begin
  For i:=1 to nraact do
  begin

{   FIND NUMBER OF GROUPS OF SIZE 32, 16...}

    ncc32:=(ncact - (i-1)) Div 32;
    ncc32r:=(ncact - (i-1))-32*ncc32;
    ncc16:=ncc32r Div 16;
    ncc16r:=ncc32r - 16*ncc16;
    ncc8:=ncc16r Div 8;
    ncc8r:=ncc16r - 8*ncc8;
    ncc4:=ncc8r Div 4;
    ncc4r:=ncc8r - 4*ncc4;

{   FIND ENTRY FOR VECTOR C }

    sum:=zero;
    k:=i-1;
    IF ncc32 > 0 THEN
      For kk:=1 to ncc32 do
      begin
        k:=k+32;
        sum:=sum + amat^[i,k-31]*bvec^[k-31] + amat^[i,k-30]*bvec^[k-30] +
          amat^[i,k-29]*bvec^[k-29] + amat^[i,k-28]*bvec^[k-28] +
          amat^[i,k-27]*bvec^[k-27] + amat^[i,k-26]*bvec^[k-26] +
          amat^[i,k-25]*bvec^[k-25] + amat^[i,k-24]*bvec^[k-24];
        sum:=sum + amat^[i,k-23]*bvec^[k-23] + amat^[i,k-22]*bvec^[k-22] +
          amat^[i,k-21]*bvec^[k-21] + amat^[i,k-20]*bvec^[k-20] +
          amat^[i,k-19]*bvec^[k-19] + amat^[i,k-18]*bvec^[k-18] +
          amat^[i,k-17]*bvec^[k-17] + amat^[i,k-16]*bvec^[k-16];
        sum:=sum + amat^[i,k-15]*bvec^[k-15] + amat^[i,k-14]*bvec^[k-14] +
          amat^[i,k-13]*bvec^[k-13] + amat^[i,k-12]*bvec^[k-12] +
          amat^[i,k-11]*bvec^[k-11] + amat^[i,k-10]*bvec^[k-10] +
          amat^[i,k-9]*bvec^[k-9] + amat^[i,k-8]*bvec^[k-8];
        sum:=sum + amat^[i,k-7]*bvec^[k-7] + amat^[i,k-6]*bvec^[k-6] +
          amat^[i,k-5]*bvec^[k-5] + amat^[i,k-4]*bvec^[k-4] + amat^[i,k-3]*bvec^[k-3] +
          amat^[i,k-2]*bvec^[k-2] + amat^[i,k-1]*bvec^[k-1] + amat^[i,k]*bvec^[k]
      end;

    IF ncc16 > 0 THEN
      For kk:=1 to ncc16 do
      begin
        k:=k+16;
        sum:=sum + amat^[i,k-15]*bvec^[k-15] + amat^[i,k-14]*bvec^[k-14] +
          amat^[i,k-13]*bvec^[k-13] + amat^[i,k-12]*bvec^[k-12] +
          amat^[i,k-11]*bvec^[k-11] + amat^[i,k-10]*bvec^[k-10] +
          amat^[i,k-9]*bvec^[k-9] + amat^[i,k-8]*bvec^[k-8];
        sum:=sum + amat^[i,k-7]*bvec^[k-7] + amat^[i,k-6]*bvec^[k-6] +
          amat^[i,k-5]*bvec^[k-5] + amat^[i,k-4]*bvec^[k-4] + amat^[i,k-3]*bvec^[k-3] +
          amat^[i,k-2]*bvec^[k-2] + amat^[i,k-1]*bvec^[k-1] + amat^[i,k]*bvec^[k]
      end;

    IF ncc8 > 0 THEN
      For kk:=1 to ncc8 do
      begin
        k:=k+8;
        sum:=sum + amat^[i,k-7]*bvec^[k-7] + amat^[i,k-6]*bvec^[k-6] +
          amat^[i,k-5]*bvec^[k-5] + amat^[i,k-4]*bvec^[k-4] + amat^[i,k-3]*bvec^[k-3] +
          amat^[i,k-2]*bvec^[k-2] + amat^[i,k-1]*bvec^[k-1] + amat^[i,k]*bvec^[k]
      end;

    IF ncc4 > 0 THEN
      For kk:=1 to ncc4 do
      begin
        k:=k+4;
        sum:=sum + amat^[i,k-3]*bvec^[k-3] + amat^[i,k-2]*bvec^[k-2] +
                amat^[i,k-1]*bvec^[k-1] + amat^[i,k]*bvec^[k]
      end;

    IF ncc4r > 0 THEN
      For kk:=1 to ncc4r do
      begin
        k:=k+1;
        sum:=sum + amat^[i,k]*bvec^[k]
      end;

    cvec^[i]:=sum
  end; {i loop}

End; {uvmul}

END. {Unit Nnes}