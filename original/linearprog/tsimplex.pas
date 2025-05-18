{***************************************************************
*           LINEAR PROGRAMMING: THE SIMPLEX METHOD             *
* ------------------------------------------------------------ *
* SAMPLE RUN:                                                  *
* Maximize z = x1 + x2 + 3x3 -0.5x4 with conditions:           * 
*          x1  + 2x3 <= 740                                    *
*          2x2 - 7x4 <= 0                                      *
*          x2  - x3 + 2x4 >= 0.5                               *
*          x1 + x2 + x3 +x4 = 9                                *
*          and all x's >=0.                                    *
*                                                              *
* Number of variables in E.F.: 4                               *
* Number of <= inequalities..: 2                               *
* Number of >= inequalities..: 1                               *
* Number of = equalities.....: 1                               *
* Input Economic Function:                                     *
* Coefficient # 1: 1                                           *
* Coefficient # 2: 1                                           *
* Coefficient # 3: 3                                           *
* Coefficient # 4: -0.5                                        *
* Constant term..: 0                                           *
* Input constraint # 1:                                        *
* Coefficient # 1: 1                                           *
* Coefficient # 2: 0                                           *
* Coefficient # 3: 2                                           *
* Coefficient # 4: 0                                           *
* Constant term..: 740                                         *
* Input constraint # 2:                                        *
* Coefficient # 1: 0                                           *
* Coefficient # 2: 2                                           *
* Coefficient # 3: 0                                           *
* Coefficient # 4: -7                                          *
* Constant term..: 0                                           *
* Input constraint # 3:                                        *
* Coefficient # 1: 0                                           *
* Coefficient # 2: 1                                           *
* Coefficient # 3: -1                                          *
* Coefficient # 4: 2                                           *
* Constant term..: 0.5                                         *
* Input constraint # 4:                                        *
* Coefficient # 1: 1                                           *
* Coefficient # 2: 1                                           *
* Coefficient # 3: 1                                           *
* Coefficient # 4: 1                                           *
* Constant term..: 9                                           *
*                                                              *
* Input Table:                                                 *
*    0.00    1.00    1.00    3.00   -0.50                      *
*  740.00   -1.00    0.00   -2.00    0.00                      *
*    0.00    0.00   -2.00    0.00    7.00                      *
*    0.50    0.00   -1.00    1.00   -2.00                      *
*    9.00   -1.00   -1.00   -1.00   -1.00                      *
*                                                              *
* Maximum of E.F. =    17.02500                                *
*  X1 =    0.000000                                            *
*  X2 =    3.325000                                            *
*  X3 =    4.725000                                            *
*  X4 =    0.950000                                            *
*                                                              *
* ------------------------------------------------------------ *
* Reference: "Numerical Recipes By W.H. Press, B. P. Flannery, *
*             S.A. Teukolsky and W.T. Vetterling, Cambridge    *
*             University Press, 1986" [BIBLI 08].              *
*                                                              *
*                       TPW Release 1.0 By J-P Moreau, Paris   *
*                                (www.jpmoreau.fr)             *
***************************************************************} 
Program Test_Simplex;

Uses WinCrt;

Label 3;

Const
       MMAX=25; NMAX=25;

Type
       MAT  = Array[1..MMAX,1..NMAX] of REAL;
       VEC  = Array[1..MMAX] of REAL;
       IVEC =  Array[1..MMAX] of Integer;

Var
       A: MAT;
       IPOSV, IZROV: IVEC;
       i,j,ICASE,N,M,M1,M2,M3: Integer;
       R: REAL;

       Procedure simp1(var a:MAT; mm:integer; ll:IVEC; nll, iabf: integer; var kp: integer;
                       var bmax:REAL);  Forward;
       Procedure simp2(var a:MAT; m, n:integer; l2:IVEC; nl2:integer; var ip:integer;
                       kp:integer; var q1:REAL); Forward;
       Procedure simp3(var a:MAT; i1,k1,ip,kp:integer); Forward;


Procedure simplx(var a:MAT; m, n, m1, m2, m3: Integer; var icase:Integer; var izrov, iposv:IVEC);
 {----------------------------------------------------------------------------------------- 
 USES simp1,simp2,simp3. 
 Simplex method for linear programming. Input parameters a, m, n, m1, m2, and m3, and 
 output parameters a, icase, izrov, and iposv are described above (see reference). 
 Constants: MMAX is the maximum number of constraints expected; NMAX is the maximum number 
 of variables expected; EPS is the absolute precision, which should be adjusted to the 
 scale of your variables.                                    
 -----------------------------------------------------------------------------------------}
Label 1,2,10,20,30, return;
Var
  i,ip,ir,is,k,kh,kp,m12,nl1,nl2: Integer;
  l1, l2, l3: IVEC; 
  bmax,q1,EPS: REAL;
Begin
  EPS:=1e-6;
  if m <> m1+m2+m3 then
  begin
    writeln(' Bad input constraint counts in simplx.');
    goto return
  end;	
  nl1:=n;
  for k:=1 to n do
  begin 
    l1[k]:=k;     {Initialize index list of columns admissible for exchange.}
    izrov[k]:=k   {Initially make all variables right-hand.}
  end; 
  nl2:=m;
  for i:=1 to m do
  begin 
    if a[i+1,1] < 0.0 then
    begin
      writeln(' Bad input tableau in simplx, Constants bi must be nonnegative.');
      goto return
    end;
    l2[i]:=i;
    iposv[i]:=n+i
{-------------------------------------------------------------------------------------------------
 Initial left-hand variables. m1 type constraints are represented by having their slackv ariable 
 initially left-hand, with no artificial variable. m2 type constraints have their slack 
 variable initially left-hand, with a minus sign, and their artificial variable handled implicitly 
 during their first exchange. m3 type constraints have their artificial variable initially 
 left-hand.     
-------------------------------------------------------------------------------------------------} 
  end;
  for i:=1 to m2 do l3[i]:=1;
  ir:=0;
  if m2+m3 = 0 then goto 30;  {The origin is a feasible starting solution. Go to phase two.}
  ir:=1;
  for k:=1 to n+1 do          {Compute the auxiliary objective function.}
  begin
    q1:=0.0;
    for i:=m1+1 to m do  q1 := q1 + a[i+1,k];
    a[m+2,k]:=-q1
  end; 
10: simp1(a,m+1,l1,nl1,0,kp,bmax);     {Find max. coeff. of auxiliary objective fn}
  if(bmax <= EPS) and (a[m+2,1] < -EPS) then
  begin
    icase:=-1;    {Auxiliary objective function is still negative and can’t be improved,}
    goto return   {hence no feasible solution exists.}
  end	
  else if (bmax <= EPS) and (a[m+2,1] <= EPS) then
{ Auxiliary objective function is zero and can’t be improved; we have a feasible starting vector. 
  Clean out the artificial variables corresponding to any remaining equality constraints by 
  goto 1’s and then move on to phase two by goto 30.  }
  begin
    m12:=m1+m2+1;
    if m12 <= m then
      for ip:=m12 to m do
        if iposv[ip] = ip+n then     {Found an artificial variable for an equalityconstraint.}
        begin
          simp1(a,ip,l1,nl1,1,kp,bmax);
          if bmax > EPS then goto 1; {Exchange with column corresponding to maximum} 
        end;                         {pivot element in row.} 
    ir:=0;
    m12:=m12-1;
    if m1+1 > m12 then goto 30;
    for i:=m1+1 to m1+m2 do           {Change sign of row for any m2 constraints}
      if l3[i-m1] = 1 then            {still present from the initial basis.}
        for k:=1 to n+1 do
          a[i+1,k] := -1.0 * a[i+1,k];
    goto 30                           {Go to phase two.} 
  end;

  simp2(a,m,n,l2,nl2,ip,kp,q1);       {Locate a pivot element (phase one). }
                                         
  if ip = 0 then                     {Maximum of auxiliary objective function is}
  begin                              {unbounded, so no feasible solution exists.}
    icase:=-1;                       
    goto return 
  end; 
1:  simp3(a,m+1,n,ip,kp);
{ Exchange a left- and a right-hand variable (phase one), then update lists.} 
  if iposv[ip] >= n+m1+m2+1 then     {Exchanged out an artificial variable for an}
  begin                              {equality constraint. Make sure it stays 
                                      out by removing it from the l1 list. } 
    for k:=1 to nl1 do
      if l1[k] = kp then goto 2;
2:  nl1:=nl1-1;
    for is:=k to nl1 do  l1[is]:=l1[is+1];
  end	
  else
  begin
    if iposv[ip] < n+m1+1 then goto 20;
    kh:=iposv[ip]-m1-n;
    if l3[kh] = 0 then goto 20; {Exchanged out an m2 type constraint.}
    l3[kh]:=0                   {If it’s the first time, correct the pivot column or the
                                 minus sign and the implicit artificial variable. } 
  end;   
  a[m+2,kp+1] := a[m+2,kp+1] + 1.0;
  for i:=1 to m+2 do  a[i,kp+1] := -1.0 * a[i,kp+1];
20: is:=izrov[kp];              {Update lists of left- and right-hand variables.}
  izrov[kp]:=iposv[ip];
  iposv[ip]:=is;
  if ir <> 0 then goto 10;       {if still in phase one, go back to 10.
  End of phase one code for finding an initial feasible solution. Now, in phase two, optimize it.} 
30: simp1(a,0,l1,nl1,0,kp,bmax); {Test the z-row for doneness.} 
  if bmax <= EPS then            {Done. Solution found. Return with the good news.}
  begin
    icase:=0;
    goto return 
  end;
  simp2(a,m,n,l2,nl2,ip,kp,q1);  {Locate a pivot element (phase two).} 
  if ip = 0 then            {Objective function is unbounded. Report and return.}
  begin
    icase:=1;
    goto return 
  end; 
  simp3(a,m,n,ip,kp);       {Exchange a left- and a right-hand variable (phase two),} 
  goto 20;                  {update lists of left- and right-hand variables and
                            {return for another iteration. }
return: End;

{The preceding routine makes use of the following utility subroutines: }

Procedure simp1(var a:MAT; mm:integer; ll:IVEC; nll, iabf: integer; var kp: integer;
                var bmax:REAL);
{ Determines the maximum of those elements whose index is contained in the supplied list 
  ll, either with or without taking the absolute value, as flagged by iabf. }
Label return;
Var
  k: integer; 
  test: REAL;
Begin   
  kp:=ll[1];
  bmax:=a[mm+1,kp+1];
  if nll < 2 then goto return; 
  for k:=2 to nll do
  begin
    if iabf = 0 then
      test:=a[mm+1,ll[k]+1]-bmax
    else
      test:=abs(a[mm+1,ll[k]+1])-abs(bmax);
    if test > 0.0 then
    begin 
      bmax:=a[mm+1,ll[k]+1];
      kp:=ll[k]
    end 
  end; 
return: End; 

Procedure simp2(var a:MAT; m, n:integer; l2:IVEC; nl2:integer; var ip:integer;
                kp:integer; var q1:REAL);
Label 2,6, return;
Var EPS: REAL;
    i,ii,k: integer;
    q,q0,qp: REAL;
Begin
  EPS:=1e-6;
{ Locate a pivot element, taking degeneracy into account.} 
  ip:=0;
  if nl2 < 1 then goto return;
  for i:=1 to nl2 do
    if a[i+1,kp+1] < -EPS then goto 2;
  goto return;  {No possible pivots. Return with message.} 
2: q1:=-a[l2[i]+1,1]/a[l2[i]+1,kp+1];
  ip:=l2[i];
  if i+1 > nl2 then goto return; 
  for i:=i+1 to nl2 do
  begin
    ii:=l2[i];
    if a[ii+1,kp+1] < -EPS then
    begin
      q:=-a[ii+1,1]/a[ii+1,kp+1];
      if q <  q1 then
      begin 
        ip:=ii;
        q1:=q
      end
      else if q = q1 then  {We have a degeneracy.}
      begin
	for k:=1 to n do
        begin
          qp:=-a[ip+1,k+1]/a[ip+1,kp+1];
          q0:=-a[ii+1,k+1]/a[ii+1,kp+1];
          if q0 <> qp then goto 6
        end; 
6:      if q0 < qp then ip:=ii
      end 
    end 
  end; 
return: End; 

Procedure simp3(var a:MAT; i1,k1,ip,kp:integer);
{ Matrix operations to exchange a left-hand and right-hand variable (see text).}
Var 
  ii,kk:integer; 
  piv:REAL;
Begin   
  piv:=1.0/a[ip+1,kp+1];
  if i1 >= 0 then
    for ii:=1 to i1+1 do
    begin
      if ii-1 <> ip then
      begin
        a[ii,kp+1] := a[ii,kp+1] * piv;
        for kk:=1 to k1+1 do
          if kk-1 <>  kp then
            a[ii,kk] := a[ii,kk] - a[ip+1,kk]*a[ii,kp+1]
      end
    end; 
  for kk:=1 to k1+1 do
    if kk-1 <> kp then a[ip+1,kk] :=-a[ip+1,kk]*piv;
  a[ip+1,kp+1]:=piv
End; 


{main program}
BEGIN

  writeln;
  write(' Number of variables in E.F.: '); readln(N);
  write(' Number of <= inequalities..: '); readln(M1);
  write(' Number of >= inequalities..: '); readln(M2);
  write(' Number of = equalities.....: '); readln(M3);

  M:=M1+M2+M3;   {Total number of constraints}

  for i:=1 to M+2 do
    for j:=1 to N+1 do
      A[i,j]:=0.0;

  writeln(' Input Economic Function:');
  for i:=2 to N+1 do
  begin
    write(' Coefficient #',i-1,': ');
    readln(A[1,i])
  end;
  write(' Constant term : ');
  readln(A[1,1]);
{ input constraints }    
  for i:=1 to M do
  begin
    writeln(' Input constraint #',i,':');
    for j:=2 to N+1 do
    begin
      write(' Coefficient #',j-1,': ');
      readln(R);
      A[i+1,j] := -R
    end;
    write(' Constant term : ');
    readln(A[i+1,1])
  end;

  writeln;
  writeln(' Input Table:');
  for i:=1 to M+1 do
  begin
    for j:=1 to N+1 do write(A[i,j]:8:2);
    writeln
  end;

  simplx(A,M,N,M1,M2,M3,ICASE,IZROV,IPOSV);

  if ICASE=0  then  {result ok.}
  begin
    writeln;
    writeln(' Maximum of E.F. = ', A[1,1]:12:6);

    for i:=1 to N do
    begin
      for j:=1 to M do
	if IPOSV[j] = i then
        begin
          writeln('  X',i,' = ', A[j+1,1]:12:6);
	  goto 3;
	end;
      writeln('  X',i,' = ', 0.0:12:6);
3:  end
  end
  else
    writeln(' No solution (error code = ', ICASE,').');
  
  writeln;
  ReadKey;
  DoneWinCrt

END.

{ end of file tsimplex.pas}