{*******************************************************************
* SMOOTHING AN ARRAY OF N ORDINATES Y's (ASCENDING ORDER ABCISSAS) *
*           using Savitzky-Golay filter coefficients               *
* ---------------------------------------------------------------- *
* Description:                                                     *
* This program uses the procedure SAVGOL for smoothing an array of *
* given ordinates (y's) that are in order of increasing abscissas  *
* (x's), but without using the abscissas themselves supposed to be *
* equally spaced. It first calculates the filter coefficients and  *
* then applies the numerical filter to low-pass filter the data.   *
* The user-specified parameter are: nl, nr and m  to enter "the    *
* amount of smoothing", given roughly as the number of points over *
* which the data should be smoothed.                               *
* ---------------------------------------------------------------- *
* SAMPLE RUN:                                                      *
* Input data file contains:                                        *
* 1024                                                             *
* 0.00000000000000E+0000 7.50000000000000E-0001                    *
* 9.21288168207468E-0003 7.77368637910513E-0001                    *
* 1.84257633641494E-0002 8.34466556277221E-0001                    *
* 2.76386450462383E-0002 9.03071871110114E-0001                    *
* 3.68515267282987E-0002 9.92958153417021E-0001                    *
* 4.60644084103592E-0002 1.09195646826811E+0000                    *
* 5.52772900924197E-0002 1.15230452277865E+0000                    *
* 6.44901717745370E-0002 1.06763022290215E+0000                    *
* 7.37030534565974E-0002 1.34541171127239E+0000                    *
* 8.29159351386579E-0002 1.48611048393104E+0000                    *
* 9.21288168207184E-0002 1.09349703210864E+0000                    *
* 1.01341698502779E-0001 1.72386602840743E+0000                    *
* 1.10554580184839E-0001 1.14317464708984E+0000                    *
* ---------------------- ----------------------                    *
* 9.37871355209791E+0000 2.43969819122867E+0001                    *
* 9.38792643377383E+0000 2.42468007203424E+0001                    *
* 9.39713931544975E+0000 2.42436619192304E+0001                    *
* 9.40635219712567E+0000 2.42829449073179E+0001                    *
* 9.41556507880159E+0000 2.42980085689633E+0001                    *
* 9.42477796047751E+0000 2.43119449022633E+0001                    *
*                                                                  *
* Output file contains (here nl=nr=5, m=4):                        *
*                                                                  *
*       Time          Y        Smoothed Y                          *
* ----------------------------------------                         *
*     0.000000     0.750000     0.504764                           *
*     0.009213     0.777369     0.739296                           *
*     0.018426     0.834467     0.890450                           *
*     0.027639     0.903072     0.973755                           *
*     0.036852     0.992958     0.966455                           *
*     0.046064     1.091956     1.028790                           *
*     0.055277     1.152305     1.162248                           *
*     0.064490     1.067630     1.173897                           *
*     0.073703     1.345412     1.281085                           *
*     0.082916     1.486110     1.351451                           *
*     0.092129     1.093497     1.426208                           *
*     0.101342     1.723866     1.348183                           *
*     0.110555     1.143175     1.332568                           *
*     --------     --------     --------                           *
*     9.341862    24.185782    24.167925                           *
*     9.351075    24.135968    24.154880                           *
*     9.360288    24.252124    24.215144                           *
*     9.369501    24.183153    24.268185                           *
*     9.378714    24.396982    24.279420                           *
*     9.387926    24.246801    24.246801                           *
*     9.397139    24.243662    24.243662                           *
*     9.406352    24.282945    24.282945                           *
*     9.415565    24.298009    24.298009                           *
*     9.424778    24.311945    24.311945                           *
*                                                                  *
* Note: the last nr points are unchanged.                          *
* ---------------------------------------------------------------- *
* Reference:  "Numerical Recipes By W.H. Press, B. P. Flannery,    *
*              S.A. Teukolsky and W.T. Vetterling, Cambridge       *
*              University Press, 1986 - 1992" [BIBLI 08].          *
*                                                                  *
*                              TPW Release By J-P Moreau, Paris    *
*                                 (with graphic option)            *
*                                   (www.jpmoreau.fr)              *
*******************************************************************}
PROGRAM tsavgol;

Uses  WinCrtMy, Type_def, Graph_2d;

Const NP=50; MMAX = 6;

Type
      pVEC1 = ^VEC1;
      VEC1 = Array[1..NP] of REAL;
      pIVEC = ^IVEC;
      IVEC = Array[1..NP] of Integer;

      pMAT = ^MAT;
      MAT  = Array[1..MMAX+1,1..MMAX+1] of REAL;
      pVEC2 = ^VEC2;
      VEC2  = Array[1..MMAX+1] of REAL;
      IVEC2 = Array[1..MMAX+1] of Integer;

Var
      signal, ysave: RV;  {see Type_def.pas}
      c: pVEC1;
      index: pIVEC;
      i,itmp,j,m,ndata,nl,nr: Integer;
      dt,t,tbegin,temp,tend: REAL;

      fp_in,fp_out: TEXT;

Function IMin(ia,ib:Integer):Integer;
Begin
  if ia<=ib then IMin:=ia
            else IMin:=ib
End;

Function IPower(x:REAL; n:integer): REAL;
{Integer power n of x (n>=0) }
var i: Integer; result: REAL;
begin
  result :=1.0;
  if n=0 then
  begin
    IPower:=result;
    exit;
  end
  else
    for i:=1 to n do result := x * result;
  IPower:=result;
end;

  Procedure LUDCMP(Var A:pMAT; N,NP: Integer;Var INDX: IVEC2; Var D,CODE:Integer); Forward;
  Procedure LUBKSB(Var A:pMAT; N,NP:Integer; INDX:IVEC2; Var B:pVEC2); Forward;

Procedure savgol(VAR c:pVEC1; np,nl,nr,ld,m:Integer);
{-------------------------------------------------------------------------------------------- 
 USES lubksb,ludcmp given below. 
 Returns in c(np), in wrap-around order (see reference) consistent with the argument respns 
 in routine convlv, a set of Savitzky-Golay filter coefficients. nl is the number of leftward 
 (past) data points used, while nr is the number of rightward (future) data points, making 
 the total number of data points used nl +nr+1. ld is the order of the derivative desired 
 (e.g., ld = 0 for smoothed function). m is the order of the smoothing polynomial, also 
 equal to the highest conserved moment; usual values are m = 2 or m = 4. 
--------------------------------------------------------------------------------------------}
Label return;
Var   d,icode,imj,ipj,j,k,kk,mm: Integer;
      indx:IVEC2;
      fac,sum: REAL;
      a:pMAT;
      b: pVEC2;
Begin
  New(a); New(b);
  if (np<nl+nr+1)or(nl<0)or(nr<0)or(ld>m)or(m>MMAX)or(nl+nr<m) then
  begin
    writeln(' Bad args in savgol.');
    goto return
  end;

  For i:=1 to MMAX+1 do
  begin
    for j:=1 to MMAX+1 do a^[i,j]:=0.0;
    b^[i]:=0.0;
    indx[i]:=0
  end;

  for ipj:=0 to 2*m do          {Set up the normal equations of the desired leastsquares fit}
  begin
    sum:=0.0;
    if ipj=0 then sum:=1.0;
    for k:=1 to nr do sum:=sum + IPower(1.0*k,ipj);
    for k:=1 to nl do sum:=sum + IPower(-1.0*k,ipj);
    mm:=IMin(ipj,2*m-ipj);
    imj:=-mm;
    repeat
      a^[1+(ipj+imj) div 2,1+(ipj-imj) div 2]:=sum;
      Inc(imj,2)
    until imj>mm; 
  end;

  ludcmp(a,m+1,MMAX+1,indx,d,icode);    {Solve them: LU decomposition}
 
  for j:=1 to m+1 do b^[j]:=0.0;
  b^[ld+1]:=1.0;   {Right-hand side vector is unit vector, depending on which derivative we want}

  lubksb(a,m+1,MMAX+1,indx,b);        {Backsubstitute, giving one row of the inverse matrix}
  
  for kk:=1 to np do c^[kk]:=0.0;     {Zero the output array (it may be bigger than the number}
                                      {of coefficients}

  for k:=-nl to nr do                 {Each Savitzky-Golay coefficient is the dot product}
  begin                               {of powers of an integer with the inverse matrix row}
    sum:=b^[1];                       
    fac:=1.0;
    for mm:=1 to m do
    begin
      fac:=fac*k;
      sum:=sum+b^[mm+1]*fac
    end; 
    kk:=((np-k) mod np) + 1;            {Store in wrap-around order}
    c^[kk]:=sum
  end;
return: Dispose(a); Dispose(b)
End;

{**************************************************************
* Given an N x N matrix A, this routine replaces it by the LU *
* decomposition of a rowwise permutation of itself. A and N   *
* are input. INDX is an output vector which records the row   *
* permutation effected by the partial pivoting; D is output   *
* as -1 or 1, depending on whether the number of row inter-   *
* changes was even or odd, respectively. This routine is used *
* in combination with LUBKSB to solve linear equations or to  *
* invert a matrix. Return code is 1, if matrix is singular.   *
**************************************************************}
Procedure LUDCMP(Var A:pMAT; N,NP: Integer;Var INDX: IVEC2; Var D,CODE:Integer);
Label return;
Const NMAX=100; TINY=1E-12;
Type  pV1 = ^V1;
       V1 = Array[1..NMAX] of REAL;
Var   AMAX,DUM, SUM:REAL;
      VV:pV1;
      I,IMAX,J,K: Integer;
Begin

 New(VV);

 D:=1; CODE:=0;

 For I:=1 to N do
 begin
   AMAX:=0.0;
   For J:=1 to N do
     IF ABS(A^[I,J]) > AMAX THEN AMAX:=ABS(A^[I,J]);
   IF AMAX < TINY THEN
   begin
     CODE := 1;
     GOTO return
   end;
   VV^[I] := 1.0 / AMAX
 end;

 For J:=1 to N do
 begin
   For I:=1 to J-1 do
   begin
     SUM := A^[I,J];
     For K:=1 to I-1 do
       SUM := SUM - A^[I,K]*A^[K,J];
     A^[I,J] := SUM
   end;
   AMAX := 0.0;
   For I:=J to N do
   begin
     SUM := A^[I,J];
     For K:=1 to J-1 do
       SUM := SUM - A^[I,K]*A^[K,J];
     A^[I,J] := SUM;
     DUM := VV^[I]*ABS(SUM);
     IF DUM >= AMAX THEN
     begin
       IMAX := I;
       AMAX := DUM
     end
   end;  
   
   IF J <> IMAX THEN
   begin
     For K:=1 to N do
     begin
       DUM := A^[IMAX,K];
       A^[IMAX,K] := A^[J,K];
       A^[J,K] := DUM
     end;
     D := -D;
     VV^[IMAX] := VV^[J]
   end;

   INDX[J] := IMAX;
   IF ABS(A^[J,J]) < TINY THEN A^[J,J] := TINY;

   IF J <> N THEN
   begin
     DUM := 1.0 / A^[J,J];
     For I:=J+1 to N do A^[I,J] := A^[I,J]*DUM
   end 
 end; {j loop}

 return: Dispose(VV)
 End;


{*****************************************************************
* Solves the set of N linear equations A . X := B.  Here A is    *
* input, not as the matrix A but rather as its LU decomposition, *
* determined by the routine LUDCMP. INDX is input as the permuta-*
* tion vector returned by LUDCMP. B is input as the right-hand   *
* side vector B, and returns with the solution vector X. A, N and*
* INDX are not modified by this routine and can be used for suc- *
* cessive calls with different right-hand sides. This routine is *
* also efficient for plain matrix inversion.                     *
*****************************************************************}
Procedure LUBKSB(Var A:pMAT; N,NP:Integer; INDX:IVEC2; Var B:pVEC2);
Var SUM: REAL;
    I,II,J,LL:Integer;
Begin

 II := 0;

 For I:=1 to N do
 begin
   LL := INDX[I];
   SUM := B^[LL];
   B^[LL] := B^[I];
   IF II <> 0 THEN
     For J:=II to I-1 do
       SUM := SUM - A^[I,J]*B^[J]
   ELSE IF SUM <> 0.0 THEN
     II := I;
   B^[I] := SUM
 end;

 For I:=N Downto 1 do
 begin
   SUM := B^[I];
   IF I < N THEN
     For J:=I+1 to N do
       SUM := SUM - A^[I,J]*B^[J];
   B^[I] := SUM / A^[I,I]
 end

End;


{main program}
BEGIN

  WinCrtInit(' TSAVGOL');

  New(signal); New(ysave); New(c); New(index);

  For i:=1 to Size do signal^[i]:=0.0;

  {open input and output file}
  Assign(fp_in,'smooth.dat'); Reset(fp_in);
  Assign(fp_out,'tsavgol.lst'); Rewrite(fp_out);

  {read number of input signal points in input file}
  read(fp_in, ndata);

  {read ndata couples T(i), Y(i) in input data file}
  for i:=1 to ndata do
  begin
    readln(fp_in, temp, signal^[i]);
    if i = 1 then  tbegin:=temp;
    if i = ndata then  tend:=temp
  end;
  close(fp_in);

  For i:=1 to size do ysave^[i]:=signal^[i];   {save unsmoothed signal}

  nl:=5; nr:=5; m:=4;                          {see savgol}
  
{ seek shift index for given case nl, nr, m (see savgol) } 
  index^[1]:=0;
{ example: case nl=nr=5
  index(2)=-1; index(3)=-2; index(4)=-3; index(5)=-4; index(6)=-5 }
  j:=3;
  for i:=2 to nl+1 do
  begin
    index^[i]:=i-j;
    j:=j+2
  end;
{ index(7)= 5; index(8)= 4; index(9)= 3; index(10)=2; index(11)=1 }
  j:=2;
  for i:=nl+2 to nl+nr+1 do
  begin
    index^[i]:=i-j;
    j:=j+2
  end;

{ calculate Savitzky-Golay filter coefficients }
  savgol(c,nl+nr+1,nl,nr,0,m);

  writeln;
  writeln(' Number of left points .......: ', nl);
  writeln(' Number of right points ......: ', nr);
  writeln(' Order of smoothing polynomial: ', m);
  writeln;
  writeln(' Savitzky-Golay Filter Coefficients:');
  for i:=1 to nl+nr+1 do write(c^[i]:10:6);
  writeln;

{ Apply filter to input data }
  for i:=1 to ndata-nr do
  begin
    signal^[i]:=0.0;
    for j:=1 to nl+nr+1 do
    begin
      itmp:=i+index^[j];
      if itmp>0 then        {skip left points that do not exist}
        signal^[i]:=signal^[i]+c^[j]*ysave^[itmp]
    end
  end;

{ write results to output file }
  dt:=(tend-tbegin)/(ndata-1);
  t:=tbegin-dt;
  writeln(fp_out,'      Time          Y        Smoothed Y ');
  writeln(fp_out,'----------------------------------------');
  For i:=1 to ndata do
  begin
    t:=t+dt;
    writeln(fp_out, '   ', t:10:6, '   ', ysave^[i]:10:6, '   ', signal^[i]:10:6)
  end;
  close(fp_out);

  writeln;
  writeln(' Results in file tsavgol.lst.');

  Dispose(c); Dispose(index);

  writeln;
  write(' Do you want a graph of signals (y/n) ? '); readln(rep);

  if rep='y' then
  begin
    clrscr;
    {manual scaling}
    Ech_auto:=False;              {see Graph_2d.pas}
    X_mini:=0.0; X_maxi:=tend;
    Y_mini:=0.0; Y_maxi:=20.0;
    CourbeXY(CrtDc,ndata,5,ysave,tbegin,tend);
    Legendes(CrtDC,'Input Data','X','Y');
    CourbeXY(CrtDc,ndata,3,signal,tbegin,tend);
    Legendes(CrtDC,'Smoothed Data','X','Y');
    SortieGraphique
  end;
  
  Dispose(signal); Dispose(ysave); 
  DoneWinCrt

END.

{end of file tsavgol.pas}