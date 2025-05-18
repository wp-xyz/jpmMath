{****************************************************
* This program finds the simple elements of a poly- *
* nomial fraction the denominator of which is given *
* as a product of irreducible polynomials.          *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUNS:                                      *
*                                                   *
*                      2x9 +1                       *
* EXAMPLE:  F(x) = ---------------                  *
*                  x3 (x2 +x +1)^2                  *
*                                                   *
*                                                   *
* FIND THE SIMPLE ELEMENTS OF A POLYNOMIAL FRACTION:*
*                                                   *
* Enter numerator: 2x9 +1                           *
*                                                   *
* Enter factor 1                                    *
* Power (0 to exit): 3                              *
* X2 coefficient: 0                                 *
* X  coefficient: 1                                 *
* Constant      : 0                                 *
* Enter factor 2                                    *
* Power (0 to exit): 2                              *
* X2 coefficient: 1                                 *
* X  coefficient: 1                                 *
* Constant      : 1                                 *
* Enter factor 3                                    *
* Power (0 to exit): 0                              *
*                                                   *
* Integer part:                                     *
*                                                   *
*       2                                           *
* 2.00 X  - 4.00 X + 2.00                           *
*                                                   *
* Factor 1                                          *
* A(1, 1) =          1.0000000000                   *
* A(1, 2) =         -2.0000000000                   *
* A(1, 3) =          1.0000000000                   *
*                                                   *
* Factor 2                                          *
* A(2, 1) =         -3.0000000000                   *
* B(2, 1) =          2.9999999999                   *
* A(2, 2) =          3.0000000001                   *
* A(2, 2) =          0.0000000002                   *
*                                                   *
*                                                   *
* So the result is:                                 *
*                                                   *
*                        1   -2   1      -3 +3x     *
* F(x) = (2x2 -4x +2) + (- + -- + --) + (------ +   *
*                        x   x2   x3     x2+x+1     *
*       <integer part>   <first factor>  <second    *
*                                                   *
*   3 + 0x                                          *
* ----------)                                       *
* (x2+x+1)^2                                        *
* factor>                                           *
*                                                   *
*                       TPW version by J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
*****************************************************
Explanations:
------------
1. Theory
   ------
Let us consider a polynomial fraction F(x)=P(x)/Q(x) the denominator
of which is given as a product of irreducible polynomials. Such a
polynomial has a degree=1, or =2 but with no real roots.

So the polynomial Q(x) must be given in the form:

   Q(x) = H1(x)^a1.H2(x)^a2...Hn(x)^an   (^ stands for power)

where the Hi(x) polynomials are irreducible.

We can then prove (not done here) that there exists a unique decomposition
of fraction F(x) in the form:

   F(x) = E(x) + B1(x) + B2(x) + ... + Bn(x)

where E(x) is a polynomial called "integer part" of F(x) and the Bi(x)
take the form:

           j=aj
   Bi(x) = Sum (Lij(x)/Hj(x))  
           j=1

where degree of Lij(x) polynomial strictly < degree of Hj(x).

So, for "1st type" elements, of degree=1, Lij(x) is a scalar.
For "2nd type" elements, Lij(x) has a degree <= 1 and has the form:

   aij + bij x.

Let us take an example:

             2x9 + 1
   F(x) = --------------
          x^3 (x2+x+1)^2

The irreducible polynomials are: H1(x)=x with a1=3
and H2(x)=x2+x+1 with a2=2.

So, according to the theory, F(x) can be written in a unique way as:
 
                 a11   a12   a13   a21 + b21x   a22 + b22x
   F(x) = E(x) + --- + --- + --- + ---------- + ----------             
                  x    x2    x3      x2+x+1     (x2+x+1)^2

                 [    B1(x)    ]   [        B2(x)         ]

Finding the simple elements of F(x) is nothing more than to find the
integer part E(x) and all the coefficients aij and bij corresponding
to the different factors.

2. The program
   -----------

The fraction F(x) is no longer given as two polynomials P(x)/Q(x) as
in the other programs but as 3 parameters:

   - numer, numerator polynomial of F(x),
   - nbfactor, number of denominator factors,
   - factor, table of ten variables of pascal type ONEFACTOR.

The type ONEFACTOR is defined as follows:  TYPE  ONEFACTOR = Record
                                                               a,b,c: REAL;
                                                               power: INTEGER
                                                             End;
where a, b, c are the real coefficients of Hi(x)=ax2+bx+c
(a=0 if degree=1) and power is ai.

The output of Function SimpleElements() is:

    - intpol, polynomial giving the "integer part" E(x),
    - U, pointer to a vector containing all the coefficients
      aij and bij, in the very order in which they appear
      in the above decomposition.

The rsult of the Function SimpleElements is TRUE if the calculation
is successful, else FALSE.

The program first verifies if input data are valid, chiefly the
irreducibility of denominator polynomial factors.

If data are correct, the program identifies initial fraction to the
sum of simple elements that must be obtained, according to theory,
for several x values with an arbitrary step of 0.25. So we obtain a
linear system to be solved by a classical Gauss-Jordan method in
double precision.

As output, the program gives the integer part and, for each factor,
the corresponding coefficients a(i,j) and b(i,j).

For accuracy problems, the size MX of the linear system is limited
to 20. 
               (adapted from above reference by J-P Moreau, Paris).
--------------------------------------------------------------------}
PROGRAM SIMPLE_ELEMENTS;
Uses WinCrt,Polynoms,PolFract;


CONST  
       MAXC = 40;      {maximum number of lines of square matrix M}

TYPE   ONEFACTOR = Record
                     a,b,c: REAL;
                     power: INTEGER
                   End;

       TENFACTORS = Array[1..10] of ONEFACTOR;


       pMAT   = ^MATRIX;
       MATRIX = Array[1..MAXC,1..MAXC] of REAL;

       pVEC   = ^VECTOR;
       VECTOR = Array[1..MAXC] of REAL;

VAR    numer,int: POLYNOM;
       factor   : TENFACTORS;
       nbfactor : INTEGER;
       U        : pVEC;
       i,j,k    : INTEGER;



  { P(X) * Q(X) = R(X) }
  Function MultPolynom(P,Q:POLYNOM;VAR R:POLYNOM): Boolean;
  Var i,j, n: INTEGER;
      u     : NUMBER;
  Begin
    MultPolynom:=TRUE;
    fillchar(R,sizeof(R),0);  {set R polynomial to zero}
    {verify that P and Q are not void}
    if (P.degree=0) and (P.coeff[0].value=0) then exit;
    if (Q.degree=0) and (Q.coeff[0].value=0) then exit;
    MultPolynom:=FALSE;
    R.degree:=P.degree+Q.degree;
    if R.degree>MAXPOL then exit;  {R degree is too big}
    for n:=0 to R.degree do
    begin
      if Not SetNumber(R.coeff[n],'0') then exit;
      for i:=0 to P.degree do
      begin
        j:=n-i;
        if (j>=0) and (j<=Q.degree) then
        begin
          if Not MultNumber(P.coeff[i],Q.coeff[j],u) then exit;
          if Not AddNumber(R.coeff[n],u,R.coeff[n]) then exit
        end
      end
    end;
    MultPolynom:=TRUE
  End;

  {Euclidian division of two polynoms}
  Function DivPolynom(P,Q:POLYNOM; VAR H,R:POLYNOM): Boolean;
  Var  i,j: INTEGER;
       u  : NUMBER;
  Begin
    DivPolynom:=FALSE;
    {The Q polynomial must be <> zero}
    if (Q.degree=0) and (Q.coeff[0].value=0) then exit;
    R:=P; H.degree:=P.degree - Q.degree;
    if H.degree<0 then
    begin
      H.degree:=0; if Not SetNumber(H.coeff[0],'0') then exit;
    end
    else
    begin
      for i:=H.degree downto 0 do
      begin
        if Not DivNumber(R.coeff[R.degree],Q.coeff[Q.degree],
                         H.coeff[i]) then exit;
        for j:=i to R.degree do
        begin
          if Not MultNumber(H.coeff[i],Q.coeff[j-i], u) then exit;
          u.p:=-u.p; u.value:=-u.value;
          if Not AddNumber(R.coeff[j],u, R.coeff[j]) then exit
        end;
        if R.degree > 0 then R.degree:=R.degree-1
      end;
      While (abs(R.coeff[R.degree].value) < SMALL) and
               (R.degree>0) do  R.degree:=R.degree-1
    end;
    DivPolynom:=TRUE
  End;

{******************************************
*  SOLVING A LINEAR MATRIX SYSTEM AX = B  *
*  with Gauss-Jordan method using full    *
*  pivoting at each step. During the pro- *
* cess, original AA matrix and BB vector  *
* are destroyed to spare storage location.*
* --------------------------------------- *
* INPUTS:    N    SIZE OF SYSTEM          *
*            AA   MATRIX N*N              *                                     
*            BB   VECTOR N                *                                     
* --------------------------------------- *                                     
* OUTPUS:    AA   INVERSE OF AA N*N       *                                     
*            DET  DETERMINANT OF AA       *                                     
*            BB   SOLUTION VECTOR N       *
* --------------------------------------- *
* NOTE: memory allocations for AA and BB  *
* must be made by the calling program.    *
*                                         *
*     Adapted from FORTRAN by J-P Moreau. *
******************************************}
Procedure MATINV(N:integer; VAR AA:pMAT; VAR BB:pVEC;VAR DET:DOUBLE);
CONST EPSMACH = 1.2e-16;
VAR   PC,PL,CS : pVEC;
      PV,PAV,temp,TT : DOUBLE;
      I,IK,J,JK,K : integer;

LABEL fin, fin1, fin2;

BEGIN
New(PC);New(PL);New(CS);                                            
{Initializations}                                                             
      DET := 1.0; 
      for I:=1 to N do
      begin                                                                
        PC^[I] := 0.0;                                                                
        PL^[I] := 0.0;                                                                
        CS^[I] := 0.0;
      end;              

{Main loop}                                                                             
      for K:=1 to N do
      begin                                                                  
{Searching greatest pivot}                                               
        PV:=AA^[K][K];                                                              
        IK:=K;                                                                    
        JK:=K;                                                                    
        PAV:=abs(PV);                                                            
        for I:=K to N do
          for J:=K to N do
          begin     
            temp := abs(AA^[I][J]);                                                        
            if temp > PAV then
            begin                                      
              PV:=AA^[I][J];
              PAV:=abs(PV);
              IK:=I;                                                              
              JK:=J
            end                                                               
          end;                                                                 
                                                                               
{Search terminated, the pivot is in location I=IK, J=JK.                     
 Memorizing pivot location:  }
                                           
        PC^[K]:=JK;
        PL^[K]:=IK;

{DETERMINANT DET is actualised                                              
 If DET=0, ERROR MESSAGE and STOP                                           
 Machine dependant EPSMACH equals here 1D-20  }                                         
        if IK<>K then DET:=-DET;                                                   
        if JK<>K then DET:=-DET;                                                   
        DET:=DET*PV;  
        temp:= abs(DET);
        if temp < EPSMACH then
        begin                                          
          writeln;                                                  
          writeln('  The determinant equals ZERO !!!');
          halt                                                                  
        end;

                                                                               
{POSITIONNING PIVOT IN K,K:}
        if IK<>K then                                                         
          for I:=1 to N do
          begin                                                              
{EXCHANGE LINES IK and K of matrix AA:}                                            
            TT:=AA^[IK][I];                                                         
            AA^[IK][I]:=AA^[K][I];
            AA^[K][I]:=TT
          end;                                                                 
                                                                           
        TT:=BB^[IK];                                                           
        BB^[IK]:=BB^[K];                                                      
        BB^[K]:=TT;

{PIVOT is at correct line}
        if JK<>K then                                                         
          for I:=1 to N do
          begin                                                              
{EXCHANGE COLUMNS JK and K of matrix AA:}                                          
            TT:=AA^[I][JK];                                                         
            AA^[I][JK]:=AA^[I][K];                                                    
            AA^[I][K]:=TT                                                          
          end;
                                                                           
{The PIVOT is at correct column.                                              
 and is located in K,K.                                                   
                                                                               
 Column K of matrix AA is stored in CS vector                             
 then column K is set to zero. }                                             
        for I:=1 to N do
        begin                                                                
          CS^[I]:=AA^[I][K];                                                         
          AA^[I][K]:= 0.0                                                          
        end;
                                                                               
        CS^[K]:= 0.0;                                                                
        AA^[K][K]:= 1.0;
{Line K of matrix AA is modified:}  
        temp := abs(PV);                                          
        if temp < EPSMACH then
        begin
          writeln;                                        
          writeln('  PIVOT TOO SMALL - STOP');                               
          exit
        end;                                                                   
        for I:=1 to N do                                                                
          AA^[K][I]:=AA^[K][I]/PV;                                                    
                                                                           
        BB^[K]:=BB^[K]/PV;
                                                                           
{Other lines of matrix AA are modified:}                                        
        for J:=1 to N do
        begin                                                                
          if J=K then goto fin;                                                  
          for I:=1 to N do
{Line J of matrix AA is modified:}                                            
            AA^[J][I]:=AA^[J][I]-CS^[J]*AA^[K][I];
          BB^[J]:=BB^[J]-CS^[J]*BB^[K];
        fin: end                                                                   
{Line K is ready.}
      end; {of K loop}

{MATRIX AA INVERSION IS DONE - REARRANGEMENT OF MATRIX AA
                                                                               
   EXCHANGE LINES                }                                                            
      for I:=N downto 1 do
      begin                                                               
        IK:=Round(PC^[I]);
        if IK=I then goto fin1;                                                   
{EXCHANGE LINES I and PC(I) of matrix AA:}                                         
        for J:=1 to N do
        begin                                                                
          TT:=AA^[I][J];                                                            
          AA^[I][J]:=AA^[IK][J];                                                      
          AA^[IK][J]:=TT                                                           
        end;
        TT:=BB^[I];                                                          
        BB^[I]:=BB^[IK];                                                    
        BB^[IK]:=TT;
{NO MORE EXCHANGE is NECESSARY                                                      
 Go to next line.}                                                  
      fin1: end; {of loop i=N downto }                                                                     
                                                                               
{EXCHANGE COLUMNS  }                                                          
      
      for J:=N downto 1 do
      begin                                                                         
        JK:=Round(PL^[J]);
        if JK=J then goto fin2;                                                   
{EXCHANGE COLUMNS J ET PL(J) of matrix AA: }                                       
        for I:=1 to N do
        begin                                                                
          TT:=AA^[I][J];                                                            
          AA^[I][J]:=AA^[I][JK];                                                      
          AA^[I][JK]:=TT
        end; 
{NO MORE EXCHANGE is NECESSARY                                                      
 Go to next column.}
      fin2: end;                                                                     
{REARRANGEMENT TERMINATED.  }
Dispose(PC);Dispose(PL);Dispose(CS)
END;   {of procedure MATINV }


{find the simple elements - See explanations above}
Function SimpleElements(numer:POLYNOM; nbfactor:integer;factor:TENFACTORS;
                          VAR intpol:POLYNOM;VAR U:pVEC): BOOLEAN;

CONST  MX   = 20;     {maximum size of system to solve with acceptable accuracy}
       STEP = 0.25;   {increment step for x}

VAR    i,j,k,h : INTEGER;
       x,s,t,w : REAL;
       deter   : DOUBLE;
       y,z     : NUMBER;
       ok      : BOOLEAN;
       M       : pMAT;
       mcol,mlig,ulig,vlig: integer;
       inter,denom,R: POLYNOM;
Begin
  SimpleElements:=FALSE; ulig:=0;
  New(M); {allocate memory for M matrix}
  if (nbfactor<1) or (nbfactor>10) then exit;
  {verify that each factor is not reductible}
  for i:=1 to nbfactor do
    with factor[i] do
      if (a<>0) and (a*a >= 4*b*c) then
      begin
        writeln(' One of the factors is reductible!');
        exit
      end;
  {verify that each factor is different from others}
  for i:=1 to nbfactor do
    for j:=1 to i-1 do
      if (factor[i].a=factor[j].a) and (factor[i].b=factor[j].b)
            and (factor[i].c=factor[j].c) then
      begin
        writeln(' At least two factors are identical!');
        exit
      end;
  {verify that system size < MX}
  ulig:=0;
  for i:=1 to nbfactor do
    With factor[i] do
      if a<>0 then
        ulig:=ulig+2*power
      else
        ulig:=ulig+power;
  if ulig>MX then
  begin
    writeln(' Denominator too big!');
    exit
  end;
  {verify that fraction is not reductible}
  for i:=1 to nbfactor do
    With factor[i] do
    begin
      if a<>0 then denom.degree:=2 else denom.degree:=1;
      denom.coeff[0].is_real:=TRUE; denom.coeff[0].value:=c;
      denom.coeff[1].is_real:=TRUE; denom.coeff[1].value:=b;
      denom.coeff[2].is_real:=TRUE; denom.coeff[2].value:=a;
      if Not DivPolynom(numer,denom,inter,R) then exit;
      if (R.degree=0) and (abs(R.coeff[0].value)<SMALL) then
      begin
        writeln(' Fraction is reductible!');
        exit
      end
    end;
  {calculate denominator polynomial}
  denom.degree:=0;
  denom.coeff[0].is_real:=TRUE; denom.coeff[0].value:=1;
  for i:=1 to nbfactor do
    With factor[i] do
    begin
      if a<>0 then inter.degree:=2 else inter.degree:=1;
      inter.coeff[0].is_real:=TRUE; inter.coeff[0].value:=c;
      inter.coeff[1].is_real:=TRUE; inter.coeff[1].value:=b;
      inter.coeff[2].is_real:=TRUE; inter.coeff[2].value:=a;
      for j:=1 to power do
        if Not MultPolynom(denom,inter,denom) then exit
    end;

  {calculate integer part}
  if Not DivPolynom(numer,denom,intpol,R) then exit;

  {create linear system MU = V to be solved}
  mlig:=ulig; mcol:=mlig; vlig:=ulig;
  x:=-vlig*STEP/2 -STEP;
  for i:=1 to vlig do
  begin
    {eliminate an x for which denominator =0}
    j:=0;
    Repeat
      x:=x+STEP; ok:=TRUE;
      for k:=1 to nbfactor do
        with factor[k] do
          if abs(c+x*(b+x*a)) < SMALL then ok:=FALSE
    Until ok;
    {fill first line of matrix M}
    for k:=1 to nbfactor do
      with factor[k] do
      begin
        y.is_real:=TRUE; y.value:=x;
        if Not EvaluatePolynom(denom,y,z) then exit;
        w:=z.value; t:=c+x*(b+a*x); s:=1;
        for h:=1 to power do
        begin
          s:=s*t; Inc(j); M^[i,j]:=w/s;
          if a<>0 then
          begin
            Inc(j); M^[i,j]:=w*x/s
          end
        end
      end;
      y.is_real:=TRUE; y.value:=x;
      if Not EvaluatePolynom(R,y,z) then exit;
      U^[i]:=z.value
    end;
    {solve linear system MU=V by Gauss-Jordan method}
    MATINV(mlig,M,U,deter);
    {free memory for M matrix}
    Dispose(M);
    {Process succesful} 
    SimpleElements:=TRUE
End;


{main program}
BEGIN
  New(U);  {allocate memory for U vector}
  writeln;
  writeln(' FIND THE SIMPLE ELEMENTS OF A POLYNOMIAL FRACTION:');
  writeln;
  if Not EnterPolynom(' Enter numerator: ',numer) then exit;
  writeln; i:=0;
  Repeat
    Inc(i);
    writeln(' Enter factor ',i);
    write(' Power (0 to exit): '); readln(factor[i].power);
    if factor[i].power<>0 then
    begin
      write(' X2 coefficient: '); readln(factor[i].a);
      write(' X  coefficient: '); readln(factor[i].b);
      write(' Constant      : '); readln(factor[i].c)
    end
  Until (i=0) or (factor[i].power=0);
  nbfactor:=i-1; writeln;
  if SimpleElements(numer,nbfactor,factor,int,U) then
  begin
    writeln(' Integer part:');
    DisplayPolynom(int); j:=0;
    writeln;
    for i:=1 to nbfactor do
    begin
      writeln; writeln(' Factor ',i);
      for k:=1 to factor[i].power do
      begin
        Inc(j);
        writeln(' A(',i,', ',k,') = ', U^[j]:22:10);
        if factor[i].a<>0 then
        begin
          Inc(j);
          writeln(' B(',i,', ',k,') = ', U^[j]:22:10)
        end
      end
    end
  end
  else
    writeln(' Error in finding simple elements.');
  writeln;
  Dispose(U); {free memory for U vector}
  Readkey;  
  DoneWinCrt
END.

{end of file simpelem.pas}