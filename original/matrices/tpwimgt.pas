{*************************************************************
* Calculate the smallest eigenvalue lambda of a real square  *
* matrix and the associated eigenvector. The method consists *
* to inverse the matrix by a classical Gauss-Jordan method,  *
* then to calculate the greatest eigenvalue gamma of the in- *
* verse matrix by the power method, and then we have:        *
* lambda=1/gamma (assuming gamma <> 0).                      *
* ---------------------------------------------------------- *
* Ref.: "Algèbre - Algorithmes et programmes en Pascal       *
*        By Jean-Louis Jardrin, Dunod Editor, Paris, 1988"   *
*        [BIBLI 10].                                         * 
* ---------------------------------------------------------- *
* SAMPLE RUN:                                                *
*                                                            *
* ----------------------------------------------             *
*  Calculate the smallest eigenvalue of a real               *
*  square matrix and the associated eigenvector              *
*  by inversion and the power method.                        *
* ----------------------------------------------             *
*                                                            *
*  Size of matrix (maximum 25): 4                            *
*                                                            *
*  Line 1                                                    *
*    Element 1: 1                                            *
*    Element 2: 3                                            *
*    Element 3: 0                                            *
*    Element 4: 0                                            *
*                                                            *
*  Line 2                                                    *
*    Element 1: 4                                            *
*    Element 2: 2                                            *
*    Element 3: 0                                            *
*    Element 4: 0                                            *
*                                                            *
*  Line 3                                                    *
*    Element 1: 1                                            *
*    Element 2: -1                                           *
*    Element 3: 5                                            *
*    Element 4: -3                                           *
*                                                            *
*  Line 4                                                    *
*    Element 1: 2                                            *
*    Element 2: 0                                            *
*    Element 3: 4                                            *
*    Element 4: -2                                           *
*                                                            *
*  Precision: 1e-10                                          *
*  Epsilon  : 1e-10                                          *
*  Maximum number of iterations: 32                          *
*                                                            *
*    DET= -2.00000000000000E+0001                            *
*                                                            *
*    Eigenvalue: 9.99999999638076E-0001                      *
*                                                            *
*    Eigenvector:                                            *
*    -3.09880974276851E-0011                                 *
*     3.09880974276487E-0011                                 *
*     7.49999999925013E-0001                                 *
*     1.00000000000000E+0000                                 *
*                                                            *
*                 English TPW Version By J-P Moreau, Paris.  *
*                            (www.jpmoreau.fr)               *
**************************************************************
  Exact values are: lambda = 1
                    eigenvector = (0,0,0.75,1)
  NOTE: the method fails if the result is an opposite value
        of another eigenvalue or if the matrix A is singular
        (not inversible) or if lambda is not real.
-------------------------------------------------------------}                                  
Program Test_PWIMGT;
Uses WinCrt;

CONST
     NMAX = 25;

TYPE
     MAT = Array[1..NMAX,1..NMAX] of Double;
     VEC = Array[1..NMAX] of Double;
     VECI = Array[1..NMAX] of Integer;

VAR
     i,it,m,n: Integer;
     dta,eps,lambda: Double;
     A: MAT;
     X: VEC;

{***********************************************************
* calculate greatest eigenvalue and associated eigenvector *
* by the power method                                      *
* -------------------------------------------------------- *
* INPUTS:                                                  *
*         eps   : smallest number in double precision      *
*         dta   : required precision                       *
*         m     : maximum number of iterations             *
*         n     : size of real square matrix A(n,n)        *
*         A     : real square matrix A(n,n)                *
* OUTPUTS:                                                 *
*         it    : error indicator: -1=no convergence,      *
*                 0=method cannot be applied,              *
*                 1=convergence ok.                        *
*         gamma : greatest eigenvalue (in absolute value)  *
*                 of input matrix A(n,n)                   *
*         X1    : associated eigenvector                   *
***********************************************************}
Procedure PWM(eps,dta:Double; m,n:Integer; VAR A:MAT; VAR it:Integer;
               VAR gamma:double; VAR X1:VEC);
Var  i,j,l: Integer;
     phi,s: Double;
     X0: VEC;
Begin
  for i:=1 to n do X0[i]:=1.0/SQRT(I);
  it:=-1; l:=1;
  While (it=-1) and (l<=m) do
  begin
    gamma:=0.0;
    for i:=1 to n do
    begin
      X1[i]:=0.0;
      for j:=1 to n do X1[i]:=X1[i]+A[i,j]*X0[j];
      if ABS(X1[i])>ABS(gamma) then gamma:=X1[i]
    end;
    if ABS(gamma) < eps then it:=0
    else
    begin
      for i:=1 to n do X1[i]:=X1[i]/gamma;
      phi:=0.0;
      for i:=1 to n do
      begin
        s:=ABS(X1[i]-X0[i]);
        if s>phi then phi:=s
      end;
      if phi<dta then it:=1
      else
      begin
        X0:=X1;
        Inc(l)
      end
    end
  end
End; {of PWM}                                                                                                     

{*****************************************************************************
*  SOLVING A LINEAR MATRIX SYSTEM AX = B with Gauss-Jordan method using full *
*  pivoting at each step. During the process, original AA and BB matrices are*
* destroyed to spare storage location.                                       *
* -------------------------------------------------------------------------- *
* INPUTS:                                                                    *
*            n    : size of matrix AA                                        *
*            m    : number of simultaneous second members                    *
*            AA   : MATRIX N*N of type MAT                                   *
*            BB   : MATRIX N*M of type MAT (not used here)                   *
* -------------------------------------------------------------------------- *
* OUTPUTS:   AA   : INVERSE OF AA of type MAT                                *
*            DET  : DETERMINANT OF AA (Double)                               *
*            BB   : SOLUTION MATRIX N*M of type MAT (not used here)          *
* -------------------------------------------------------------------------- *
* NOTA:      If M=0 inversion of AA matrix only (BB is not used here).       *
*                                                                            *
*                               TPW english version by J-P Moreau, Paris.    *
*****************************************************************************}
Procedure MATINV(N,M:integer; VAR AA:MAT; VAR BB:MAT;VAR DET:DOUBLE);
CONST EPSMACH=1.2e-16;
VAR   PC,PL,CS : VEC;                                                  
      PV,PAV,temp,TT : DOUBLE;
      I,IK,J,JK,K : integer;

LABEL fin, fin1, fin2;

BEGIN                                                         
{Initializations}                                                             
      DET := 1.0; 
      for I:=1 to N do
      begin                                                                
        PC[I] := 0.0;                                                                
        PL[I] := 0.0;                                                                
        CS[I] := 0.0;
      end;              

{Main loop}                                                                             
      for K:=1 to N do
      begin                                                                  
{Searching greatest pivot}                                               
        PV:=AA[K][K];                                                              
        IK:=K;                                                                    
        JK:=K;                                                                    
        PAV:=abs(PV);                                                            
        for I:=K to N do
          for J:=K to N do
          begin     
            temp := abs(AA[I][J]);                                                        
            if temp > PAV then
            begin                                      
              PV:=AA[I][J];                                                        
              PAV:=abs(PV);
              IK:=I;                                                              
              JK:=J
            end                                                               
          end;                                                                 
                                                                               
{Search terminated, the pivot is in location I=IK, J=JK.                     
 Memorizing pivot location:  }
                                           
        PC[K]:=JK;
        PL[K]:=IK;

{DETERMINANT DET is actualised                                              
 If DET=0, ERROR MESSAGE and STOP                                           
 Machine dependent EPSMACH equals here 1D-20  }                                         
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
            TT:=AA[IK][I];                                                         
            AA[IK][I]:=AA[K][I];
            AA[K][I]:=TT                                                          
          end;                                                                 
                                                                           
        if M<>0 then                                                           
          for I:=1 to M do
          begin                                                               
            TT:=BB[IK][I];                                                           
            BB[IK][I]:=BB[K][I];                                                      
            BB[K][I]:=TT                                                            
          end;                                                                   

{PIVOT is at correct line}
        if JK<>K then                                                         
          for I:=1 to N do
          begin                                                              
{EXCHANGE COLUMNS JK and K of matrix AA:}                                          
            TT:=AA[I][JK];                                                         
            AA[I][JK]:=AA[I][K];                                                    
            AA[I][K]:=TT                                                          
          end;
                                                                           
{The PIVOT is at correct column.                                              
 and is located in K,K.                                                   
                                                                               
 Column K of matrix AA is stored in CS vector                             
 then column K is set to zero. }                                             
        for I:=1 to N do
        begin                                                                
          CS[I]:=AA[I][K];                                                         
          AA[I][K]:= 0.0                                                          
        end;
                                                                               
        CS[K]:= 0.0;                                                                
        AA[K][K]:= 1.0;                                                              
{Line K of matrix AA is modified:}  
        temp := abs(PV);                                          
        if temp < EPSMACH then
        begin
          writeln;                                        
          writeln('  PIVOT TOO SMALL - STOP');
          halt                                                                  
        end;                                                                   
        for I:=1 to N do                                                                
          AA[K][I]:=AA[K][I]/PV;                                                    
                                                                           
        if M<>0 then                                                         
          for I:=1 to M do                                                             
            BB[K][I]:=BB[K][I]/PV;
                                                                           
{Other lines of matrix AA are modified:}                                        
        for J:=1 to N do
        begin                                                                
          if J=K then goto fin;                                                  
          for I:=1 to N do
{Line J of matrix AA is modified:}                                            
            AA[J][I]:=AA[J][I]-CS[J]*AA[K][I];
          if M<>0 then                                                       
            for I:=1 to M do                                                          
              BB[J][I]:=BB[J][I]-CS[J]*BB[K][I];
        fin: end                                                                   
{Line K is ready.}
      end; {of K loop}

{MATRIX AA INVERSION IS DONE - REARRANGEMENT OF MATRIX AA
                                                                               
   EXCHANGE LINES                }                                                            
      for I:=N downto 1 do
      begin                                                               
        IK:=Round(PC[I]);
        if IK=I then goto fin1;                                                   
{EXCHANGE LINES I and PC(I) of matrix AA:}                                         
        for J:=1 to N do
        begin                                                                
          TT:=AA[I][J];                                                            
          AA[I][J]:=AA[IK][J];                                                      
          AA[IK][J]:=TT                                                           
        end;
        if M<>0 then                                                         
          for J:=1 to M do
          begin
            TT:=BB[I][J];                                                          
            BB[I][J]:=BB[IK][J];                                                    
            BB[IK][J]:=TT                                                         
          end;
{NO MORE EXCHANGE is NECESSARY                                                      
 Go to next line.}                                                  
      fin1: end; {of loop i=N downto }                                                                     
                                                                               
{EXCHANGE COLUMNS  }                                                          
      
      for J:=N downto 1 do
      begin                                                                         
        JK:=Round(PL[J]);                                                                
        if JK=J then goto fin2;                                                   
{EXCHANGE COLUMNS J ET PL(J) of matrix AA: }                                       
        for I:=1 to N do
        begin                                                                
          TT:=AA[I][J];                                                            
          AA[I][J]:=AA[I][JK];                                                      
          AA[I][JK]:=TT                                                           
        end; 
{NO MORE EXCHANGE is NECESSARY                                                      
 Go to next column.}
      fin2: end;                                                                     
{REARRANGEMENT TERMINATED.  }                                                        
END;   {of procedure MATINV }


{***********************************************************
* calculate smallest eigenvalue and associated eigenvector *
* by the Gauss-Jordan and the power methods                *
* -------------------------------------------------------- *
* INPUTS:                                                  *
*         eps   : smallest number in double precision      *
*         dta   : required precision                       *
*         m     : maximum number of iterations             *
*         n     : size of real square matrix A(n,n)        *
*         A     : real square matrix A(n,n)                *
* OUTPUTS:                                                 *
*         it    : error indicator: -1=no convergence,      *
*                 0=method cannot be applied,              *
*                 1=convergence ok.                        *
*        lambda : smallest eigenvalue (in absolute value)  *
*                 of input matrix A(n,n)                   *
*         X     : associated eigenvector                   *
***********************************************************}
Procedure PWMIMGT(eps,dta:Double; m,n:Integer; VAR A:MAT; VAR it:Integer;
                    VAR lambda:double; VAR X:VEC);
Var  det,gamma: Double;
     AM1: MAT;  {dummy matrix for MATINV}
Begin
  MATINV(N,0,A,AM1,det); {inverse of A is now in A}
  writeln;
  writeln('   DET=', det);
  if ABS(det)>eps then
  begin
    PWM(eps,dta,m,n,A,it,gamma,X);
    if it=1 then lambda:=1.0/gamma
  end
End; {of PWMIMGT}                                                                                                     

{input data from screen}
Procedure Read_data;
Var i,j: Integer;
Begin
  writeln;
  write(' Size of matrix (maximum ',NMAX,'): '); readln(n);
  for i:=1 to n do
  begin
    writeln;
    writeln(' Line ',i);
    for j:=1 to n do
    begin
      write('  Element ',j,': '); readln(A[i,j])
    end
  end;
  writeln;
  write(' Precision: '); readln(dta);
  write(' Epsilon  : '); readln(eps);
  write(' Maximum number of iterations: '); readln(m)
End;

{main program}
BEGIN
  writeln(' ----------------------------------------------');
  writeln('  Calculate the smallest eigenvalue of a real  ');
  writeln('  square matrix and the associated eigenvector ');
  writeln('  by inversion and the power method.           ');
  writeln(' ----------------------------------------------');

  Read_Data;

  PWMIMGT(eps,dta,m,n,A,it,lambda,X);

  Case it+1 of
    0: writeln('  No convergence !');
    1: writeln('  Method does not apply.');
    2: begin
         writeln;
         writeln('  Eigenvalue: ', lambda);
         writeln;
         writeln('  Eigenvector:');
         for i:=1 to n do writeln('  ',X[i])
       end
  End;
  writeln;
  Readkey; DoneWinCrt
END.

{end of file tpwimgt.pas}