{****************************************************
*     Solving a complex linear system A Z = B,      *
*     where A is complex matrix of size N*N,        *
*     B is a complex matrix of size N*M,            *
*     Z is a complex solution matrix of size N*M.   *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal by M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* Solve complex linear system:                      *
*                                                   *
*   z1 +  z2 + iz3 = 5                              *
*   z1 + iz2 +  z3 = -2i                            *
*  iz1 +  z2 +  z3 = 1                              *
*                                                   *
*  Here: N=3, M=1                                   *
*                                                   *
*  System solution:                                 *
*                                                   *
*  z1 =   1.500000  i  -0.500000                    *
*  z2 =   1.000000  i   1.000000                    *
*  z3 =  -0.500000  i  -2.500000                    *
*                                                   *
*  DET =  2.00000000000000E+0001                    *
*                                                   *
*****************************************************
 Explanations:
 The original complex linear system to solve is:
     [ a11 a12 ... a1N ] [ z1 ]   [ b1 ]
     [ a21 a22 ... a2N ] [ z2 ]   [ b2 ]
     [ ... ... ... ... ] [... ] = [... ]
     [ aN1 aN2 ... aNN ] [ zN ]   [ bN ]
     a k,l = c k,l + i d k,l
     z k = x k + i y k
     b k = e k + i f k
     where c,d,e,f,x,y are real numbers.
 The system is replaced by the following REAL system of size 2*N:
     [ c11 -d11 c12 -d12 ... c1N -d1N ] [ x1 ]   [ e1 ]
     [ d11  c11 d12  c12 ... d1N  c1N ] [ y1 ]   [ f1 ]
     [ c21 -d21 c22 -d22 ... c2N -d2N ] [ x2 ]   [ e2 ]
     [ d21  c21 d22  c22 ... d2N  c2N ] [ y2 ] = [ f2 ]
     [ ... .... ... .... ... ... .... ] [ ...]   [ ...]
     [ cN1 -dN1 cN2 -dN2 ... cNN -dNN ] [ xN ]   [ eN ]
     [ dN1  cN1 dN2  cN2 ... dNN  cNN ] [ yN ]   [ fN ]
 The real system is then solved by a classic Gauss-Jordan method.         
}
PROGRAM Zmatsys;
Uses Wincrt;

Const
       EPSMACH = 2E-16;
       NMAX = 3;
       MMAX = 1;

Type
       MAT  = Array[1..NMAX,1..NMAX] of DOUBLE;
       VEC  = Array[1..NMAX] of DOUBLE;
       VEC2 = Array[1..2*NMAX] of DOUBLE;
       MAT2 = Array[1..2*NMAX,1..2*NMAX] of DOUBLE;
       MAT3 = Array[1..2*NMAX,1..MMAX] of DOUBLE;

Var
       Xm, Ym : MAT;
       Xv, Yv : VEC;
       M      : MAT2;
       V      : MAT3;

       i,j,li,co : integer;
       det       : DOUBLE;


{******************************************
*  SOLVING A LINEAR MATRIX SYSTEM AX = B  *
*  with Gauss-Jordan method using full    *
*  pivoting at each step. During the pro- *
* cess, original AA and BB matrices are   *
* destroyed to spare storage location.    *
* --------------------------------------- *
* INPUTS:    AA   MATRIX N*N              *                                     
*            BB   MATRIX N*M              *                                     
* --------------------------------------- *                                     
* OUTPUS:    AA   INVERSE OF AA N*N       *                                     
*            DET  DETERMINANT OF AA       *                                     
*            BB   SOLUTION MATRIX N*M     *                                     
* --------------------------------------- *                                     
* NOTA - If M=0 inversion of AA matrix    *
*        only.                            *                                     
******************************************}
Procedure MATINV(N,M:integer; VAR AA:MAT2; VAR BB:MAT3; VAR DET:DOUBLE);
VAR   PC,PL,CS : VEC2;
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
 Machine dependant EPSMACH equals here 1D-20  }                                         
        if IK<>K then DET:=-DET;                                                   
        if JK<>K then DET:=-DET;                                                   
        DET:=DET*PV;  
        temp:= abs(DET);
        if temp < EPSMACH then
        begin                                          
          writeln;
          writeln(' The determinant equals ZERO !!!');
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
          writeln(' PIVOT TOO SMALL - STOP');
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


BEGIN
  li:=NMAX; co:=MMAX;
  {read Xm matrix}                        {real part of A matrix: }
  Xm[1,1]:= 1; Xm[1,2]:=1; Xm[1,3]:=0;    {  1  1  0              }
  Xm[2,1]:= 1; Xm[2,2]:=0; Xm[2,3]:=1;    {  1  0  1              }
  Xm[3,1]:= 0; Xm[3,2]:=1; Xm[3,3]:=1;    {  0  1  1              } 
  {read Ym matrix}                        {imaginary part of A matrix}
  Ym[1,1]:= 0; Ym[1,2]:=0; Ym[1,3]:=1;    {  0  0  1                 }
  Ym[2,1]:= 0; Ym[2,2]:=1; Ym[2,3]:=0;    {  0  1  0                 }
  Ym[3,1]:= 1; Ym[3,2]:=0; Ym[3,3]:=0;    {  1  0  0                 }
  {read Xv vector}                        {real part of 1st column of B matrix}
  Xv[1]:=5; Xv[2]:=0; Xv[3]:=1;           {  5  0  1                          }
  {read Yv vector}                        {imaginary part of 1st column of B matrix}
  Yv[1]:=0; Yv[2]:=-2; Yv[3]:=0;          {  0 -2  0                               }

  {The complex system AZ=B is replaced by real system MX=V of size 2*NMAX}
  {create M matrix}
  for i:=1 to NMAX do
    for j:=1 to NMAX do
    begin
      M[2*i-1,2*j-1]:=Xm[i,j]; M[2*i,2*j]:=Xm[i,j];
      M[2*i,2*j-1]:=Ym[i,j]; M[2*i-1,2*j]:=-Ym[i,j]
    end;
  {create V matrix}
  for i:=1 to NMAX do
  begin
    V[2*i-1,1]:=Xv[i]; V[2*i,1]:=Yv[i];
  end;

  {solve real linear system MX=V of zize 2*NMAX
   by Gauss-Jordan method - See program sysmat.pas}
  MATINV(li*2, co, M, V, det);

  {write complex solutions stored in vector V in the form:
   x1,y1,x2,y2,...,xn,yn, with zk = xk + i yk... }
  writeln;
  writeln('  System solution:');
  writeln;
  for i:=1 to li do
    writeln('  z',i,' = ',V[2*i-1,1]:10:6,'  i ',V[2*i,1]:10:6);
  writeln;
  writeln('  DET = ',det);
  writeln;
  ReadKey; DoneWinCrt

END.

{end of file zmatsys.pas}