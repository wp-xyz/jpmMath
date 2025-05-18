{****************************************************************
* This program seeks the equation of a conical passing through  *
* 5 given points.                                               *
* ------------------------------------------------------------- *
* SAMPLE RUN:                                                   *
*                                                               *
*    SEEK A CONICAL PASSING THROUGH 5 POINTS                    *
*                                                               *
*        X1  Y1 = 0 2                                           *
*        X2  Y2 = 1 5                                           *
*        X3  Y3 = 5 77                                          *
*        X4  Y4 = -4 50                                         *
*        X5  Y5 = 2 14                                          *
*                                                               *
*                                                               *
*  The coefficients of the conical ax2+bxy+cy2+dx+ey+f = 0 are: *
*                                                               *
*   a =   1.500000                                              *
*   b =   0.000000                                              *
*   c =  -0.000000                                              *
*   d =  -0.000000                                              *
*   e =  -0.500000                                              *
*   f =   1.000000                                              *
*                                                               *
*                           TPW version by J-P Moreau, Paris.   *
*                                  (www.jpmoreau.fr)            *
*****************************************************************
 Explanations:
 ------------
 We seek the cartesian equation ax^2+bxy+cy^2+dx+ey+f=0 of the
 conical passing through 5 points: (x1,y1),(x2,y2)... to (x5,y5).

 The unknowns are a,b,c,d,e,f. Since the equation can be multiplied
 by any real number, the actual number of independant unknowns is 5.

 That is the reason why only five different points are necassary to
 determine the conical. The linear system to solve is a system of
 five homogeneous equations with 6 unknowns.

 The resolution of the system will be made by successive steps:

    Let us suppose that f<>0 and let us fix f=1, then the other
    parameters a,b,c,d,e are solutions of the linear system:

    x1^2 a + x1y1 b + y1^2 c + x1 d + y1 e = -1
    x2^2 a + x2y2 b + y2^2 c + x2 d + y2 e = -1
    -------------------------------------------
    x5^2 a + x5y5 b + y5^2 c + x5 d + y5 e = -1

    This linear system of size 5 is solved by a classical Gauss-Jordan
    method found in procedure SolveSystem (see also program MATRICES/sysmat.pas).

    if the system has a solution (det <> 0) the problem is finished, else
    the hypothesis f<>0 was wrong and we then try e<>0 (that is e=1).

    The linear system to solve is then:

    x1^2 a + x1y1 b + y1^2 c + x1 d + f = -y1
    x2^2 a + x2y2 b + y2^2 c + x2 d + f = -y2
    -------------------------------------------
    x5^2 a + x5y5 b + y5^2 c + x5 d + f = -y5

    if the system has a solution (det <> 0) the problem is finished, else
    the hypothesis e<>0 was again wrong and we successively explore d<>0,
    c<>0... to a<>0.

    So in the most unfavourable case, that is if only a degenerated conical
    answers the problem, we can have up to six linear systems of size 5
    to solve.

    Usually (f<>0 or e<>0), we only solve one or two linear systems.  

----------------------------------------------------------------------------------}
PROGRAM Seek_a_conical;
Uses WinCrt;


Type  MAT5 = array[1..5,1..5] of DOUBLE;
      VEC5 = array[1..5,1..1] of DOUBLE;

Var   x1,y1,x2,y2,x3,y3,x4,y4,x5,y5:DOUBLE;
      a,b,c,d,e,f:DOUBLE;
      det: DOUBLE;


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
*                                         *
*           Pascal version by J-P Moreau  *
*          (here limited to N=5 and M=1). *
******************************************}
Procedure SolveSystem(N,M:integer; VAR AA:MAT5; VAR BB:VEC5;VAR DET:DOUBLE);
LABEL 50;
CONST EPSMACH=2E-16;
VAR   PC,PL,CS : array[1..5] of DOUBLE;
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
        if temp < EPSMACH then goto 50;

                                                                               
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
        if temp < EPSMACH then goto 50;
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
50:END;   {of procedure SOLVESYSTEM}



Function SeekConical(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5:DOUBLE;VAR a,b,c,d,e,f:DOUBLE):Boolean;
Label 100;
Var
    tm      : MAT5;
    tv      : VEC5;
    x,y     : Array[1..5] of DOUBLE;
    i,j,iter: INTEGER;

    {internal procedure}
    Procedure Solve(i:INTEGER);
    begin
      SolveSystem(5,1,tm,tv,det);
      if det<>0 then
      begin
        case i of
          1:begin {case f<>0}
              a:=tv[1,1]; b:=tv[2,1]; c:=tv[3,1]; d:=tv[4,1]; e:=tv[5,1]; f:=1
            end;
          2:begin {case e<>0}
              a:=tv[1,1]; b:=tv[2,1]; c:=tv[3,1]; d:=tv[4,1]; e:=1; f:=tv[5,1]
            end;
          3:begin {case d<>0}
              a:=tv[1,1]; b:=tv[2,1]; c:=tv[3,1]; d:=1; e:=tv[4,1]; f:=tv[5,1]
            end;
          4:begin {case c<>0}
              a:=tv[1,1]; b:=tv[2,1]; c:=1; d:=tv[3,1]; e:=tv[4,1]; f:=tv[5,1]
            end;
          5:begin {case b<>0}
              a:=tv[1,1]; b:=1; c:=tv[2,1]; d:=tv[3,1]; e:=tv[4,1]; f:=tv[5,1]
            end;
          6:begin {case a<>0}
              a:=1; b:=tv[1,1]; c:=tv[2,1]; d:=tv[3,1]; e:=tv[4,1]; f:=tv[5,1]
            end
        end;
        SeekConical:=TRUE
      end
    end;

Begin  {SeekConical}
  SeekConical:=FALSE;
  x[1]:=x1; y[1]:=y1;
  x[2]:=x2; y[2]:=y2;
  x[3]:=x3; y[3]:=y3;
  x[4]:=x4; y[4]:=y4;
  x[5]:=x5; y[5]:=y5;
  {verify that the 5 given points are different}
  for i:=2 to 5 do
    for j:=1 to Pred(i) do
      if (x[i]=x[j]) and (y[i]=y[j]) then exit;
  {try f<>0}
  for i:=1 to 5 do
  begin
    tm[i,1]:=x[i]*x[i]; tm[i,2]:=x[i]*y[i]; tm[i,3]:=y[i]*y[i];
    tm[i,4]:=x[i]; tm[i,5]:=y[i]; tv[i,1]:=-1
  end;
  Solve(1);
  if det<>0 then goto 100;
  {try e<>0}
  for i:=1 to 5 do
  begin
    tm[i,1]:=x[i]*x[i]; tm[i,2]:=x[i]*y[i]; tm[i,3]:=y[i]*y[i];
    tm[i,4]:=x[i]; tm[i,5]:=1; tv[i,1]:=-y[i]
  end;
  Solve(2);
  if det<>0 then goto 100;
  {try d<>0}
  for i:=1 to 5 do
  begin
    tm[i,1]:=x[i]*x[i]; tm[i,2]:=x[i]*y[i]; tm[i,3]:=y[i]*y[i];
    tm[i,4]:=y[i]; tm[i,5]:=1; tv[i,1]:=-x[i]
  end;
  Solve(3);
  if det<>0 then goto 100;
  {try c<>0}
  for i:=1 to 5 do
  begin
    tm[i,1]:=x[i]*x[i]; tm[i,2]:=x[i]*y[i]; tm[i,3]:=x[i];
    tm[i,4]:=y[i]; tm[i,5]:=1; tv[i,1]:=-y[i]*y[i]
  end;
  Solve(4);
  if det<>0 then goto 100;
  {try b<>0}
  for i:=1 to 5 do
  begin
    tm[i,1]:=x[i]*x[i]; tm[i,2]:=y[i]*y[i]; tm[i,3]:=x[i];
    tm[i,4]:=y[i]; tm[i,5]:=1; tv[i,1]:=-x[i]*y[i]
  end;
  Solve(5);
  if det<>0 then goto 100;
  {last try a<>0}
  for i:=1 to 5 do
  begin
    tm[i,1]:=x[i]*y[i]; tm[i,2]:=y[i]*y[i]; tm[i,3]:=x[i];
    tm[i,4]:=y[i]; tm[i,5]:=1; tv[i,1]:=-x[i]*x[i]
  end;
  Solve(6);
100: End;      {det<>0 means success}


{main program}
BEGIN
  writeln;
  writeln('    SEEK A CONICAL PASSING THROUGH 5 POINTS');
  writeln;
  write('        X1  Y1 = ');  read(x1,y1);
  write('        X2  Y2 = ');  read(x2,y2);
  write('        X3  Y3 = ');  read(x3,y3);
  write('        X4  Y4 = ');  read(x4,y4);
  write('        X5  Y5 = ');  read(x5,y5);
  writeln;
  writeln;
  writeln;
  if SeekConical(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,a,b,c,d,e,f) then
  begin
    writeln('  The coefficients of the conical ax2+bxy+cy2+dx+ey+f = 0 are:');
    writeln;
    writeln('   a=',a:10:6);
    writeln('   b=',b:10:6);
    writeln('   c=',c:10:6);
    writeln('   d=',d:10:6);
    writeln('   e=',e:10:6);
    writeln('   f=',f:10:6)
  end
  else
    writeln('  No conical found.');
  writeln;

  ReadKey; DoneWinCrt

END.

{end of file conical1.pas}