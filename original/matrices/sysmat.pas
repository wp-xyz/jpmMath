{********************************************************
*          SOLVING A LINEAR MATRIX SYSTEM AX=B          *
*                by Gauss-Jordan method                 *
* ----------------------------------------------------- *
*                 Pascal version by J-P Moreau, Paris   *
*                         (From [BIBLI 13]).            *
*                         (www.jpmoreau.fr)             *
*                                                       *
* INPUTS  from example file sysmat.dat,                 *
* OUTPUTS to sysmat.lst.                                *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
* Input file                                            *
* 3                                                     *
* 2.0  -1.0   1.0                                       *
* 1.0   5.0  -2.0                                       *
* 3.0  -2.0   3.0                                       *
* 4                                                     *
* 5.0  -1.0   1.0  2.0                                  *
* 1.0   2.0   3.0  4.0                                  *
* 3.0   7.556 4.0  4.0                                  *
*                                                       *
* Output file                                           *
* ------------------------------------------------------*
*  SOLVING THE MATRIX LINEAR SYSTEM AX = B              *
* ------------------------------------------------------*
*     (BY GAUSS-JORDAN METHOD)                          *
*                                                       *
*  N=3                                                  *
*                                                       *
*  MATRIX A:                                            *
*                                                       *
*   2.000000  -1.000000   1.000000                      *
*   1.000000   5.000000  -2.000000                      *
*   3.000000  -2.000000   3.000000                      *
*                                                       *
*  M=4                                                  *
*                                                       *
*  MATRIX B:                                            *
*                                                       *
*   5.000000  -1.000000   1.000000   2.000000           *
*   1.000000   2.000000   3.000000   4.000000           *
*   3.000000   7.556000   4.000000   4.000000           *
*                                                       *
*  INVERSE OF MATRIX A:                                 *
*                                                       *
*   0.785714   0.071429  -0.214286                      *
*  -0.642857   0.214286   0.357143                      *
*  -1.214286   0.071429   0.785714                      *
*                                                       *
*  SOLUTION MATRIX X:                                   *
*                                                       *
*   3.357143  -2.262000   0.142857   1.000000           *
*  -1.928571   3.770000   1.428571   1.000000           *
*  -3.642857   7.294000   2.142857   1.000000           *
*                                                       *
*  DETERMINANT= 14.000000                               *
*                                                       *
*                                                       *
*  VERIFICATION OF AX = B:                              *
*                                                       *
*   5.000000  -1.000000   1.000000   2.000000           *
*   1.000000   2.000000   3.000000   4.000000           *
*   3.000000   7.556000   4.000000   4.000000           *
*                                                       *
*                                                       *
* ------------------------------------------------------*
*                                                       *
********************************************************}
Uses WinCrt;

CONST

  EPSMACH = 2E-16;

  NMAX = 10;
  MMAX = 5;

TYPE

  MAT = Array[1..NMAX,1..NMAX] of DOUBLE;
  MAT1= Array[1..NMAX,1..MMAX] of DOUBLE;

  VECT= ARRAY[1..NMAX] of DOUBLE;
  
             
VAR
        
  fp1, fp2 : TEXT;
  A,A1,D : MAT;                                             
  B,C : MAT1;
  DETER, temp : DOUBLE;
  I,J,M,N : integer;                                             

  {write empty line to output file}
  Procedure ligne;
  begin
    writeln(fp2)
  end;


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
Procedure MATINV(N,M:integer; VAR AA:MAT; VAR BB:MAT1;VAR DET:DOUBLE);
VAR   PC,PL,CS : VECT;                                                  
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
          ligne;                                                  
          writeln(fp2,'  The determinant equals ZERO !!!');                                                              
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
          ligne;                                        
          writeln(fp2,'  PIVOT TOO SMALL - STOP');                               
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



{******************************************                                     
*    MULTIPLICATION OF TWO SQUARE REAL    *                                     
*    MATRICES                             *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*N                *                                     
*            B  MATRIX N*N                *                                     
*            N  INTEGER                   *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*N PRODUCT A*B    *                                     
*                                         *
******************************************}
Procedure MATMUL(A, B:MAT; VAR C: MAT; N: integer);
VAR
      SUM : DOUBLE;
      I,J,K : integer;
BEGIN                                               
      for I:=1 to N do                                                                  
        for J:=1 to N do
        begin                                                                
          SUM:= 0.0;                                                                
          for K:=1 to N do                                                              
            SUM:=SUM+A[I][K]*B[K][J];                                               
          C[I][J]:=SUM                                                            
        end                                                                   
END;                                                                       

Procedure MATMUL1(A:MAT; B:MAT1; VAR C:MAT1; N,M: integer);
{******************************************                                     
*   MULTIPLICATION OF TWO REAL MATRICES   *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*N                *                                     
*            B  MATRIX N*M                *                                     
*            N  INTEGER                   *
*            M  INTEGER                   *                                     
* --------------------------------------- *
* OUTPUTS:   C  MATRIX N*M PRODUCT A*B    *                                     
*                                         *                                     
******************************************}
VAR 
      SUM : DOUBLE;  
      I,J,K : integer;

BEGIN                                              
      for I:=1 to N do                                                                  
        for J:=1 to M do
        begin                                                                
          SUM:= 0.0;                                                                
          for K:=1 to N do                                                              
            SUM:=SUM+A[I][K]*B[K][J];                                               
          C[I][J]:=SUM
        end                                                                   
END;


Procedure FWriteHead(VAR f:TEXT; nom:PChar);
Begin
  Writeln(f,'------------------------------------------------------');
  Writeln(f,nom);
  Writeln(f,'------------------------------------------------------')
End;

Procedure FWriteEnd(VAR f:TEXT; nom:PChar);
Begin
  Writeln(f,'------------------------------------------------------')
End;


{***** main programme *****}
BEGIN

{  READ LEFT HAND MATRIX FROM INPUT FILE AND ECHO TO OUTPUT FILE}

      Assign(fp1,'sysmat.dat'); reset(fp1);
      Assign(fp2,'sysmat.lst'); rewrite(fp2);

      read(fp1,N);

      FWriteHead(fp2,'  SOLVING THE MATRIX LINEAR SYSTEM AX = B');
      writeln(fp2,'     (BY GAUSS-JORDAN METHOD)');
      ligne;
      writeln(fp2,'  N=',N);
      ligne;

      for I:=1 to N do
        for J:=1 to N do
        begin
          read(fp1,temp);
          A[I][J] := temp
        end;

      writeln(fp2,'  MATRIX A:');
      ligne;

      for I:=1 to N do
        for J:=1 to N do
        begin
          write(fp2,'  ',A[I][J]:9:6);
          if J=N then ligne
        end;
                                                                               
{READ FROM INPUT FILE RIGHT HAND MATRIX IF M>0}                                              

      read(fp1,M);

      ligne;
      writeln(fp2,'  M=',M);
      ligne;
                                                                             
      if M<>0 then
      begin                                                           
        for I:=1 to N do
          for J:=1 to M do
          begin
            read(fp1,temp); 
            B[I][J] := temp
          end;     

        writeln(fp2,'  MATRIX B:');
        ligne;

        for I:=1 to N do
          for J:=1 to M do
          begin
            write(fp2,'  ',B[I][J]:9:6);
            if J=M then ligne
          end
      end;
{END SECTION READ DATA}      
      close(fp1);                                                                     
                                                                               
{STORE MATRIX A IN MATRIX A1}

      for I:=1 to N do
        for J:=1 to N do
          A1[I][J]:=A[I][J];        
                                                                               
{CALL MATRIX INVERSION ROUTINE}                                               
                                                                               
      MATINV(N,M,A,B,DETER);
                                                                               
{PRINT RESULTS TO OUTPUT FILE}                                                          

      ligne;                                                                         
      writeln(fp2,'  INVERSE OF MATRIX A:');
      ligne;   
      for I:=1 to N do
        for J:=1 to N do
        begin
          write(fp2,'  ',A[I][J]:9:6);   
          if J=N then ligne
        end;

      if M<>0 then
      begin
        ligne;
        writeln(fp2,'  SOLUTION MATRIX X:');
        ligne;
        for I:=1 to N do
          for J:=1 to M do
          begin
            write(fp2,'  ',B[I][J]:9:6);
            if J=M then ligne
          end
      end;

      ligne;
      writeln(fp2,'  DETERMINANT= ',DETER:9:6);
      ligne;

{VERIFY THAT NEW PRODUCT A*B EQUALS ORIGINAL B MATRIX
   (OPTIONAL)     }

      if M<>0 then
      begin
        ligne;                                                       
        writeln(fp2,'  VERIFICATION OF AX = B:');
        ligne;                                                            
{CALL MULTIPLICATION ROUTINE}
        MATMUL1(A1,B,C,N,M);
        for I:=1 to N do
          for J:=1 to M do
          begin
            write(fp2,'  ',C[I][J]:9:6);   
            if J=M then ligne
          end
      end                                                                 
      else
      begin
        ligne;                                                                    
        writeln(fp2,'  VERIFICATION OF A1 * A = I:');
        ligne;
{CALL MULTIPLICATION ROUTINE}                                         
        MATMUL(A1,A,D,N);   
        for I:=1 to N do
          for J:=1 to N do
          begin
            write(fp2,'  ',D[I][J]:9:6);   
            if J=N then ligne
          end  
      end;                                                         

      ligne;
      ligne;
      FWriteEnd(fp2,'  End of results file.');
      close(fp2);
      writeln;
      writeln('  Results in file sysmat.lst...');
      readln;
      DoneWinCrt
END.                                                                          

{end of file sysmat.pas j-p moreau july 1997}