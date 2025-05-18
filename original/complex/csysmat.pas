{**********************************************************************************************
*            SOLVING A COMPLEX LINEAR MATRIX SYSTEM AX=B BY GAUSS-JORDAN METHOD               *
* ------------------------------------------------------------------------------------------- *
*                                                       Pascal version by J-P Moreau, Paris   *
*                                                                (www.jpmoreau.fr)            *
* SAMPLE RUN:                                                                                 *
* INPUTS  from example file csysmat.dat,                                                      *
* OUTPUTS to file csysmat.lst.                                                                *
* ------------------------------------------------------------------------------------------- *
* SAMPLE RUN:                                                                                 *
* Input file                                                                                  *
* 4                                                                                           *
* 47.0 -15.0  62.0 5.0    0.0 -72.0 61.0  20.0                                                *
*  6.0  14.0 -17.0 3.0 -102.0  91.0  7.0 -12.0                                                *
* 13.0 -55.0  32.0 8.0   41.0   7.0 25.0   1.0                                                *
*111.0  25.0  40.0 0.0   12.0 -82.0 58.0 -30.0                                                *
* 1                                                                                           *
* 629.0 988.0                                                                                 *
*-180.0 825.0                                                                                 *
* 877.0 441.0                                                                                 *
* 734.0 -88.0                                                                                 *
*                                                                                             *
* Output file                                                                                 *
*------------------------------------------------------                                       *
*  SOLVING THE COMPLEX MATRIX LINEAR SYSTEM AX = B                                            *
*------------------------------------------------------                                       *
*    (BY GAUSS-JORDAN METHOD WITH FULL PIVOTING)                                              *
*                                                                                             *
*  N=4                                                                                        *
*                                                                                             *
*  COMPLEX MATRIX A:                                                                          *
*                                                                                             *
* (47.000000,-15.000000) (62.000000, 5.000000) ( 0.000000,-72.000000) (61.000000,20.000000)   *
* ( 6.000000,14.000000) (-17.000000, 3.000000) (-102.000000,91.000000) ( 7.000000,-12.000000) *
* (13.000000,-55.000000) (32.000000, 8.000000) (41.000000, 7.000000) (25.000000, 1.000000)    *
* (111.000000,25.000000) (40.000000, 0.000000) (12.000000,-82.000000) (58.000000,-30.000000)  *
*                                                                                             *
*  M=1                                                                                        *
*                                                                                             *
*  COMPLEX MATRIX B:                                                                          *
*                                                                                             *
* (629.000000,988.000000)                                                                     *
* (-180.000000,825.000000)                                                                    *
* (877.000000,441.000000)                                                                     *
* (734.000000,-88.000000)                                                                     *
*                                                                                             *
*  INVERSE OF COMPLEX MATRIX A:                                                               *
*                                                                                             *
* ( 0.003941,-0.007988) (-0.008374,-0.004058) (-0.009490, 0.017938) ( 0.000210,-0.001660)     *
* ( 0.036781, 0.021820) ( 0.014793,-0.027076) (-0.033797,-0.035129) (-0.005015,-0.016173)     *
* (-0.004225,-0.005071) (-0.007674,-0.001821) ( 0.003482, 0.004805) ( 0.000692, 0.003649)     *
* (-0.019727,-0.016588) (-0.001490, 0.018806) ( 0.033739, 0.015368) ( 0.005363, 0.017239)     *
*                                                                                             *
*  SOLUTION MATRIX X:                                                                         *
*                                                                                             *
* (-1.000000, 3.000000)                                                                       *
* ( 2.000000,10.000000)                                                                       *
* ( 7.000000,-5.000000)                                                                       *
* (17.000000, 6.000000)                                                                       *
*                                                                                             *
*  DETERMINANT= (17416564.000000,-10598320.000000)                                            *
*                                                                                             *
*                                                                                             *
*  VERIFICATION OF A1 * A = I:                                                                *
*                                                                                             *
* ( 1.000000, 0.000000) (-0.000000, 0.000000) ( 0.000000, 0.000000) ( 0.000000, 0.000000)     *
* (-0.000000, 0.000000) ( 1.000000, 0.000000) ( 0.000000,-0.000000) ( 0.000000,-0.000000)     *
* ( 0.000000, 0.000000) ( 0.000000,-0.000000) ( 1.000000,-0.000000) ( 0.000000,-0.000000)     *
* (-0.000000, 0.000000) ( 0.000000, 0.000000) ( 0.000000,-0.000000) ( 1.000000, 0.000000)     *
*                                                                                             *
*  VERIFICATION OF AX = B:                                                                    *
*                                                                                             *
* (629.000000,988.000000)                                                                     *
* (-180.000000,825.000000)                                                                    *
* (877.000000,441.000000)                                                                     *
* (734.000000,-88.000000)                                                                     *
*                                                                                             *
* --------------------------------------------------------------------------------------------*
*                                                                                             *
**********************************************************************************************}
PROGRAM CSYSMAT;
Uses WinCrt;

CONST

  EPSMACH = 2.2E-16;

  NMAX = 15;
  MMAX = 5;

Type
     Complex = Record
       R, I: Real
     End;

     Matc = Array[1..NMAX,1..NMAX] of Complex;
     Matc1= Array[1..NMAX,1..MMAX] of Complex;
     Vecc = Array[1..NMAX] of Complex;
     Veci = Array[1..NMAX] of Integer;

VAR
        
  fp1, fp2 : TEXT;
  A,A1,D : MATC;                                             
  B,C : MATC1;
  DETER, temp : Complex;
  I,J,M,N : integer;                                             


  Function CABS(Z:Complex): Real;
  Begin
    CABS := sqrt(Z.R*Z.R+Z.I*Z.I)
  End;

  Procedure CADD(Z1,Z2:Complex; Var Z:Complex);
  Begin
    Z.R := Z1.R + Z2.R;
    Z.I := Z1.I + Z2.I
  End;

  Procedure CDIF(Z1,Z2:Complex; Var Z:Complex);
  Begin
    Z.R := Z1.R - Z2.R;
    Z.I := Z1.I - Z2.I
  End;

  Procedure CMUL(Z1,Z2:Complex; Var Z:Complex);
  Begin
    Z.R := Z1.R*Z2.R - Z1.I*Z2.I;
    Z.I := Z1.R*Z2.I + Z1.I*Z2.R
  End;

  Procedure CDIV(Z1,Z2:Complex; Var Z:Complex);
  Var d:Real; C:Complex;
  Begin
    d := Z2.R*Z2.R+Z2.I*Z2.I;
    if d<1E-10 then
      writeln(' Complex Divide by zero !')
    else
    begin
      C.R:=Z2.R; C.I:=-Z2.I;
      CMUL(Z1,C,Z);
      Z.R:=Z.R/d; Z.I:=Z.I/d
    end
  End;

  Procedure CSwap(Var Z1,Z2:Complex);
  Var C: Complex;
  Begin
    C:=Z1; Z1:=Z2; Z2:=C
  End;   

  {write empty line to output file}
  Procedure ligne;
  begin
    writeln(fp2)
  end;


{***********************************************                                     
* SOLVING A COMPLEX MATRIX LINEAR SYSTEM AX=B  *
* with Gauss-Jordan method using full pivoting *
* at each step. During the process, original   *
* AA and BB matrices are destroyed to spare    *
* storage location.                            *
* -------------------------------------------- *
* INPUTS:    AA   COMPLEC MATRIX N*N           *                                     
*            BB   COMPLEX MATRIX N*M           *                                     
* -------------------------------------------- *                                     
* OUTPUTS:   AA   INVERSE OF AA N*N            *                                     
*            DET  COMPLEX DETERMINANT OF AA    *                                     
*            BB   SOLUTION COMPLEX MATRIX N*M  *
* -------------------------------------------- *                                     
* NOTA - If M=0 inversion of AA matrix only.   *                                     
***********************************************}
Procedure CMATINV(N,M:integer; VAR AA:MATC; VAR BB:MATC1;VAR DET:Complex);
VAR   PC,PL: VECI; CS:Vecc;                                                  
      CZERO,PV,temp: Complex;
      I,IK,J,JK,K : integer;
      PAV,tmp:Real;

LABEL fin, fin1, fin2;

BEGIN                                                         
{Initializations}
      CZERO.R:=0.0; CZERO.I:=0.0;
      DET.R := 1.0; DET.I:=0.0; 
      for I:=1 to N do
      begin                                                                
        PC[I] := 0;                                                                
        PL[I] := 0;                                                                
        CS[I] := CZERO
      end;              

{Main loop}                                                                             
      for K:=1 to N do
      begin                                                                  
{Searching greatest pivot}                                               
        PV:=AA[K,K];                                                              
        IK:=K;                                                                    
        JK:=K;                                                                    
        PAV:=CABS(PV);                                                            
        for I:=K to N do
          for J:=K to N do
          begin     
            tmp := CABS(AA[I,J]);
            if tmp > PAV then
            begin                                      
              PV:=AA[I,J];                                                        
              PAV:=CABS(PV);
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
 Machine dependant EPSMACH equals here 2.2D-16 }                                           
        if IK<>K then begin DET.R:=-DET.R; DET.I:=-DET.I end;
        if JK<>K then begin DET.R:=-DET.R; DET.I:=-DET.I end;
        CMUL(DET,PV,DET);  
        tmp:= CABS(DET);
        if tmp < EPSMACH then
        begin                                          
          ligne;                                                  
          writeln(fp2,'  The complex determinant equals ZERO !');                                                              
          halt                                                                  
        end;

                                                                               
{POSITIONNING PIVOT IN K,K:}
        if IK<>K then
{EXCHANGE LINES IK and K of matrix AA:}
          for I:=1 to N do CSwap(AA[K,I],AA[IK,I]);
                                                                           
        if M<>0 then                                                           
          for I:=1 to M do CSwap(BB[K,I],BB[IK,I]);

{PIVOT is at correct line}
        if JK<>K then
{EXCHANGE COLUMNS JK and K of matrix AA:}
          for I:=1 to N do CSwap(AA[I,K],AA[I,JK]);
                                                                           
{The PIVOT is at correct column, and is located in K,K.
                                                                               
 Column K of matrix AA is stored in CS vector                             
 then column K is set to zero. }                                             
        for I:=1 to N do
        begin                                                                
          CS[I]:=AA[I,K];                                                         
          AA[I,K]:= CZERO                                                          
        end;
                                                                               
        CS[K]:= CZERO;                                                                
        AA[K,K].R:= 1.0; AA[K,K].I:= 0.0;
{Line K of matrix AA is modified:}  
        tmp := CABS(PV);                                          
        if tmp < EPSMACH then
        begin
          ligne;                                        
          writeln(fp2,'  COMPLEX PIVOT TOO SMALL - STOP');                               
          halt                                                                  
        end;                                                                   
        for I:=1 to N do                                                                
          CDIV(AA[K,I],PV,AA[K,I]);                                                    
                                                                           
        if M<>0 then                                                         
          for I:=1 to M do                                                             
            CDIV(BB[K,I],PV,BB[K,I]);
                                                                           
{Other lines of matrix AA are modified:}                                        
        for J:=1 to N do
        begin                                                                
          if J=K then goto fin;                                                  
          for I:=1 to N do
          begin
{Line J of matrix AA is modified:                                            
            AA[J,I]:=AA[J,I]-CS[J]*AA[K,I]; }
            CMUL(CS[J],AA[K,I],temp);
            CDIF(AA[J,I],temp,AA[J,I])
          end;
          if M<>0 then                                                       
            for I:=1 to M do
            begin                                                
              {BB[J,I]:=BB[J,I]-CS[J]*BB[K,I]; }
              CMUL(CS[J],BB[K,I],temp);
              CDIF(BB[J,I],temp,BB[J,I])
            end;
        fin: end                                                                   
{Line K is ready.}
      end; {of K loop}

{MATRIX AA INVERSION IS DONE - REARRANGEMENT OF COMPLEX MATRIX AA
                                                                               
 EXCHANGE LINES                }                                                            
      for I:=N downto 1 do
      begin                                                               
        IK:=Round(PC[I]);
        if IK=I then goto fin1;                                                   
{EXCHANGE LINES I and PC(I) of matrix AA:}                                         
        for J:=1 to N do CSwap(AA[I,J],AA[IK,J]);
        if M<>0 then                                                         
          for J:=1 to M do CSwap(BB[I,J],BB[IK,J]);
{NO MORE EXCHANGE IS NECESSARY                                                      
 Go to next line.}                                                  
      fin1: end; {of loop i=N downto }                                                                     
                                                                               
{EXCHANGE COLUMNS  }                                                          
      
      for J:=N downto 1 do
      begin                                                                         
        JK:=Round(PL[J]);                                                                
        if JK=J then goto fin2;                                                   
{EXCHANGE COLUMNS J ET PL(J) of matrix AA: }                                       
        for I:=1 to N do CSwap(AA[I,J],AA[I,JK]);
{NO MORE EXCHANGE IS NECESSARY                                                      
 Go to next column.}
      fin2: end;                                                                     
{REARRANGEMENT TERMINATED.  }                                                        
END;   {of procedure CMATINV }                                                                     


Procedure MATMUL(A, B:MATC; VAR C: MATC; N: integer);
{******************************************                                     
*  MULTIPLICATION OF TWO SQUARE COMPLEX   *                                     
*  MATRICES                               *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*N                *                                     
*            B  MATRIX N*N                *                                     
*            N  INTEGER                   *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*N PRODUCT A*B    *                                     
*                                         *
******************************************}
VAR
      SUM,PROD: Complex;
      I,J,K : integer;
BEGIN                                               
      for I:=1 to N do                                                                  
        for J:=1 to N do
        begin                                                                
          SUM.R:=0.0; SUM.I:=0.0;
          for K:=1 to N do
          begin                                                              
            {SUM:=SUM+A[I,K]*B[K,J] }
            CMUL(A[I,K],B[K,J],PROD);
            CADD(SUM,PROD,SUM)
          end;                                                 
          C[I,J]:=SUM                                                            
        end                                                                   
END;

Procedure MATMUL1(A:MATC; B:MATC1; VAR C:MATC1; N,M: integer);
{******************************************                                     
* MULTIPLICATION OF TWO COMPLEX MATRICES  *
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
      SUM,PROD,TMP: Complex;
      I,J,K : integer;
BEGIN                                               
      for I:=1 to N do                                                                  
        for J:=1 to M do
        begin                                                                
          SUM.R:=0.0; SUM.I:=0.0;
          for K:=1 to N do
          begin                                                              
            {SUM:=SUM+A[I,K]*B[K,J] }
            CMUL(A[I,K],B[K,J],PROD);
            CADD(SUM,PROD,TMP);
            SUM:=TMP
          end;                                                 
          C[I,J]:=SUM                                                            
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

      Assign(fp1,'cmat4.dat'); reset(fp1);
      Assign(fp2,'csysmat.lst'); rewrite(fp2);

      read(fp1,N);

      FWriteHead(fp2,'  SOLVING THE COMPLEX MATRIX LINEAR SYSTEM AX = B');
      writeln(fp2,'    (BY GAUSS-JORDAN METHOD WITH FULL PIVOTING)');
      ligne;
      writeln(fp2,'  N=',N);
      ligne;

      for I:=1 to N do
        for J:=1 to N do
        begin
          read(fp1,temp.R, temp.I);
          A[I,J] := temp
        end;

      writeln(fp2,'  COMPLEX MATRIX A:');
      ligne;

      for I:=1 to N do
        for J:=1 to N do
        begin
          write(fp2,' (',A[I,J].R:9:6,',',A[I,J].I:9:6,')');
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
            read(fp1,temp.R, temp.I); 
            B[I,J] := temp
          end;     

        writeln(fp2,'  COMPLEX MATRIX B:');
        ligne;

        for I:=1 to N do
          for J:=1 to M do
          begin
            write(fp2,' (',B[I,J].R:9:6,',',B[I,J].I:9:6,')');
            if J=M then ligne
          end
      end;
{END SECTION READ DATA}      
      close(fp1);                                                                     
                                                                               
{STORE MATRIX A IN MATRIX A1}

      for I:=1 to N do
        for J:=1 to N do
          A1[I,J]:=A[I,J];
                                                                               
{CALL COMPLEX MATRIX INVERSION ROUTINE}                                               
                                                                               
     CMATINV(N,M,A,B,DETER);
                                                                               
{PRINT RESULTS TO OUTPUT FILE}                                                          

      ligne;                                                                         
      writeln(fp2,'  INVERSE OF COMPLEX MATRIX A:');
      ligne;   
      for I:=1 to N do
        for J:=1 to N do
        begin
          write(fp2,' (',A[I,J].R:9:6,',',A[I,J].I:9:6,')');
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
            write(fp2,' (',B[I,J].R:9:6,',',B[I,J].I:9:6,')');
            if J=M then ligne
          end
      end;
      ligne;
      writeln(fp2,'  DETERMINANT= (',DETER.R:9:6,',',DETER.I:9:6,')');
      ligne;

{VERIFY THAT NEW PRODUCT A*B EQUALS ORIGINAL B MATRIX
 AND THAT A * INV(A) = I  (OPTIONAL)     }

      ligne;                                                                    
      writeln(fp2,'  VERIFICATION OF A1 * A = I:');
      ligne;
{CALL MULTIPLICATION ROUTINE}                                         
      MATMUL(A1,A,D,N);   
      for I:=1 to N do
        for J:=1 to N do
        begin
          write(fp2,' (',D[I,J].R:9:6,',',D[I,J].I:9:6,')');
          if J=N then ligne
        end;  

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
            write(fp2,' (',C[I,J].R:9:6,',',C[I,J].I:9:6,')');
            if J=M then ligne
          end
      end;

      ligne;
      ligne;
      FWriteEnd(fp2,'  End of results file.');
      close(fp2);
      writeln;
      writeln('  Results in file csysmat.lst...');
      readln;
      DoneWinCrt
END.                                                                          

{end of file csysmat.pas j-p moreau Feb. 2010}