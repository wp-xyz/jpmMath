{************************************************************************
*       Response of a N d.o.f. Mass-Spring System with Damping          *
*              to a sinusoidal Force By a Direct Method                 *
* --------------------------------------------------------------------- *
* EXPLANATION:                                                          *
* Let us consider the 3 d.o.f. Mass-Spring system with viscous damping: *
*                                                                       *
*                      --> x1          --> x2    k     --> x3           *
*          \|    k     *----*    2k    *----*--/\/\/\--*----*           *
*          \|--/\/\/\--| 2m |--/\/\/\--| m  |  ---*    | 3m |           *
*          \|          *----*          *----*---| |----*----*           *
*                      --> F1(t)               ---*                     *
*                                                c                      *
*                      | 2  0  0 |                                      *
* The Mass Matrix M is | 0  1  0 |  assuming m=1.                       *
*                      | 0  0  3 |                                      *
*                                                                       *
*                           |  3 -2  0 |                                *
* The Stiffness Matrix K is | -2  3 -1 |  assuming k=1.                 *
*                           |  0 -1  1 |                                *
*                                                                       *
*                       | 0  0  0   |                                   *
* The Damping Matrix is | 0  0  0   | assuming c=0.5.                   *
*                       | 0  0  0.5 |                                   *
*                                                                       *
* The resonance fresuencies (neglecting c) are:                         *
*             F1 = w1/2pi = 0.05161 sqr(k/m)                            *
*             F2 = w2/2pi = 0.1431 sqr(k/m)                             *
*             F3 = w3/2pi = 0.3151 sqr(k/m)                             *
* with the corresponding modal vectors:                                 *
*             |   1   |         |   1    |         |    1    |          *
*      phi1 = | 1.395 |  phi2 = | 0.6914 |  phi3 = | -2.420  |          *
*             | 2.038 |         |-0.4849 |         |  0.2249 |          *
*                                                                       *
* (See program ndof01.bas).                                             *
* --------------------------------------------------------------------- *
* SAMPLE RUN:                                                           *
*                                                                       *
* How many degrees of freedom (d.o.f.): 3                               *
*                                                                       *
* Input Mass Matrix [M]:                                                *
*  M(1,1) = 2                                                           *
*  M(1,2) = 0                                                           *
*  M(1,3) = 0                                                           *
*  M(2,2) = 1                                                           *
*  M(2,3) = 0                                                           *
*  M(3,3) = 3                                                           *
*                                                                       *
* Viscous Damping...: VIS                                               *
* Structural Damping: STR                                               *
*                                                                       *
* Your choice: VIS                                                      *
*                                                                       *
* Input Stiffness Matrix [K]:                                           *
*  K(1,1) = 3                                                           *
*  K(1,2) = -2                                                          *
*  K(1,3) = 0                                                           *
*  K(2,2) = 3                                                           *
*  K(2,3) = -1                                                          *
*  K(3,3) = 1                                                           *
*                                                                       *
* Input Damping Matrix [C]:                                             *
*  C(1,1) = 0                                                           *
*  C(1,2) = 0                                                           *
*  C(1,3) = 0                                                           *
*  C(2,2) = 0                                                           *
*  C(2,3) = 0                                                           *
*  C(3,3) = 0.5                                                         *
*                                                                       *
* Input Excitation Vector:                                              *
*  F1(1) = 1000                                                         *
*  F1(2) = 0                                                            *
*  F1(3) = 0                                                            *
*                                                                       *
*  F2(1) = 0                                                            *
*  F2(2) = 0                                                            *
*  F2(3) = 0                                                            *
*                                                                       *
* Number of d.o.f. to calculate: 3                                      *
*                                                                       *
* Frequency Sweep                                                       *
* ---------------                                                       *
* Starting Frequency: 0.31                                              *
* Ending Frequency..: 0.32                                              *
*                                                                       *
* Linear frequency step.....: LIN                                       *
* Logarithmic frequency step: LOG                                       *
*                                                                       *
* Your choice: LIN                                                      *
*                                                                       *
* Frequency Step....: 0.0005                                            *
*                                                                       *
*                                                                       *
* Frequency (Hz)   Displacement   Phase (deg.)                          *
* --------------------------------------------                          *
*     0.3100          240.66        0.0                                 *
*     0.3105          264.73        0.0                                 *
*     0.3110          294.72        0.0                                 *
*     0.3115          333.08        0.0                                 *
*     0.3120          383.85        0.0                                 *
*     0.3125          454.17        0.0                                 *
*     0.3130          557.90        0.0                                 *
*     0.3135          725.80        0.0                                 *
*     0.3140         1041.65        0.0                                 *
*     0.3145         1824.42        0.0                                 *
*     0.3150         4347.56        0.0  <-- 3rd resonance              *
*     0.3155         2248.74       -0.0                                 *
*     0.3160         1151.46       -0.0                                 *
*     0.3165          757.64       -0.0                                 *
*     0.3170          560.61       -0.0                                 *
*     0.3175          443.06       -0.0                                 *
*     0.3180          365.13       -0.0                                 *
*     0.3185          309.74       -0.0                                 *
*     0.3190          268.36       -0.0                                 *
*     0.3195          236.29       -0.0                                 *
*     0.3200          210.72       -0.0                                 *
*                                                                       *
* --------------------------------------------------------------------- *
* REFERENCE:  "Mécanique des vibrations linéaires By M. Lalanne,        *
*              P. Berthier, J. Der Hagopian, Masson, Paris 1980"        *
*              [BIBLI 16].                                              *
*                                                                       *
*                                 Pascal Release By J-P Moreau, Paris.  *
*                                          (www.jpmoreau.fr)            *
************************************************************************}
PROGRAM NDOF03;

Uses WinCrt;

Label 600,900,1400,1660,1750,1760,1800;

Const NMAX = 12;
      EPSMACH = 1.2E-16;

Type  MAT = Array[1..NMAX,1..NMAX] of Double;
      MAT1 = Array[1..2*NMAX,1..2] of Double;
      MAT2 = Array[1..2*NMAX,1..2*NMAX] of Double;
      VEC2 = Array[1..2*NMAX] of Double;

Var
     M, K1, K2: MAT; 
     A: MAT2;
     X, B: VEC2;
     B1: MAT1;

     C9,I,I1,J,J1,N,N3: Integer;
     DET,F,F1,F2,F3, M1,O1,S1,SUM,T1,W: Double;
     ANS:String[3];

     {y^x}
     Function Power(y,x:double): double;
     Begin
       Power:=0.0;
       IF x<0 THEN EXIT;
       Power:=Exp(x*Ln(y))
     End;


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
Procedure MATINV(N,M:integer; VAR AA:MAT2; VAR BB:MAT1;VAR DET:DOUBLE);
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


{main program}
BEGIN

  Writeln;
  Write(' How many degrees of freedom (d.o.f.): '); Readln(N);
  Writeln;
  Writeln(' Input Mass Matrix [M]:');
  Writeln;
  FOR I := 1 TO N DO
    FOR J := I TO N DO
    begin
      Write('  M(', I, ',', J, ') = '); Readln(M[I, J])
    end;

  Writeln;
  Writeln(' Viscous Damping...: VIS');
  Writeln(' Structural Damping: STR');
  Writeln;
  Write(' Your choice: '); Readln(ANS);

  IF ANS = 'VIS' THEN GOTO 600;

  {case structural damping}
  C9 := 2;
  Writeln;
  Writeln(' Input stiffness Matrix [K1]:');
  FOR I := 1 TO N DO
    FOR J := I TO N DO
    begin
      Write('  K1(', I, ',', J, ') = '); Readln(K1[I, J])
    end;

  Writeln;
  Writeln(' Input Damping Matrix [K2]:');
  FOR I := 1 TO N DO
    FOR J := I TO N DO
    begin
      Write('  K2(', I, ',', J, ') = '); Readln(K2[I, J])
    end;
  Writeln;
  GOTO 900;

  600: {case viscous damping}
  C9 := 1;
  Writeln;
  Writeln(' Input stiffness Matrix [K]:');
  FOR I := 1 TO N DO
    FOR J := I TO N DO
    begin
      Write('  K(', I, ',', J, ') = '); Readln(K1[I, J])
    end;
  Writeln;
  Writeln(' Input Damping Matrix [C]:');
  FOR I := 1 TO N DO
    FOR J := I TO N DO
    begin
      Write('  C(', I, ',', J, ') = '); Readln(K2[I, J])
    end;
  Writeln;

900:  {Complete Matrices by symmetry}
  FOR I := 1 TO N DO
    FOR J := I TO N DO
    begin
      K1[J, I] := K1[I, J];
      K2[J, I] := K2[I, J];
      M[J, I] := M[I, J]
    end;

  Writeln;
  Writeln(' Input Excitation Vector:');
  FOR I := 1 TO N DO
  begin
    Write('  F1(', I, ') = '); Readln(B[I])
  end;
  Writeln;
  FOR I := N + 1 TO 2*N DO
  begin
    J := I - N;
    Write('  F2(', J, ') = '); Readln(B[I])
  end;

  Writeln;
  Write(' Number of d.o.f. to calculate: '); Readln(N3);
  Writeln;

  {calculate numerical response of d.o.f. #N3}
  Writeln;
  Writeln(' Frequency Sweep');
  Writeln(' ---------------');
  Writeln;
1400: Write(' Starting frequency: '); Readln(F1);
  Write(' Ending frequency..: '); Readln(F2);
  Writeln;
  Writeln(' Linear frequency step.....: LIN');
  Writeln(' Logarithmic frequency step: LOG');
  Writeln;
  Write(' Your choice: '); Readln(ANS);
  Writeln;

  IF ANS = 'LOG' THEN GOTO 1660;

  {case linear step}
  Write(' Frequency step: '); Readln(F3);
  GOTO 1750;

1660:  {case log. step}
  IF F1 = 0.0 THEN
  begin
    Writeln(' CAUTION - Starting frequency must be > 0 for log. step.');
    GOTO 1400
  end
  ELSE
  begin
    Write(' Number of constant log. increments from F1 to F2: '); Readln(F3);
    F3 := Power(10.0, ((Ln(F2) - Ln(F1)) / Ln(10.0) / F3))
  end;

1750:  {Writeln; header}
  Writeln;
  Writeln(' Frequency [Hz)   Displacement   Phase [deg.)');
  Writeln(' --------------------------------------------');

  F := F1;  {Initialize frequency}

1760: W := 2.0 * PI * F;
  FOR I := 1 TO N DO
  begin
    I1 := I + N;
    FOR J := 1 TO N DO
    begin
      J1 := J + N;
      A[I, J] := K1[I, J] - W * W * M[I, J];
      A[I1, J1] := A[I, J];
      IF C9 = 1 THEN
        A[I1, J] := K2[I, J] * W
      ELSE
        A[I1, J] := K2[I, J];
      A[I, J1] := -A[I1, J]
    end
  end;

  {A=INV(A) }
  MATINV(2*N,0,A,B1,DET);

  {X:=A MPY B}
  FOR I := 1 TO 2*N DO
  begin
    SUM := 0.0;
    FOR J := 1 TO 2*N DO
      SUM := SUM + A[I, J] * B[J];
    X[I] := SUM
  end;

  M1 := SQRT(X[N3]*X[N3] + X[N3+N]*X[N3+N]);

  IF M1 = 0.0 THEN T1 := 0.0
  ELSE
  begin
    O1 := X[N3] / M1;
    S1 := -X[N3+N] / M1;
    IF O1 = -1.0 THEN T1 := PI
    ELSE IF O1 = 1.0 THEN T1 := 0.0
    ELSE
      IF S1 >= 0.0 THEN
        {T1:=ACOS(O1) }
        T1 := ArcTan(SQRT(1.0 - O1*O1) / O1)
      ELSE
        {T1=2*PI-ACOS[O1) }
        T1 := 2.0 * PI - ArcTan(SQRT(1.0 - O1*O1) / O1)
  end;

1800: T1 := T1 / PI / 180;  {Convert phase in degrees}
  {Writeln; line of results}
  Writeln('   ',F:8:4,'       ',M1:9:2,'      ',T1:5:1);

  IF ANS = 'LIN' THEN
    F := F + F3
  ELSE
    F := F * F3;

  IF F <= F2 THEN GOTO 1760;   {next frequency}

  Writeln;
  ReadKey;
  DoneWinCrt

END.

{end of file ndof03.pas}