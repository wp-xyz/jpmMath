(********************************************************************
*   Frequencies and eigenmodes, masses and modal stiffnesses of a   *
*   Mass-Spring undamped system represented by Motion Equation:     *
*                 [M] . [X]" + [K] . [X] = [0}                      *
* ----------------------------------------------------------------- *
* EXPLANATION:                                                      *
* Let us take the 3 D.O.F. Mass-Spring System below:                *
*                                                                   *
*     \|    k     ------    2k    ------    k     ------            *
*     \|--/\/\/\--| 2m |--/\/\/\--| m  |--/\/\/\--| 3m |            *
*     \|          ------          ------          ------            *
*                                                                   *
*                      | 2  0  0 |                                  *
* The Mass Matrix M is | 0  1  0 |  assuming m=1.                   *
*                      | 0  0  3 |                                  *
*                                                                   *
*                           |  3 -2  0 |                            *
* The Stiffness Matrix K is | -2  3 -1 |  assuming k=1.             *
*                           |  0 -1  1 |                            *
*                                                                   *
* The resonance fresuencies are:                                    *
*             F1 = w1/2pi = 0.05161 sqr(k/m)                        *
*             F2 = w2/2pi = 0.1431 sqr(k/m)                         *
*             F3 = w3/2pi = 0.3151 sqr(k/m)                         *
* with the corresponding modal vectors:                             *
*             |   1   |         |   1    |         |    1    |      *
*      phi1 = | 1.395 |  phi2 = | 0.6914 |  phi3 = | -2.420  |      *
*             | 2.038 |         |-0.4849 |         |  0.2249 |      *
*                                                                   *
* (See program ndof01.bas).                                         *
* ----------------------------------------------------------------- *
* SAMPLE RUN:                                                       *
*                                                                   *
* Order of system: 3                                                *
*                                                                   *
* M(1,1) = 2                                                        *
* M(1,2) = 0                                                        *
* M(1,3) = 0                                                        *
* M(2,2) = 1                                                        *
* M(2,3) = 0                                                        *
* M(3,3) = 3                                                        *
*                                                                   *
* K(1,1) = 3                                                        *
* K(1,2) = -2                                                       *
* K(1,3) = 0                                                        *
* K(2,2) = 3                                                        *
* K(2,3) = -1                                                       *
* K(3,3) = 1                                                        *
*                                                                   *
* Number of eigenvectors asked for: 3                               *
*                                                                   *
* What precision for eigenvalues: 1e-5                              *
*                                                                   *
* Maximum number of iterations: 100                                 *
*                                                                   *
*                                                                   *
* Eigenvector #1, Convergence after 10 iterations.                  *
* Eigenvector #2, Convergence after 9 iterations.                   *
* Eigenvector #3, Convergence after 2 iterations.                   *
*                                                                   *
*                                                                   *
* Eigenvalues:                                                      *
*                                                                   *
* L(1) =  1.05173321566663E-0001                                    *
* L(2) =  8.08612063101068E-0001                                    *
* L(3) =  3.91945398307909E+0000                                    *
*                                                                   *
* Pulsations:                                                       *
*                                                                   *
* W(1) =  3.24304365629980E-0001                                    *
* W(2) =  8.99228593351584E-0001                                    *
* W(3) =  1.97976109242481E+0000                                    *
*                                                                   *
* Frequencies:                                                      *
*                                                                   *
* F(1) =  5.16146428562926E-0002                                    *
* F(2) =  1.43116675601476E-0001                                    *
* F(3) =  3.15088764000420E-0001                                    *
*                                                                   *
* Eigenvectors:                                                     *
*                                                                   *
* E.V.  1                                                           *
*              1.00000000000000E+0000                               *
*              1.39482667843334E+0000                               *
*              2.03778199949639E+0000                               *
*                                                                   *
* E.V.  2                                                           *
*              1.00000000000000E+0000                               *
*              0.69138346664719E+0000                               *
*             -0.48490130871242E+0000                               *
*                                                                   *
* E.V.  3                                                           *
*              1.00000000000000E+0000                               *
*             -2.41959287715014E+0000                               *
*              2.24903464725100E+0000                               *
*                                                                   *
* Modal Mass Matrix:                                                *
*                                                                   *
*  16.40321  -0.00001   0.00000                                     *
*  -0.00001   3.18340  -0.00003                                     *
*   0.00000  -0.00003   8.00617                                     *
*                                                                   *
* Modal Stiffness Matrix:                                           *
*                                                                   *
*   1.72517   0.00000  -0.00000                                     *
*   0.00000   2.57413   0.00001                                     *
*  -0.00000   0.00001  31.38059                                     *
*                                                                   *
* ----------------------------------------------------------------- *
* REFERENCE: "Mécanique des vibrations linéaires By M. Lalanne,     *
*            P. Berthier, J. Der Hagopian, Masson, Paris 1980"      *
*            [BIBLI 16].                                            *
*                                                                   *
*                             Pascal Release By J-P Moreau, Paris.  *
*                                      (www.jpmoreau.fr)            *
********************************************************************)
Program NDOF02;

Uses WinCrt;

Label 800,820,970,1070,1160,1290,1490;

Const

  EPSMACH = 2E-16;
  NMAX = 25;

Type

  pMAT = ^MAT;
  MAT = Array[1..NMAX,1..NMAX] of DOUBLE;
  pVEC = ^VEC;
  VEC = Array[1..NMAX] of DOUBLE;

Var

  M, K: pMAT;
  A, B, C, D, K1, S, R, R1, R2, S1, T, T1, V: pMAT;
  L, X, F, F1: pVEC;

  A0,A1,A2,D1,DET,M1,P,SUM,W2: Double;
  I,J,J1,K0,K8,K9,N,N1,N2,N3: Integer;


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
Procedure MATINV(N,M:integer; VAR AA:pMAT; VAR BB:pMAT;VAR DET:DOUBLE);
VAR   PC,PL,CS : pVEC;                                                  
      PV,PAV,temp,TT : DOUBLE;
      I,IK,J,JK,K : integer;

LABEL fin, fin1, fin2;

BEGIN                                                         
{Initializations}
      New(PC); New(PL); New(CS);
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
          Dispose(PC); Dispose(PL); Dispose(CS);
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
                                                                           
        if M<>0 then                                                           
          for I:=1 to M do
          begin                                                               
            TT:=BB^[IK][I];
            BB^[IK][I]:=BB^[K][I];
            BB^[K][I]:=TT
          end;                                                                   

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
          Dispose(PC); Dispose(PL); Dispose(CS);
          halt                                                                  
        end;                                                                   
        for I:=1 to N do                                                                
          AA^[K][I]:=AA^[K][I]/PV;
                                                                           
        if M<>0 then                                                         
          for I:=1 to M do                                                             
            BB^[K][I]:=BB^[K][I]/PV;
                                                                           
{Other lines of matrix AA are modified:}                                        
        for J:=1 to N do
        begin                                                                
          if J=K then goto fin;                                                  
          for I:=1 to N do
{Line J of matrix AA is modified:}                                            
            AA^[J][I]:=AA^[J][I]-CS^[J]*AA^[K][I];
          if M<>0 then                                                       
            for I:=1 to M do                                                          
              BB^[J][I]:=BB^[J][I]-CS^[J]*BB^[K][I];
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
        if M<>0 then                                                         
          for J:=1 to M do
          begin
            TT:=BB^[I][J];
            BB^[I][J]:=BB^[IK][J];
            BB^[IK][J]:=TT                                                         
          end;
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
      Dispose(PC); Dispose(PL); Dispose(CS);
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
Procedure MATMUL(A, B:pMAT; VAR C: pMAT; N: integer);
VAR
      SUM : DOUBLE;
      I,J,K : integer;
BEGIN                                               
  for I:=1 to N do                                                                  
    for J:=1 to N do
    begin                                                                
      SUM:= 0.0;                                                                
      for K:=1 to N do                                                              
        SUM:=SUM+A^[I,K]*B^[K,J];                                               
      C^[I,J]:=SUM
    end                                                                   
END;                                                                       


{main program}
BEGIN

 New(M); New(K);  New(A); New(B); New(C); New(D);
 New(K1); New(S); New(R); New(R1); New(R2); New(S1);
 New(T); New(T1); New(V);
 New(L); New(X); New(F); New(F1);

 Writeln;
 Write(' Order of system: '); Readln(N);
 Writeln;

 FOR I := 1 TO N DO
   FOR J := I TO N DO
   begin
     Write(' M(', I, ',', J, ') = '); Readln(M^[I, J])
   end;
 Writeln;

 FOR I := 1 TO N DO
   FOR J := I TO N DO
   begin
     Write(' K(', I, ',', J, ') = '); Readln(K^[I, J])
   end;

 A1 := 0.0; A2 := 0.0;
 FOR I := 1 TO N DO
   FOR J := I TO N DO
   begin
     A1 := A1 + M^[I, I];
     A2 := A2 + K^[I, I];
     M^[J, I] := M^[I, J];
     K^[J, I] := K^[I, J]
   end;

 Writeln;
 Write(' Number of eigenvectors asked for: '); Readln(N1);

 K8 := 1;

 Writeln;
 Write(' What precision for eigenvalues: '); Readln(P);
 P := P * P;
 Writeln;
 Write(' Maximum number of iterations: '); Readln(N3);
 Writeln;

 A0 := A2 / A1 / 10.0;
 FOR I := 1 TO N DO
   FOR J := 1 TO N DO
     K1^[I, J] := K^[I, J] + A0 * M^[I, J];

 {K1:=INV(K1) }
 MATINV(N,0,K1,S,DET);

 {D:=K1 MPY M }
 MATMUL(K1, M, D, N);

 {V=0 }
  FOR I := 1 TO N DO
   FOR J := 1 TO N1 DO
     V^[I, J] := 0.0;

 {S1=D }
 FOR I := 1 TO N DO
   FOR J := 1 TO N DO
     S1^[I, J] := D^[I, J];

800: FOR I := 1 TO N DO X^[I] := 1.0;

 K9 := 0;

820: K9 := K9 + 1;
 {F:=S1 MPY X}
 FOR I := 1 TO N DO
 begin
   SUM := 0.0;
   FOR J := 1 TO N DO
     SUM := SUM + S1^[I, J] * X^[J];
   F^[I] := SUM
 end;

 W2 := 1.0 / F^[1];
 L^[K8] := W2 - A0;
 FOR I := 1 TO N DO F^[I] := F^[I] * W2;
 FOR I := 1 TO N DO F1^[I] := F^[I] - X^[I];
 D1 := 0.0; M1 := 0.0;
 FOR I := 1 TO N DO
 begin
   D1 := D1 + F1^[I] * F1^[I];
   M1 := M1 + F^[I] * F^[I];
   X^[I] := F^[I]
 end;

970: IF K9 - N3 <= 0 THEN GOTO 1070;
 Writeln(' No convergence after', K9, ' Iterations.');
 Halt;

1070: IF P - D1 / M1 < 0.0 THEN GOTO 820;
 Writeln(' Eigenvector #', K8, ', Convergence after', K9, ' iterations.');
 FOR I := 1 TO N DO V^[I, K8] := X^[I];

 IF K8 >= N1 THEN GOTO 1490;  {print results}

1160: N2 := N - K8;
 {T=TRN(V) }
 FOR I := 1 TO N DO
   FOR J := 1 TO N1 DO
     T^[J, I] := V^[I, J];

 {T1:=T MPY M}
 MATMUL(T, M, T1, N);

 FOR I := 1 TO K8 DO
   FOR J := 1 TO N DO
   begin
     IF J <= K8 THEN
     begin
       R^[I, J] := T1^[I, J];
       GOTO 1290
     end;
     J1 := J - K8;
     R1^[I, J1] := T1^[I, J];
1290:end;

 {R=INV(R) }
 MATINV(K8,0,R,S,DET);

 {R2(N,N2):=R(N,K8) MPY R1(N,N2) }
 FOR I := 1 TO K8 DO
   FOR J := 1 TO N2 DO
   begin
     SUM := 0.0;
     FOR K0 := 1 TO K8 DO
     begin
       SUM := SUM + R^[I, K0] * R1^[K0, J];
       R2^[I, J] := SUM
     end
   end;

  {S=Identity Matrix}
 FOR I := 1 TO N DO
   FOR J := 1 TO N DO
     IF J = I THEN
       S^[I, J] := 1.0
     ELSE
       S^[I, J] := 0.0;

 FOR I := 1 TO K8 DO
 begin
   S^[I, I] := 0.0;
   FOR J := K8 + 1 TO N DO
   begin
     J1 := J - K8;
     S^[I, J] := -R2^[I, J1]
   end
 end;

 {S1:=D MPY S }
 MATMUL(D, S, S1, N);

 K8 := K8 + 1;
 GOTO 800;  {go to next eigenvalue/eigenvector}

{results section}
1490: Writeln;
 Writeln(' Eigenvalues:');
 Writeln;
 FOR I := 1 TO N1 DO
   Writeln(' L(', I, ') = ', L^[I]);
 Writeln;
 Writeln(' Pulsations:');
 Writeln;
 FOR I := 1 TO N1 DO 
   Writeln(' W(', I, ') = ', SQRT(L^[I]));
 Writeln;
 Writeln(' Frequencies:');
 Writeln;
 FOR I := 1 TO N1 DO
   Writeln(' F(', I, ') = ', SQRT(L^[I]) / 2.0 / PI);
 Writeln;
 ReadKey;
 Writeln;
 Writeln(' Eigenvectors:');
 Writeln;
 FOR J := 1 TO N1 DO
 begin
   Writeln;
   Writeln(' E.V. ', J);
   FOR I := 1 TO N DO
     Writeln('             ', V^[I, J])
 end;
 Writeln;
 ReadKey;

 {S1(N,N1):=M(N,N) MPY V(N,N1) }
 FOR I := 1 TO N DO
   FOR J := 1 TO N1 DO
   begin
     SUM := 0.0;
     FOR K0 := 1 TO N DO
     begin
       SUM := SUM + M^[I, K0] * V^[K0, J];
       S1^[I, J] := SUM
     end
   end;

 {T1:=TRN(V) }
 FOR I := 1 TO N DO
   FOR J := 1 TO N1 DO
     T1^[J, I] := V^[I, J];

 {R(N1,N1):=T1(N1,N) MPY S1(N,N1) }
  FOR I := 1 TO N1 DO
   FOR J := 1 TO N1 DO
   begin
     SUM := 0.0;
     FOR K0 := 1 TO N DO
     begin
       SUM := SUM + T1^[I, K0] * S1^[K0, J];
       R^[I, J] := SUM
     end
   end;

 {Writeln modal mass matrix}
 Writeln;
 Writeln(' Modal Mass Matrix:');
 Writeln;
 FOR I := 1 TO N1 DO
 begin
   FOR J := 1 TO N1 DO
     Write(R^[I, J]:10:5);
   Writeln
 end;

 {S1(N,N1):=K(N,N) MPY V(N,N1) }
 FOR I := 1 TO N DO
   FOR J := 1 TO N1 DO
   begin
     SUM := 0.0;
     FOR K0 := 1 TO N DO
     begin
       SUM := SUM + K^[I, K0] * V^[K0, J];
       S1^[I, J] := SUM
     end
   end;

 {R(N1,N1):=T1(N1,N) MPY S1(N,N1) }
 FOR I := 1 TO N1 DO
   FOR J := 1 TO N1 DO
   begin
     SUM := 0.0;
     FOR K0 := 1 TO N DO
     begin
       SUM := SUM + T1^[I, K0] * S1^[K0, J];
       R^[I, J] := SUM
     end
   end;

 {Writeln modal stiffness matrix}
 Writeln;
 Writeln(' Modal Stiffness Matrix:');
 Writeln;
 FOR I := 1 TO N1 DO
 begin
   FOR J := 1 TO N1 DO
     Write(R^[I, J]:10:5);
   Writeln
 end;
 Writeln;

 ReadKey;
 Dispose(M); Dispose(K);  Dispose(A); Dispose(B); Dispose(C); Dispose(D);
 Dispose(K1); Dispose(S); Dispose(R); Dispose(R1); Dispose(R2); Dispose(S1);
 Dispose(T); Dispose(T1); Dispose(V); Dispose(L); Dispose(X); Dispose(F);
 Dispose(F1);
 DoneWinCrt

 END.

{end of file ndof02.pas}