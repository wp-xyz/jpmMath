{****************************************************************************
*     Step by step solution of system [M] X" + [C] X' + [K] X = F(t)        *
*                    by the "Wilson-Theta" Method                           *
* ------------------------------------------------------------------------- *
* SAMPLE RUN:                                                               *
*                                                                           *
* Let us consider the 3 d.o.f. Mass-Spring system with viscous damping:     *
*                                                                           *
*                      --> x1          --> x2    k     --> x3               *
*          \|    k     *----*    2k    *----*--/\/\/\--*----*               *
*          \|--/\/\/\--| 2m |--/\/\/\--| m  |  ---*    | 3m |               *
*          \|          *----*          *----*---| |----*----*               *
*                      --> F1(t)               ---*                         *
*                                                c                          *
*                      | 2  0  0 |                                          *
* The Mass Matrix M is | 0  1  0 |  assuming m=1.                           *
*                      | 0  0  3 |                                          *
*                                                                           *
*                           |  3 -2  0 |                                    *
* The Stiffness Matrix K is | -2  3 -1 |  assuming k=1.                     *
*                           |  0 -1  1 |                                    *
*                                                                           *
*                       | 0  0  0   |                                       *
* The Damping Matrix is | 0  0  0   | assuming c=0.5.                       *
*                       | 0  0  0.5 |                                       *
*                                                                           *
* The resonance fresuencies (neglecting c) are:                             *
*             F1 = w1/2pi = 0.05161 sqr(k/m)                                *
*             F2 = w2/2pi = 0.1431 sqr(k/m)                                 *
*             F3 = w3/2pi = 0.3151 sqr(k/m)                                 *
* with the corresponding modal vectors:                                     *
*             |   1   |         |   1    |         |    1    |              *
*      phi1 = | 1.395 |  phi2 = | 0.6914 |  phi3 = | -2.420  |              *
*             | 2.038 |         |-0.4849 |         |  0.2249 |              *
*                                                                           *
* (See program ndof03.pas).                                                 *
*                                                                           *
*                                                                           *
* How many degrees of freedom (d.o.f.): 3                                   *
*                                                                           *
* Mass Matrix:                                                              *
*  M(1,1) = 2                                                               *
*  M(1,2) = 0                                                               *
*  M(1,3) = 0                                                               *
*  M(2,2) = 1                                                               *
*  M(2,3) = 0                                                               *
*  M(3,3) = 3                                                               *
*                                                                           *
* Damping Matrix:                                                           *
*  C(1,1) = 0                                                               *
*  C(1,2) = 0                                                               *
*  C(1,3) = 0                                                               *
*  C(2,2) = 0                                                               *
*  C(2,3) = 0                                                               *
*  C(3,3) = 0.5                                                             *
*                                                                           *
* Stiffness Matrix:                                                         *
*  K(1,1) = 3                                                               *
*  K(1,2) = -2                                                              *
*  K(1,3) = 0                                                               *
*  K(2,2) = 3                                                               *
*  K(2,3) = -1                                                              *
*  K(3,3) = 1                                                               *
*                                                                           *
* Starting time (sec.) = 0                                                  *
* Ending time......... = 60                                                 *
* Time increment...... = 0.1                                                *
*                                                                           *
* Starting motion X(1) = 0                                                  *
* Starting motion X(2) = 0                                                  *
* Starting motion X(3) = 0                                                  *
*                                                                           *
* Starting speed V(1) = 0                                                   *
* Starting speed V(2) = 0                                                   *
* Starting speed V(3) = 0                                                   *
*                                                                           *
* Starting acceleration G(1) = 0                                            *
* Starting acceleration G(2) = 0                                            *
* Starting acceleration G(3) = 0                                            *
*                                                                           *
* Theta = 1.4                                                               *
*                                                                           *
* Number of force componant: 1                                              *
*                                                                           *
* Force (maximum) = 1000                                                    *
* Force frequency = 0.315                                                   *
*                                                                           *
* Output file ndof04.txt contains:                                          *
*                                                                           *
* Time= 0.00000000000000E+0000 Force= 0.00000000000000E+0000                *
* Node  Displacement  Speed  Acceleration                                   *
* ---------------------------------------                                   *
* 1  0.00000000000000E+0000  0.00000000000000E+0000  0.00000000000000E+0000 *
* 2  0.00000000000000E+0000  0.00000000000000E+0000  0.00000000000000E+0000 *
* 3  0.00000000000000E+0000  0.00000000000000E+0000  0.00000000000000E+0000 *
* Time= 1.00000000000000E-0001 Force= 1.96630694615420E+0002                *
* 1  1.63063348145476E-0001  4.89190044436429E+0000  9.78380088872859E+0001 *
* 2  1.05501179771710E-0003  3.16503539315131E-0002  6.33007078630263E-0001 *
* 3  1.13432172045768E-0006  3.40296516137303E-0005  6.80593032274606E-0004 *
* Time= 2.00000000000000E-0001 Force= 3.85583992277397E+0002                *
* 1  1.34210952239292E+0000  2.06956838943306E+0001  2.18237660112040E+0002 *
* 2  1.21558298960262E-0002  2.38073481154732E-0001  3.49545546583412E+0000 *
* 3  1.68196325261618E-0005  3.68470369329933E-0004  6.00822132204945E-0003 *
* Time= 3.00000000000000E-0001 Force= 5.59482258102167E+0002                *
* 1  4.69708708325055E+0000  4.83460760314658E+0001  3.34770182630666E+0002 *
* 2  6.53333107705596E-0002  9.44404690634833E-0001  1.06311687237679E+0001 *
* 3  1.20390195018308E-0004  2.06976507000204E-0003  2.80176726913926E-0002 *
* Time= 4.00000000000000E-0001 Force= 7.11535677209285E+0002                *
* 1  1.13732336500405E+0001  8.68537358092338E+0001  4.35383012924695E+0002 *
* 2  2.35184676961172E-0001  2.67517316826032E+0000  2.39842008287418E+0001 *
* 3  5.75234806581125E-0004  8.10492457231081E-0003  9.26855173547828E-0002 *
* Time= 5.00000000000000E-0001 Force= 8.35807361368270E+0002                *
* 1  2.23645614383329E+0001  1.34263211384070E+0002  5.12806498572030E+0002 *
* 2  6.57745512897410E-0001  6.12726870012940E+0000  4.50577098086397E+0001 *
* 3  2.10380937902387E-0003  2.50131121609216E-0002  2.45478234417433E-0001 *
* Time= 6.00000000000000E-0001 Force= 9.27445153334661E+0002                *
* 1  3.84369149652071E+0001  1.88003858109483E+0002  5.62006435936224E+0002 *
* 2  1.54516272061255E+0000  1.21150933407633E+0001  7.46987830040385E+0001 *
* 3  6.34915625419126E-0003  6.50602702123068E-0002  5.55464926610270E-0001 *
* -------------------------  ----------------------  ---------------------- *
* Time= 5.96000000000006E+0001 Force=-9.88651744737742E+0002                *
* 1  8.85667274205269E+0002 -4.46171811953303E+0003 -3.52638052370999E+0003 *
* 2 -1.55015816910563E+0003  1.04323928291804E+0004  6.41528284413283E+0003 *
* 3  9.27486031299907E+0001 -9.41879364154716E+0002 -3.83685509209653E+0002 *
* Time= 5.97000000000006E+0001 Force=-9.39811951085928E+0002                *
* 1  4.24805036107209E+0002 -4.72611187769025E+0003 -1.76149463943446E+0003 *
* 2 -4.81806425736124E+0002  1.08650025005178E+0004  2.23691058261484E+0003 *
* 3 -2.71286527277457E+0000 -9.60901048313042E+0002  3.25182604315031E+0000 *
* Time= 5.98000000000006E+0001 Force=-8.54277431698698E+0002                *
* 1 -5.35496823847710E+0001 -4.81034306740717E+0003  7.68708450960193E+0001 *
* 2  6.08748956281694E+0002  1.08748109303682E+0004 -2.04074198560662E+0003 *
* 3 -9.81404593989860E+0001 -9.41188318462418E+0002  3.91002770969319E+0002 *
* Time= 5.99000000000006E+0001 Force=-7.35387860780236E+0002                *
* 1 -5.31131703638733E+0002 -4.71061804505930E+0003  1.91762960186138E+0003 *
* 2  1.67900827481359E+0003  1.04601947945007E+0004 -6.25158073174253E+0003 *
* 3 -1.89681955023360E+0002 -8.83418370354851E+0002  7.64396191182026E+0002 *
* Time= 6.00000000000006E+0001 Force=-5.87785252291540E+0002                *
* 1 -9.89652416589468E+0002 -4.43026677849652E+0003  3.68939572939423E+0003 *
* 2  2.68713644244187E+0003  9.63603447643417E+0003 -1.02316256295891E+0004 *
* 3 -2.73627842586534E+0002 -7.89759695744624E+0002  1.10877730102250E+0003 *
*                                                                           *
* Maximum Displacement (node #1) =  2.43218400402500E+0003                  *
* ------------------------------------------------------------------------- *
* REFERENCE:  "Mécanique des vibrations linéaires By M. Lalanne,            *
*              P. Berthier, J. Der Hagopian, Masson, Paris 1980"            *
*              [BIBLI 16].                                                  *
*                                                                           *
*                                     Pascal Release By J-P Moreau, Paris.  *
*                                              (www.jpmoreau.fr)            *
****************************************************************************}
Program ndof04;

Uses WinCrt;

Label 500;

Const
      NMAX = 12;
      EPSMACH = 2e-16;

Type
      MAT = Array[1..NMAX,1..NMAX] of Double;
      MAT1= Array[1..NMAX,1..2] of Double;
      VEC = Array[1..NMAX] of Double;

Var
      M, K, K1, A, C: MAT;
      B1: MAT1;
      X0, V0, G0, F0, F1, F2, X1, V1, G1, X2: VEC;

      D,DET,F,Fmax,freq,SUM,T0,T1,T2,T3: Double;
      A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,XMAX: Double;
      I,J,N,N1:Integer;

      fp: TEXT;


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


BEGIN

  Assign(fp, 'ndof4.txt'); Rewrite(fp);

  Writeln;
  Write(' How many degrees of freedom (d.o.f.): '); Readln(N);

  {Data section}
  Writeln;
  Writeln(' Mass Matrix:');
  FOR I := 1 TO N DO
    FOR J := I TO N DO
    begin
      Write(' M(', I, ',', J, ') = '); Readln(M[I, J])
    end;

  Writeln;
  Writeln(' Damping Matrix:');
  FOR I := 1 TO N DO
    FOR J := I TO N DO
    begin
      Write(' C(', I, ',', J, ') = '); Readln(C[I, J])
    end;

  Writeln;
  Writeln(' Stiffness Matrix:');
  FOR I := 1 TO N DO
    FOR J := I TO N DO
    begin
      Write(' K(', I, ',', J, ') = '); Readln(K[I, J])
    end;

  {Complete Matrices by symmetry}
  FOR I := 1 TO N DO
    FOR J := I TO N DO
    begin
      K[J, I] := K[I, J];
      C[J, I] := C[I, J];
      M[J, I] := M[I, J]
    end;

  Writeln;
  Write(' Starting time (sec.) = '); Readln(T0);
  Write(' Ending time......... = '); Readln(T3);
  Write(' Time increment...... = '); Readln(D);

  Writeln;
  FOR I := 1 TO N DO
  begin
    Write(' Starting motion X(', I, ') = '); Readln(X0[I])
  end;
  Writeln;
  FOR I := 1 TO N DO
  begin
    Write(' Starting speed V(', I, ') = '); Readln(V0[I])
  end;
  Writeln;
  FOR I := 1 TO N DO
  begin
    Write(' Starting acceleration G(', I, ') = '); Readln(G0[I])
  end;
  Writeln;
  Write(' Theta = '); Readln(T2);
  Writeln;
  Write(' Number of force componant: '); Readln(N1);
  Writeln;
  Write(' Force (maximum) = '); Readln(Fmax);
  Write(' Force frequency = '); Readln(freq);
  F := Fmax * SIN(2.0 * PI * freq * T0);
  {end of data section}

  Writeln(fp);
  Writeln(fp, ' Time=', T0, ' Force=', F);
  Writeln(fp, ' Node  Displacement  Speed  Acceleration');
  Writeln(fp, ' ---------------------------------------');
  FOR I := 1 TO N DO
    Writeln(fp, ' ', I, ' ', X0[I], ' ', V0[I], ' ', G0[I]);
  Writeln;
  FOR I := 1 TO N DO F0[I] := 0.0;

  F0[N1] := F;
  A9 := T2 * D;
  A0 := 6.0 / (A9*A9);
  A1 := 3.0 / A9;
  A2 := 2.0 * A1;
  A3 := A9 / 2.0;
  A4 := A0 / T2;
  A5 := -A2 / T2;
  A6 := 1.0 - 3.0 / T2;
  A7 := D / 2.0;
  A8 := D*D / 6.0;

{ Construct K+A0*M+A1*C
  --------------------- }
500: FOR I := 1 TO N DO F1[I] := 0.0;
  T1 := T0 + D;
  F := Fmax * SIN(2.0 * PI * freq * T1);
  Writeln(fp, ' Time=', T1, ' Force=', F);
  F1[N1] := F;
  FOR I := 1 TO N DO F2[I] := 0.0;

  FOR I := 1 TO N DO
  begin
    F2[I] := F0[I] + T2 * (F1[I] - F0[I]);
    FOR J := 1 TO N DO
    begin
      K1[I, J] := K[I, J] + A0 * M[I, J] + A1 * C[I, J];
      F2[I] := F2[I] + M[I, J] * (A0 * X0[J] + A2 * V0[J] + 2.0 * G0[J]);
      F2[I] := F2[I] + C[I, J] * (A1 * X0[J] + 2.0 * V0[J] + A3 * G0[J])
    end
  end;

  {K1:=INV(K1) }
  MATINV(N,0,K1,B1,DET);

  {X2:=K1 MPY F2 }
  FOR I := 1 TO N DO
  begin
    SUM := 0.0;
    FOR J := 1 TO N DO
      SUM := SUM + K1[I, J] * F2[J];
    X2[I] := SUM
  end;

  FOR I := 1 TO N DO
  begin
    G1[I] := A4 * (X2[I] - X0[I]) + A5 * V0[I] + A6 * G0[I];
    V1[I] := V0[I] + A7 * (G1[I] + G0[I]);
    X1[I] := X0[I] + D * V0[I] + A8 * (G1[I] + 2.0 * G0[I]);
    IF ABS(X1[1]) > XMAX THEN XMAX := X1[1]
  end;

  FOR I := 1 TO N DO
    Writeln(fp, ' ', I, ' ', X1[I], ' ', V1[I], ' ', G1[I]);
  FOR I := 1 TO N DO
  begin
    X0[I] := X1[I]; V0[I] := V1[I]; G0[I] := G1[I]
  end;
  T0 := T1;
  IF T1 < T3 THEN GOTO 500;  {next time step}

  Writeln;
  Writeln(' Results in file ndof4.txt...');
  Writeln;
  Writeln(fp);
  Writeln(fp, ' Maximum Displacement (node #1) = ', XMAX);
  Close(fp);

  ReadKey;
  DoneWinCrt

END.

{end of file ndof4.pas}