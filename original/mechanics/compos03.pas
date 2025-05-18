{**************************************************************************************
*   This program calculates the matrix linking the stresses to deformamations in a    *
*   laminated material made of n unidirectional composite layers.                     *
*                                                                                     *
* The 6x6 stiffness matrix to define has the following coefficients:                  *
*                                                                                     *
*     A11 A12 A16 B11 B12 B16                                                         *
*     A12 A22 A26 B12 B22 B26        ! A  B !                                         *
*     A16 A26 A66 B16 B26 B66    or  !      !                                         *
*     B11 B12 B16 D11 D12 D16        ! B  D !                                         *
*     B12 B22 B26 D12 D22 D26                                                         *
*     B16 B26 B66 D16 D26 D66                                                         *
*                                                                                     *
* The stress vector is (Nx,Ny,Nxy,Mx,My,Mxy)                                          *
*                                                                                     *
* The deformation vector is (EPs xx, EPs yy, gam xy, Kx, Ky, Kxy)                     *
*                                                                                     *
* ----------------------------------------------------------------------------------- *
* SAMPLE RUN:                                                                         *
*                                                                                     *
*  CALCULATE THE STIFFNESS MATRIX OF A LAMINATE                                       *
*                                                                                     *
*  Number of layers: 4                                                                *
*   Layer 1: angle  -30.0 deg. - thickness    1.0 mm                                  *
*   Layer 2: angle   15.0 deg. - thickness    1.5 mm                                  *
*   Layer 3: angle  -15.0 deg. - thickness    1.5 mm                                  *
*   Layer 4: angle   30.0 deg. - thickness    1.0 mm                                  *
*                                                                                     *
*  Material Parameters:                                                               *
*     El =   38.0    Et =    9.0                                                      *
*   NUlt =   0.32    Glt=    3.6                                                      *
*                                                                                     *
*                                                                                     *
*  BASIC REDUCED STIFFNESS MATRIX IN GPa:                                             *
*    38.9445     2.9516     0.0000                                                    *
*     2.9516     9.2237     0.0000                                                    *
*     0.0000     0.0000     3.6000                                                    *
*                                                                                     *
*  Layer #1                                                                           *
*  REDUCED STIFFNESS MATRIX IN GPa:                                                   *
*    26.2896     8.1763    -9.4512                                                    *
*     8.1763    11.4292    -3.4183                                                    *
*    -9.4512    -3.4183     8.8247                                                    *
*                                                                                     *
*  Layer #2                                                                           *
*  REDUCED STIFFNESS MATRIX IN GPa:                                                   *
*    35.2120     4.6931     6.7316                                                    *
*     4.6931     9.4731     0.6986                                                    *
*     6.7316     0.6986     5.3416                                                    *
*                                                                                     *
*  Layer #3                                                                           *
*  REDUCED STIFFNESS MATRIX IN GPa:                                                   *
*    35.2120     4.6931    -6.7316                                                    *
*     4.6931     9.4731    -0.6986                                                    *
*    -6.7316    -0.6986     5.3416                                                    *
*                                                                                     *
*  Layer #4                                                                           *
*  REDUCED STIFFNESS MATRIX IN GPa:                                                   *
*    26.2896     8.1763     9.4512                                                    *
*     8.1763    11.4292     3.4183                                                    *
*     9.4512     3.4183     8.8247                                                    *
*                                                                                     *
*  A MATRIX IN 1E6 N/M:                                                               *
*   158.2153    30.4320     0.0000                                                    *
*    30.4320    51.2776     0.0000                                                    *
*     0.0000     0.0000    33.6741                                                    *
*                                                                                     *
*  B MATRIX IN 1E3 N:                                                                 *
*     0.0000     0.0000    22.6588                                                    *
*     0.0000     0.0000    12.1012                                                    *
*    22.6588    12.1012     0.0000                                                    *
*                                                                                     *
*  D MATRIX IN N.DM:                                                                  *
*   293.9255    77.3325     0.0000                                                    *
*    77.3325   114.6529     0.0000                                                    *
*     0.0000     0.0000    84.0869                                                    *
*                                                                                     *
*  STIFFNESS R MATRIX:                                                                *
*  1.5822e+08  3.0432e+07  0.00000000  0.00000000  0.00000000  2.2659e+04             *
*  3.0432e+07  5.1278e+07  0.00000000  0.00000000  0.00000000  1.2101e+04             *
*  0.00000000  0.00000000  3.3674e+07  2.2659e+04  1.2101e+04  0.00000000             *
*  0.00000000  0.00000000  2.2659e+04  293.925540  77.3325243  0.00000000             *
*  0.00000000  0.00000000  1.2101e+04  77.3325243  114.652882  0.00000000             *
*  2.2659e+04  1.2101e+04  0.00000000  0.00000000  0.00000000  84.0868611             *
*                                                                                     *
*  END OF PROGRAM.                                                                    *
* ----------------------------------------------------------------------------------- *
* REFERENCE: "MATERIAUX COMPOSITES - Comportement mécanique et analyse des structures *
*             By J.-M. Berthelot - MASSON 1996".                                      *
*                                                                                     *
*                                        Pascal Release By J-P MOREAU - August 1997.  *
*                                                    (www.jpmoreau.fr)                *
**************************************************************************************}
PROGRAM COMPOS03;
Uses WinCrt;


{ MAXCOU = maximum number of layers }
    Const NMAX=3; MAXCOU=10;

    Type
      MAT=Array[1..NMAX,1..NMAX] of Double;
      MAT1=Array[1..MAXCOU,1..NMAX,1..NMAX] of Double;
      MAT6=Array[1..6,1..6] of Double;
      VEC=Array[1..MAXCOU] of Double;
    Var 
      EL,ET,HM,NULT,GLT: Double;
{ layer angle and thickness } 
      TH, EP, H: VEC;
      Q0, A, B, D: MAT;
      Q: MAT1;
      R: MAT6;
      I,J,K, NCOUCHES: Integer;

      PROCEDURE disp_real(l2 : double); Forward;

{ set a matrix (N,N) to zero }
    Procedure MATNUL(N:Integer; A:MAT);   
    Var I,J:Integer;
    Begin 
      For I:=1 to N do
        for J:=1 to N do
          A[I,J]:=0.0
    End;               

{ A(3,3) = B(3,3) }
    Procedure MATEGAL(K:Integer; Var A:MAT1; B:MAT);
    Var I,J:Integer;
    Begin
      for I:=1 to NMAX do
        for J:=1 to NMAX do
          A[K,I,J]:=B[I,J]
    End;
     
{ multiply a 3x3 matrixpar by a 3x1 vector }  
    Procedure MATMULT(N:Integer; A:MAT; B:VEC; Var C:VEC);
    var SUM:Double; I,K:Integer;
    Begin
      for I:=1 to N do
      begin
        SUM:= 0.0;
        for K:=1 to N do
          SUM:=SUM+A[I,K]*B[K];
        C[I]:=SUM
      end
    End;
  
{ print a NxN matrix }           
    Procedure MATPR2(title:String; N:Integer; A:MAT);
    Var I,J:Integer; 
    Begin
      writeln;
      writeln(title);
      for I:=1 to N do
      begin
        for j:=1 to N do write(' ',A[I,J]:10:4);
        writeln
      end
    End;

    Procedure MATPR1(title:String; K,N:Integer; A:MAT1);
    Var I,J:Integer; 
    Begin
      writeln;
      writeln(title);
      for I:=1 to N do
      begin
        for j:=1 to N do write(' ',A[K,I,J]:10:4);
        writeln
      end
    End;

    Procedure MATPR6(title:String; N:Integer; A:MAT6);
    Var I,J:Integer; 
    Begin
      writeln;
      writeln(title);
      for I:=1 to N do
      begin
        for j:=1 to N do disp_real(A[I,J]);
        writeln
      end
    End;

    FUNCTION dcos4(t:Double): Double;
    Var a: Double;
    Begin
      a:=cos(t);
      dcos4 := a*a*a*a
    End;

    FUNCTION dsin4(t:Double): Double;
    Var a: Double;
    Begin
      a:=sin(t);
      dsin4 := a*a*a*a
    End;
  
  
{  calculate the reduced stiffness matrix for theta=0 }
    Procedure CAL_Q0(Var Q:MAT);
    Begin   
      Q[1,1] := EL/(1.0-((ET/EL)*NULT*NULT));
      Q[2,2] := (ET/EL)*Q[1,1];
      Q[1,2] := NULT*Q[2,2];
      Q[2,1] := Q[1,2];
      Q[3,3] := GLT
    End;  
  
{  calculate the reduced stiffness Q matrix for theta <> 0 
   with the condition that the Q matrix for theta=0 is available
   called here Q2                      }
    Procedure CAL_Q(K:Integer; Q2:MAT; Var Q:MAT1; TH:Double);  
    Var TEMP:Double;
    Begin
      TEMP := SIN(TH)*SIN(TH)*COS(TH)*COS(TH);
      Q[K,1,1]:= Q2[1,1]*DCOS4(TH)+Q2[2,2]*DSIN4(TH)+2.0*(Q2[1,2]+2.0*Q2[3,3])*TEMP;
      Q[K,1,2]:= (Q2[1,1]+Q2[2,2]-4.0*Q2[3,3])*TEMP+Q2[1,2]*(DCOS4(TH)+DSIN4(TH));
      Q[K,2,1]:=Q[K,1,2];
      TEMP := SIN(TH)*COS(TH)*COS(TH)*COS(TH);
      Q[K,1,3]:=(Q2[1,1]-Q2[1,2]-2.0*Q2[3,3])*TEMP;
      TEMP := SIN(TH)*SIN(TH)*SIN(TH)*COS(TH);
      Q[K,1,3]:=Q[K,1,3]+(Q2[1,2]-Q2[2,2]+2.0*Q2[3,3])*TEMP;
      Q[K,3,1]:=Q[K,1,3];
      TEMP := SIN(TH)*SIN(TH)*COS(TH)*COS(TH);
      Q[K,2,2]:= Q2[1,1]*DSIN4(TH)+2.0*(Q2[1,2]+2.0*Q2[3,3])*TEMP+Q2[2,2]*DCOS4(TH);
      TEMP := SIN(TH)*SIN(TH)*SIN(TH)*COS(TH);
      Q[K,2,3]:=(Q2[1,1]-Q2[1,2]-2.0*Q2[3,3])*TEMP;
      TEMP := COS(TH)*COS(TH)*COS(TH)*SIN(TH);
      Q[K,2,3]:=Q[K,2,3]+(Q2[1,2]-Q2[2,2]+2.0*Q2[3,3])*TEMP;
      Q[K,3,2]:=Q[K,2,3];
      TEMP := SIN(TH)*SIN(TH)*COS(TH)*COS(TH);
      Q[K,3,3]:= (Q2[1,1]+Q2[2,2]-2.0*(Q2[1,2]+Q2[3,3]))*TEMP+Q2[3,3]*(DSIN4(TH)+DCOS4(TH))
    End;

  FUNCTION evalue_log (l : double) : INTEGER;
  {This procedure allows calculating n, the role of which is explained in
   procedure Disp_real}
  BEGIN
    IF l<1 THEN evalue_log:=trunc(ln(l)/ln(10))-1
    ELSE evalue_log:=trunc(ln(l)/ln(10)*1.0000001);
    {Multiplying by 1.000 .. 0001 allows to avoid that the Pascal language,
    calculating by default, give such a result as log(100)=1.9999999
    leading to 1 pour n, which would be erroneous}
  END;

  FUNCTION power10 (n : INTEGER) : double;
   {This simple procedure iteratively calculates a power of 10 to avoid
   the rounds up of a formula such as exp(ln(10)*n)                   }
   VAR i : INTEGER;
  temp : double;
   BEGIN
    temp:=1;
    IF n>=0 THEN FOR i:=1 TO n DO temp:=temp*10
    ELSE FOR i:=1 TO -n DO temp:=temp*0.1;
    power10:=temp;
  END;

  PROCEDURE disp_real(l2 : double);
 {--------------------------------------------------------------------------
  This procedure allows correctly displaying to screen a real number with 12 
  {characters. The pascal language does not provide an automatic display of
  big (or small) numbers in scientific notation (with exponant) which can
  cause some discomfort.
  The scope is to display a real number entirely with 12 characters beginning
  with a space and by usind if necessary the scientific notation with four
  decimals.
  --------------------------------------------------------------------------}

  {l equals abs(l2) that allows to have positive numbers}

  VAR  l: double;

  {texte_res is the 12 characters string that will be displayed to screen}
    texte_res : STRING[12];

  {nom is a temporary variable used when  converting reals -> strings}
    nom : STRING;

  {n is used for the scientific notation}
    n : INTEGER;

  BEGIN  {disp_real}

    {if the l2 coefficient equals zero, special log. treatment}

    IF l2=0 THEN texte_res:='  0.00000000'

    ELSE
    BEGIN

      {else we take the absolute value, sign will be examined later.}
      l:=abs(l2);

      {initialize n and texte_res}
      texte_res:='';
      n:=0;

      IF (l>=10000) OR (l<0.001) THEN

      {if l is very big or very small, use the scientific notation by
      calculating and multiplying l by power10(-n), this allows 
      getting a number between 1 and 10}
      BEGIN
        n:=evalue_log(l);
        l:=l*power10(-n);
      END

      {else, round up to avoid inconvenience, such as 1.99999}
      ELSE IF (l<>0) THEN IF (abs((l-round(l))/l)<0.0000001) THEN l:=round(l);
      
      IF n=0 THEN
      {In this case, display the number in a normal way, since n=0 means no
      use of scientific notation}
      BEGIN
        {code first integer part of number}
        str(trunc(l),nom);
        texte_res:=nom;
        {Then, the decimal part avoiding that 0.00123 becomes 0.123}
        texte_res:=texte_res+'.';
        str(round((l-int(l))*power10(10-length(texte_res))),nom);
        IF (10-length(texte_res))>length(nom) THEN
        WHILE NOT ((10-length(texte_res))=length(nom)) DO nom:='0'+nom;
        texte_res:=texte_res+nom;
      END
      ELSE
      BEGIN
        {First code the integer part with 1 character}
        str(trunc(l),nom);
        texte_res:=nom;
        texte_res:=texte_res+'.';
        {Then code with 4 characters the decimal part avoiding: 0.01 -> 0.1}
        str(round((l-int(l))*10000),nom);
        IF length(nom)<4 THEN WHILE length(nom)<4 DO nom:='0'+nom;
        texte_res:=texte_res+nom;
        {At last set the scientific notation using classical e}
        IF n<0 THEN texte_res:=texte_res+'e-' else texte_res:=texte_res+'e+';
        str(abs(n),nom);
        IF length(nom)=1 THEN nom:='0'+nom;
        {Two characters are reserved to code the power of 10, which is
        enough to cover the real type}

        texte_res:=texte_res+nom;

      END;

      {if the number was negative, add a minus sign and pad the string up 
      to 12 characters}

      IF l2<0 THEN texte_res:='-'+texte_res;
      WHILE NOT (length(texte_res)>=12) DO texte_res:=' '+texte_res;

    END;

    WRITE(texte_res);

  END;  {disp_real}


{main program}
BEGIN
    
{ material constants in GPa }

      EL := 38.0;
      ET :=  9.0;
      NULT :=  0.32;
      GLT := 3.6;

{ material data }

      writeln;
      writeln('  CALCULATE THE STIFFNESS MATRIX OF A LAMINATE');
      writeln;
{ number of layers }    
      NCOUCHES:=4;
{ define angles and thicknesses of layers 1-4 }      
      TH[1]:=-30.0;
      EP[1]:=1E-3;
      TH[2]:= 15.0;
      EP[2]:=1.5E-3;
      TH[3]:=-15.0;
      EP[3]:=1.5E-3;
      TH[4]:= 30.0;
      EP[4]:=1E-3;
           
      writeln('  Number of layers: ', NCOUCHES);
      for I:=1 to NCOUCHES do
      begin
        writeln('   Layer ',I,': angle ',TH[I]:6:1,' deg. - thickness ',1000*EP[I]:6:1,' mm');
        {put angles in radians}
        TH[I] := TH[I]*PI/180.0
      end; 

      writeln;
      writeln('  Material Parameters:');
      writeln('     El = ',EL:6:1,'    Et = ',ET:6:1);
      writeln('   NUlt = ',NULT:6:2,'    Glt= ',GLT:6:1);
      writeln;

{ reduced stiffness matrix (angle = 0) }
      MATNUL(NMAX,Q0); 
      CAL_Q0(Q0);  
      MATPR2 ('  BASIC REDUCED STIFFNESS MATRIX IN GPa:',NMAX,Q0);

      ReadKey;

{ calculate reduced stiffness matrices for different angles <> 0 }
      for I:=1 to NCOUCHES do
      begin
        IF (TH[I] = 0.0) THEN MATEGAL(I,Q,Q0)
        ELSE  CAL_Q(I,Q0,Q,TH[I]);
        writeln;
        write('  Layer #',I);
        MATPR1('  REDUCED STIFFNESS MATRIX IN GPa:',I,NMAX,Q)
      end;  

      ReadKey;
      Clrscr;
      writeln;

{ calculate A matrix
  A = sum of 1e3*EP(k)*Q(K)  }
      MATNUL(NMAX,A);
      for I:=1 to NMAX do
        for J:=1 to NMAX do
        begin
          for K:=1 to NCOUCHES do
            A[I,J]:=A[I,J]+1000*(EP[K]*Q[K,I,J]);
          IF A[I,J] < 1E-8 Then A[I,J]:=0.0
        end;  
      MATPR2('  A MATRIX IN 1E6 N/M:',NMAX,A);

{ calculate B matrix
  B = sum of 1e6*(1/2)*(h(k)*h(k)-h(k-1)*h(k-1))*Q(k) }
      MATNUL(NMAX,B);
      HM:=0.0;
      for I:=1 to NCOUCHES do HM:=HM+EP[I];
{ average height of stack }
      HM:=HM/2.0;
{ lower height layer #1 }
      H[1] := -HM;
      for I:=2 to NCOUCHES+1 do H[I]:=H[I-1]+EP[I-1];
      for I:=1 to NMAX do
        for J:=1 to NMAX do
        begin
          for K:=1 to NCOUCHES do
            B[I,J]:=B[I,J]+1E6*0.5*((H[K+1]*H[K+1]-H[K]*H[K])*Q[K,I,J]);
         IF B[I,J] < 1E-8 Then B[I,J]:=0.0
        end;
      MATPR2('  B MATRIX IN 1E3 N:',NMAX,B);                           
    
{ calculate D matrix
  D = sum of 1e9*(1/3)*(h(k)*h(k)*h(k)-h(k-1)*h(k-1)*h(k-1))*Q(k) }
      MATNUL(NMAX,D);
      for I:=1 to NMAX do
        for J:=1 to NMAX do
        begin
          for K:=2 to NCOUCHES+1 do
            D[I,J]:=D[I,J]+1E9*(1.0/3.0)*((H[K]*H[K]*H[K]-H[K-1]*H[K-1]*H[K-1])*Q[K-1,I,J]);
          IF D[I,J] < 1E-8 Then D[I,J]:=0.0
        end;  
      MATPR2('  D MATRIX IN N.DM:',NMAX,D);
    
{ assemble the stiffness matrix of the laminate 
           ! A B !
       R = !     !
           ! B D !           }
      for I:=1 to 6 do
        for J:=1 to 6 do
          R[I,J]:=0.0;
      for I:=1 to NMAX do
        for J:=1 to NMAX do
          R[I,J] := 1E6*A[I,J];
      for I:=1 to NMAX do
        for J:=4 to NMAX+3 do
          R[I,J] := 1E3*B[I,J-3];
      for I:=4 to NMAX+3 do
        for J:=1 to NMAX do
          R[I,J] := 1E3*B[I-3,J];
      for I:=4 to NMAX+3 do
        for J:=4 to NMAX+3 do
          R[I,J] := D[I-3,J-3];
      MATPR6('  STIFFNESS R MATRIX:',6,R);
      writeln;
      writeln('  END OF PROGRAM.');
      ReadKey;
      DoneWinCrt

END.

{ end of file compos03.pas}