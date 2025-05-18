{**************************************************************************************
* This program calculates the deformations and stresses in a laminated material made  *
* of n unidirectional composite layers, knowing the resulting imposed efforts.        *
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
* The imposed efforts vector is (Nx,Ny,Nxy,Mx,My,Mxy)                                 *
*                                                                                     *
* The deformation vector is (EPs xx, EPs yy, gam xy, Kx, Ky, Kxy)                     *
*                                                                                     *
* ----------------------------------------------------------------------------------- *
* SAMPLE RUN:                                                                         *
*   The output file Compos04.lst contains (extract):                                  *
*                                                                                     *
*  CALCULATE DEFORMATIONS AND STRESSES                                                *
*       IN A LAMINATED MATERIAL                                                       *
*                                                                                     *
*  Imposed resulting efforts:                                                         *
*     NX =   1.0000e+06  NY =   5.0000e+05  NXY =   2.5000e+05                        *
*     MX =   0.00000000  MY =   0.00000000  MXY =   0.00000000                        *
*                                                                                     *
*  Number of layers: 4                                                                *
*                                                                                     *
*   Layer 1: angle   15.0 deg. - thickness    1.5 mm                                  *
*                                                                                     *
*   Layer 2: angle  -30.0 deg. - thickness    1.0 mm                                  *
*                                                                                     *
*   Layer 3: angle  -15.0 deg. - thickness    1.5 mm                                  *
*                                                                                     *
*   Layer 4: angle   30.0 deg. - thickness    1.0 mm                                  *
*                                                                                     *
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
*  Layer #:1                                                                          *
*                                                                                     *
*  REDUCED STIFFNESS MATRIX IN GPa:                                                   *
*    35.2120     4.6931     6.7316                                                    *
*     4.6931     9.4731     0.6986                                                    *
*     6.7316     0.6986     5.3416                                                    *
*                                                                                     *
*  Layer #:2                                                                          *
*                                                                                     *
*  REDUCED STIFFNESS MATRIX IN GPa:                                                   *
*    26.2896     8.1763    -9.4512                                                    *
*     8.1763    11.4292    -3.4183                                                    *
*    -9.4512    -3.4183     8.8247                                                    *
*                                                                                     *
*  Layer #:3                                                                          *
*                                                                                     *
*  REDUCED STIFFNESS MATRIX IN GPa:                                                   *
*    35.2120     4.6931    -6.7316                                                    *
*     4.6931     9.4731    -0.6986                                                    *
*    -6.7316    -0.6986     5.3416                                                    *
*  --/--                                                                              *
*  Stresses in each layer in main axes:                                               *
*                                                                                     *
*  Layer #1                                                                           *
*   For z =     -2.50000 mm:                                                          *
*    sigma L :    282.73562  MPa                                                      *
*    sigma T :     78.07995  MPa                                                      *
*    sigma LT:     45.15936  MPa                                                      *
*   For z =     -1.00000 mm:                                                          *
*    sigma L :    288.58230  MPa                                                      *
*    sigma T :     70.53304  MPa                                                      *
*    sigma LT:     34.70107  MPa                                                      *
*                                                                                     *
*  Layer #2                                                                           *
*   For z =     -1.00000 mm:                                                          *
*    sigma L :     86.43422  MPa                                                      *
*    sigma T :    105.75930  MPa                                                      *
*    sigma LT:      5.73649  MPa                                                      *
*   For z =      0.00000 mm:                                                          *
*    sigma L :    111.92641  MPa                                                      *
*    sigma T :     96.96499  MPa                                                      *
*    sigma LT:      8.38896  MPa                                                      *
*                                                                                     *
*  Layer #3                                                                           *
*   For z =      0.00000 mm:                                                          *
*    sigma L :    151.46586  MPa                                                      *
*    sigma T :     90.07486  MPa                                                      *
*    sigma LT:     21.12949  MPa                                                      *
*   For z =      1.50000 mm:                                                          *
*    sigma L :    192.64459  MPa                                                      *
*    sigma T :     76.37100  MPa                                                      *
*    sigma LT:     19.34600  MPa                                                      *
*                                                                                     *
*  Layer #4                                                                           *
*   For z =      1.50000 mm:                                                          *
*    sigma L :    333.21160  MPa                                                      *
*    sigma T :     51.87584  MPa                                                      *
*    sigma LT:      8.77294  MPa                                                      *
*   For z =      2.50000 mm:                                                          *
*    sigma L :    317.90586  MPa                                                      *
*    sigma T :     50.19097  MPa                                                      *
*    sigma LT:      1.40860  MPa                                                      *
*                                                                                     *
* ----------------------------------------------------------------------------------- *
* REFERENCE: "MATERIAUX COMPOSITES - Comportement mécanique et analyse des structures *
*             By J.-M. Berthelot - MASSON 1996".                                      *
*                                                                                     *
*                                        Pascal Release By J-P MOREAU - August 1997.  *
*                                                    (www.jpmoreau.fr)                *
**************************************************************************************}
PROGRAM COMPOS04;
Uses WinCrt;

{ MAXCOU = maximum number of layers }
Const NMAX=3; MAXCOU=10;

Type
      MAT = Array[1..NMAX,1..NMAX] of Double;
      MAT1 = Array[1..MAXCOU,1..NMAX,1..NMAX] of Double;
      MAT2 = Array[1..MAXCOU,1..NMAX] of Double;
      MAT6 = Array[1..6,1..6] of Double;
      VEC = Array[1..MAXCOU] of Double;
      VEC1 = Array[1..NMAX] of Double;
      VEC6 = Array[1..6] of Double;
Var
      EL,ET,HM,NULT,GLT: Double;

{ angle and thickness of layers } 
      TH, EP, H: VEC;

      A,B,D, Q0, T: MAT;  Q: MAT1;
      R, X:MAT6; N, E0, EDEB, EFIN: VEC6;
      DET,NX,NY,NXY,MX,MY,MXY: Double;
      AK, AK1, BK, BK1: MAT2;

      TEMP, TEMP1: VEC1;
      I,J,K, NCOUCHES: Integer;

      fp:TEXT;

      PROCEDURE disp_real(l2 : double); Forward;

{ set NxN matrix to zero }
    Procedure MATNUL(A:MAT; N:Integer);   
    Var I,J: Integer;
    Begin 
      for I:=1 to N do
        for J:=1 to N do
          A[I,J]:=0.0
    End;               

{ 3x3 A matrix = 3x3 B matrix }
    Procedure MATEGAL(K,N:Integer; Var A:MAT1; B:MAT);   
    Var I,J: Integer;
    Begin
      for I:=1 to N do
        for J:=1 to N do
          A[K,I,J]:=B[I,J]
    End;
     
{ multiply a NxN matrix by a Nx1 vector }  
    Procedure MATMUL(A:MAT; B:VEC1; Var C:VEC1; N:Integer);
    Var I,K:Integer; SUM:Double;
    Begin 
      for I:=1 to N do
      begin
        SUM:= 0.0;
        for K:=1 to N do SUM:=SUM+A[I,K]*B[K];
        C[I]:=SUM
      end                                                                   
    End;

    Procedure MATMUL6(A:MAT6; B:VEC6; Var C:VEC6; N:Integer);
    Var I,K:Integer; SUM:Double;
    Begin 
      for I:=1 to N do
      begin
        SUM:= 0.0;
        for K:=1 to N do SUM:=SUM+A[I,K]*B[K];
        C[I]:=SUM
      end                                                                   
    End;

{ print a NxN matrix with caption }           
    Procedure MATPR2(title:String; N:Integer; A:MAT);
    Var I,J:Integer; 
    Begin
      writeln(fp);
      writeln(fp,title);
      for I:=1 to N do
      begin
        for j:=1 to N do write(fp,' ',A[I,J]:10:4);
        writeln(fp)
      end
    End;

    Procedure MATPR1(title:String; K,N:Integer; A:MAT1);
    Var I,J:Integer; 
    Begin
      writeln(fp);
      writeln(fp,title);
      for I:=1 to N do
      begin
        for j:=1 to N do write(fp,' ',A[K,I,J]:10:4);
        writeln(fp)
      end
    End;

    Procedure MATPR6(title:String; N:Integer; A:MAT6);
    Var I,J:Integer; 
    Begin
      writeln(fp);
      writeln(fp,title);
      for I:=1 to N do
      begin
        for j:=1 to N do disp_real(A[I,J]);
        writeln(fp)
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
Procedure MATINV(N,M:integer; VAR AA:MAT6; VAR BB:MAT6; VAR DET:DOUBLE);
CONST EPSMACH = 1E-15;
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
          writeln(fp,'  The determinant equals ZERO !!!');                                                              
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
          writeln(fp,'  PIVOT TOO SMALL - STOP');                               
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


{ Matrix to change axes (x,y) for deformations }

    Procedure INITT(Var T:MAT; TH:Double);
    Begin
      T[1,1]:=COS(TH)*COS(TH);
      T[1,2]:=SIN(TH)*SIN(TH);
      T[1,3]:=SIN(TH)*COS(TH);
      T[2,1]:=SIN(TH)*SIN(TH);
      T[2,2]:=COS(TH)*COS(TH);
      T[2,3]:=-SIN(TH)*COS(TH);
      T[3,1]:=-2.0*SIN(TH)*COS(TH);
      T[3,2]:=2.0*SIN(TH)*COS(TH);
      T[3,3]:=COS(TH)*COS(TH)-SIN(TH)*SIN(TH)
    End;

{ Matrix to change axes (x,y) for stresses }

    Procedure INITT1(Var T:MAT; TH:Double);
    Begin
      T[1,1]:=COS(TH)*COS(TH);
      T[1,2]:=SIN(TH)*SIN(TH);
      T[1,3]:=2.0*SIN(TH)*COS(TH);
      T[2,1]:=SIN(TH)*SIN(TH);
      T[2,2]:=COS(TH)*COS(TH);
      T[2,3]:=-2.0*SIN(TH)*COS(TH);
      T[3,1]:=-SIN(TH)*COS(TH);
      T[3,2]:=SIN(TH)*COS(TH);
      T[3,3]:=COS(TH)*COS(TH)-SIN(TH)*SIN(TH)
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

    WRITE(fp, texte_res);

  END;  {disp_real}


{main program}
BEGIN
    
{ Material constants in GPa }

      EL := 38.0;
      ET :=  9.0;
      NULT :=  0.32;
      GLT := 3.6;

{ imposed efforts }
      NX:=1E6;
      NY:=0.5E6;
      NXY:=0.25E6;
      MX:=0.0;
      MY:=0.0;
      MXY:=0.0;

      N[1]:=1E-9*NX;
      N[2]:=1E-9*NY;
      N[3]:=1E-9*NXY;
      N[4]:=1E-9*MX;
      N[5]:=1E-9*MY;
      N[6]:=1E-9*MXY;

      for I:=1 to 6 do
      begin
        E0[I]:=0.0;
        EDEB[I]:=0.0;
        EFIN[I]:=0.0
      end;

{ print data to output file }
      Assign(fp,'compos04.lst'); Rewrite(fp);
      WRITELN(fp);
      WRITELN(fp,'  CALCULATE DEFORMATIONS AND STRESSES');
      WRITELN(fp,'       IN A LAMINATED MATERIAL');    
      WRITELN(fp); 
      WRITELN(fp,'  Imposed resulting efforts:');
      write(fp,'     NX = '); disp_real(NX);
      write(fp,'  NY = '); disp_real(NY);
      write(fp,'  NXY = '); disp_real(NXY); writeln(fp);
      write(fp,'     MX = '); disp_real(MX);
      write(fp,'  MY = '); disp_real(MY);
      write(fp,'  MXY = '); disp_real(MXY); writeln(fp);

{ number of layers in laminate }    
      NCOUCHES:=4;

{ define angles and thicknesses of layers 1-4 }      
      TH[1]:=15.0;
      EP[1]:=1.5E-3;
      TH[2]:= -30.0;
      EP[2]:=1E-3;
      TH[3]:=-15.0;
      EP[3]:=1.5E-3;
      TH[4]:= 30.0;
      EP[4]:=1E-3;
           
      WRITELN(fp);
      WRITELN(fp,'  Number of layers: ', NCOUCHES);
      WRITELN(fp);
      for I:=1 to NCOUCHES do
      begin
        writeln(fp,'   Layer ',I,': angle ',TH[I]:6:1,' deg. - thickness ',1000*EP[I]:6:1,' mm');
        writeln(fp);
        {put angles in radians}
        TH[I] := TH[I]*PI/180.0
      end;

      writeln(fp);
      writeln(fp,'  Material Parameters:');
      writeln(fp,'     El = ',EL:6:1,'    Et = ',ET:6:1);
      writeln(fp,'   NUlt = ',NULT:6:2,'    Glt= ',GLT:6:1);
      writeln(fp);
           
{ reduced stiffness matrix (angle = 0)  }
      MATNUL(Q0,NMAX); 
      CAL_Q0(Q0);  
      MATPR2 ('  BASIC REDUCED STIFFNESS MATRIX IN GPa:', NMAX, Q0);
    
{ calculate the reduced stiffness matrices for different angles <> 0 }
      for I:=1 to NCOUCHES do
      begin
        IF TH[I] = 0.0 THEN 
          MATEGAL(I,NMAX,Q,Q0)
        ELSE 
          CAL_Q(I,Q0,Q,TH[I]);                 
        WRITELN(fp);
        WRITELN(fp,'  Layer #:',I);
        MATPR1('  REDUCED STIFFNESS MATRIX IN GPa:', I, NMAX, Q)
      end;
     
{ calculate A matrix
  A = sum of 1e3*EP(k)*Q(K) }
      MATNUL(A,NMAX);
      for I:=1 to NMAX do
        for J:=1 to NMAX do
        begin
          for K:=1 to NCOUCHES do
            A[I,J]:=A[I,J]+1000*(EP[K]*Q[K,I,J]);
          IF ABS(A[I,J]) < 1E-8 Then A[I,J]:=0.0
        end;  
      MATPR2('  A MATRIX IN 1E6 N/M:', NMAX, A);
                      
{ calculate B matrix
  B = sum of 1e6*(1/2)*(h(k)*h(k)-h(k-1)*h(k-1))*Q(k) }
      MATNUL(B,NMAX);
      HM:=0.0;
      for I:=1 to NCOUCHES do HM:=HM+EP[I];
{ average height of stack }
      HM:=HM/2.0;
{ lower height of layer #1 }
      H[1] := -HM;
      for I:=2 to NCOUCHES+1 do H[I]:=H[I-1]+EP[I-1];
      for I:=1 to NMAX do
        for J:=1 to NMAX do
        begin
          for K:=1 to NCOUCHES do
            B[I,J]:=B[I,J]+1E6*0.5*((H[K+1]*H[K+1]-H[K]*H[K])*Q[K,I,J]);
         IF ABS(B[I,J]) < 1E-8 Then B[I,J]:=0.0
        end;
      MATPR2('  B MATRIX IN 1E3 N:',NMAX,B);
    
{ calculate D matrix
  D = sum of 1e9*(1/3)*(h(k)*h(k)*h(k)-h(k-1)*h(k-1)*h(k-1))*Q(k) }
      MATNUL(D,NMAX);
      for I:=1 to NMAX do
        for J:=1 to NMAX do
        begin
          for K:=2 to NCOUCHES+1 do
            D[I,J]:=D[I,J]+1E9*(1.0/3.0)*((H[K]*H[K]*H[K]-H[K-1]*H[K-1]*H[K-1])*Q[K-1,I,J]);
          IF ABS(D[I,J]) < 1E-8 Then D[I,J]:=0.0
        end;  
      MATPR2('  D MATRIX IN N.DM:',NMAX,D);
    
{ assemble the stiffness matrix of laminate 
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
      MATPR6('  STIFFNESS R MATRIX:', 6, R);

{ Direct inversion of R matrix, stored in R }
      MATINV(6,0,R,X,DET);
      WRITELN(fp);
      WRITE(fp,'  DETERMINANT OF R: '); disp_real(DET); writeln(fp);
      for I:=1 to 6 do
        for J:=1 to 6 do
          R[I,J]:=1E9*R[I,J];
      MATPR6('  INVERSE MATRIX OF R:', 6, R);

{ Calculate déformations eps0xx, eps0yy, gam0xy,kx,ky,kxy
  in main axes, stored in matrix E0(6)   }

      MATMUL6(R,N,E0,6);
      WRITELN(fp);
      WRITELN(fp,'  Membrane deformations (in mm) and curbatures:');
      WRITELN(fp);
      writeln(fp,'    eps0xx: ', 1E3*E0[1]:12:5);
      writeln(fp,'    eps0yy: ', 1E3*E0[2]:12:5);
      writeln(fp,'    eps0xy: ', 1E3*E0[3]:12:5);
      writeln(fp,'    kx    : ', E0[4]:12:5);
      writeln(fp,'    ky    : ', E0[5]:12:5);
      writeln(fp,'    kxy   : ', E0[6]:12:5);

{ Deformations in reference axes (x,y)
      ! epsxx !   ! eps0xx !    ! kx  !
      ! epsyy ! = ! eps0yy ! +z ! ky  !
      ! gamxy !   ! gam0xy !    ! kxy !    }

      WRITELN(fp);
      WRITELN(fp,'  Deformations in reference axes (x,y):');
      writeln(fp,'   For z = ',1E3*H[1]:12:5,' mm:');
      for I:=1 to NMAX do EDEB[I]:=E0[I]+H[1]*E0[I+3];
      writeln(fp,'    epsxx: ',1E3*EDEB[1]:12:5,' mm');
      writeln(fp,'    epsyy: ',1E3*EDEB[2]:12:5,' mm');
      writeln(fp,'    epsxy: ',1E3*EDEB[3]:12:5,' mm');
      writeln(fp,'   For z = ',1E3*H[NCOUCHES+1]:12:5,' mm:');
      for I:=1 to NMAX do EFIN[I]:=E0[I]+H[NCOUCHES+1]*E0[I+3];
      writeln(fp,'    epsxx: ',1E3*EFIN[1]:12:5,' mm');
      writeln(fp,'    epsyy: ',1E3*EFIN[2]:12:5,' mm');
      writeln(fp,'    epsxy: ',1E3*EFIN[3]:12:5,' mm');

{ Deformations in each layer in main axes of layer }
      WRITELN(fp);
      WRITELN(fp,'  Deformations in each layer in main axes:');
      for K:=1 to NCOUCHES do
      begin
        INITT(T,TH[K]);
        for I:=1 to NMAX do TEMP[I]:=E0[I];
        MATMUL(T,TEMP,TEMP1,3);
        for I:=1 to NMAX do AK[K,I]:=TEMP1[I];
        for I:=1 to NMAX do TEMP[I]:=E0[I+3];
        MATMUL(T,TEMP,TEMP1,3);
        for I:=1 to NMAX do BK[K,I]:=TEMP1[I];
        WRITELN(fp);
        WRITELN(fp,'  Layer #', K);
        writeln(fp,'   For z = ',1E3*H[K]:12:5,' mm:');
        for I:=1 to NMAX do EDEB[I]:=AK[K,I]+H[K]*BK[K,I];
        writeln(fp,'    epsxx: ',1E3*EDEB[1]:12:5,' mm');
        writeln(fp,'    epsyy: ',1E3*EDEB[2]:12:5,' mm');
        writeln(fp,'    epsxy: ',1E3*EDEB[3]:12:5,' mm');
        writeln(fp,'   For z = ',1E3*H[K+1]:12:5,' mm:');
        for I:=1 to NMAX do EFIN[I]:=AK[K,I]+H[K+1]*BK[K,I];
        writeln(fp,'    epsxx: ',1E3*EFIN[1]:12:5,' mm');
        writeln(fp,'    epsyy: ',1E3*EFIN[2]:12:5,' mm');
        writeln(fp,'    epsxy: ',1E3*EFIN[3]:12:5,' mm')
      end;

{ Stresses in each layer in reference axes }

      WRITELN(fp);
      WRITELN(fp,'  Stresses in each layer in reference axes:');
      for K:=1 to NCOUCHES do
      begin
        for I:=1 to NMAX do TEMP[I]:=E0[I];
        for I:=1 to NMAX do
          for J:=1 to NMAX do
            A[I,J]:=Q[K,I,J];
        MATMUL(A,TEMP,TEMP1,3);
        for I:=1 to NMAX do AK[K,I]:=TEMP1[I];
        for I:=1 to NMAX do TEMP[I]:=E0[I+3];
        MATMUL(A,TEMP,TEMP1,3);
        for I:=1 to NMAX do BK[K,I]:=TEMP1[I];
        WRITELN(fp);
        WRITELN(fp,'  Layer #', K);
        writeln(fp,'   For z = ',1E3*H[K]:12:5,' mm:');
        for I:=1 to NMAX do EDEB[I]:=1E3*(AK[K,I]+H[K]*BK[K,I]);
        writeln(fp,'    sigma L : ',EDEB[1]:12:5,'  MPa');
        writeln(fp,'    sigma T : ',EDEB[2]:12:5,'  MPa');
        writeln(fp,'    sigma LT: ',EDEB[3]:12:5,'  MPa');
        writeln(fp,'   For z = ',1E3*H[K+1]:12:5,' mm:');
        for I:=1 to NMAX do EFIN[I]:=1E3*(AK[K,I]+H[K+1]*BK[K,I]);
        writeln(fp,'    sigma L : ',EFIN[1]:12:5,'  MPa');
        writeln(fp,'    sigma T : ',EFIN[2]:12:5,'  MPa');
        writeln(fp,'    sigma LT: ',EFIN[3]:12:5,'  MPa')
      end;

{ Stresses in each layer in main axes }

      WRITELN(fp);
      WRITELN(fp,'  Stresses in each layer in main axes:');
      for K:=1 to NCOUCHES do
      begin
        INITT1(T,TH[K]);
        for I:=1 to NMAX do TEMP[I]:=AK[K,I];
        MATMUL(T,TEMP,TEMP1,3);
        for I:=1 to NMAX do AK1[K,I]:=TEMP1[I];
        for I:=1 to NMAX do TEMP[I]:=BK[K,I];
        MATMUL(T,TEMP,TEMP1,3);
        for I:=1 to NMAX do BK1[K,I]:=TEMP1[I];
        WRITELN(fp);
        WRITELN(fp,'  Layer #', K);
        writeln(fp,'   For z = ',1E3*H[K]:12:5,' mm:');
        for I:=1 to NMAX do EDEB[I]:=1E3*(AK1[K,I]+H[K]*BK1[K,I]);
        writeln(fp,'    sigma L : ',EDEB[1]:12:5,'  MPa');
        writeln(fp,'    sigma T : ',EDEB[2]:12:5,'  MPa');
        writeln(fp,'    sigma LT: ',EDEB[3]:12:5,'  MPa');
        writeln(fp,'   For z = ',1E3*H[K+1]:12:5,' mm:');
        for I:=1 to NMAX do EFIN[I]:=1E3*(AK1[K,I]+H[K+1]*BK1[K,I]);
        writeln(fp,'    sigma L : ',EFIN[1]:12:5,'  MPa');
        writeln(fp,'    sigma T : ',EFIN[2]:12:5,'  MPa');
        writeln(fp,'    sigma LT: ',EFIN[3]:12:5,'  MPa')
      end;
     
      WRITELN(fp);
      Close(fp);
      writeln('  END OF PROGRAM.');

      ReadKey;
      DoneWinCrt

END.  {of main program}

{ End of file Compos04.pas}