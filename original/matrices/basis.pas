 {**********************************************************
  *    UNIT basis.pas - Pascal version of basis_r.c,        *
  *    reduced and modified version of basis.c              *
  *                                  By J-P Moreau, Paris   *
  *                                    (www.jpmoreau.fr)    *
  * ------------------------------------------------------- *
  *    Collection of basic routines used by programmes      *
  *    concerning matrices.                                 *                          
  * ------------------------------------------------------- *
  * Reference of original basis.c:                          *
  *                                                         *
  * "Numerical Algorithms with C, By Gisela Engeln-Muellges *
  *  and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].    *
  **********************************************************}
  UNIT basis;     

  INTERFACE
  Uses WinCrt;

  CONST
  SIZE      = 2047;
  MACH_EPS  = 1E-12;
  Separator =
  '--------------------------------------------------------------------';

  TYPE
        pVECT  = ^VECT;
        VECT   = Array[0..SIZE] of REAL;
        pIVECT = ^IVECT;
        IVECT  = Array[0..SIZE] of INTEGER;

  VAR   fp : TEXT;

  Function  norm_max(VAR vektor : pVECT; n : INTEGER): REAL;
  Procedure copy_vector(VAR ziel : pVECT; VAR quelle : pVECT; n : INTEGER);
  Procedure LogError(text:STRING);
  Procedure LogError1(VAR fp:TEXT; text:STRING);
  Procedure ReadVec (VAR fp:TEXT; n, ic:INTEGER; VAR x:pVECT);
  Procedure ReadVec1(VAR fp:TEXT; n, ic:INTEGER; VAR x:pVECT);
  Procedure SetVec  (n:INTEGER; VAR x:pVECT; val:REAL);
  Procedure WriteVec(VAR fp:TEXT; n, ic:INTEGER; VAR x:pVECT);
  Procedure WriteVec1(VAR fp:TEXT; n, ic:INTEGER; VAR x:pVECT);
  Procedure ReadMat (VAR fp:TEXT; m, n, ic:INTEGER; VAR a:pVECT);
  Procedure ReadMat1(VAR fp:TEXT; m, n, ic:INTEGER; VAR a:pVECT);
  Procedure WriteMat(VAR fp:TEXT; m, n, ic:INTEGER; VAR a:pVECT);
  Procedure WriteMat1(VAR fp:TEXT; m, n, ic:INTEGER; VAR a:pVECT);
  Procedure SetMat  (m, n:INTEGER; VAR a:pVECT; val:REAL);
  Procedure WriteHead (VAR fp:TEXT; nom : string);
  Procedure WriteEnd  (VAR fp:TEXT);
  Function  min(a,b:INTEGER):INTEGER;
  Function  max(a,b:INTEGER):INTEGER;
  Procedure Swap(VAR a:REAL; VAR b:REAL);
  PROCEDURE f_aff_reel(VAR fp:TEXT; l2 : REAL);

  IMPLEMENTATION     

  {print a message to srdout in case of error}
  Procedure LogError(text:STRING);
  begin
    Writeln(' ERROR - ',text);
    Readln;
    Close(fp);
    DoneWinCrt
  end;

  {print a message to output file in case of error}
  Procedure LogError1(VAR fp:TEXT; text:STRING);
  begin
    Writeln(fp,' ERROR - ',text);
  end;

  FUNCTION evalue_log (l : REAL) : INTEGER;

  {This procedure is called by f_aff_reel}

  BEGIN

    IF l<1 THEN evalue_log:=trunc(ln(l)/ln(10))-1
    ELSE evalue_log:=trunc(ln(l)/ln(10)*1.0000001);

    {Multiplication by 1.000 .. 0001 avoid to get a result such as
     log(100)=1.9999999}
    
  END;

  FUNCTION puissance10 (n : INTEGER) : REAL;

  {This simple procedure calculates iteratively a number power 10 to
   avoid round results such as in formula exp(ln(10)*n).}

  VAR i : INTEGER;
  temp : REAL;

  BEGIN

    temp:=1;

    IF n>=0 THEN FOR i:=1 TO n DO temp:=temp*10
    ELSE FOR i:=1 TO -n DO temp:=temp*0.1;

    puissance10:=temp;

  END;

  Function norm_max    { Find the maximum norm of a REAL vector .........}
             (
              VAR vektor:pVECT;            { vector ................. }
              n:INTEGER                    { length of vector ....... }
             ) : REAL;                     { Maximum norm ........... }

  {***********************************************************************
  * Return the maximum norm of a [0..n-1] vector  v.                     *
  *                                                                      *
  * global names used:                                                   *
  * ==================                                                   *
  * REAL, ABS, ZERO                                                    *
  ***********************************************************************}
  VAR
    norm : REAL;                                             { local max }
    betrag : REAL;                            { magnitude of a component }
    i : INTEGER;
  Begin
    norm:=0.0;
    for i:=0 to n-1 do
    begin
      betrag:=abs(vektor^[i]);
      if betrag > norm then norm := betrag
    end;   
    norm_max := norm
  End;
  
  Procedure copy_vector     { copy a REAL vector ........................}
                (
                 VAR ziel:pVECT;           { copied vector ............}
                 VAR quelle:pVECT;         { original vector ..........}
                 n:INTEGER                 { length of vector .........}
                );
  {***********************************************************************
  * copy the n elements of the vector quelle into the vector ziel.       *
  *                                                                      *
  * global names used:                                                   *
  * ==================                                                   *
  * REAL                                                               *
  ***********************************************************************}
  Var i : INTEGER;
  Begin
    for i:=0 to n-1 do
      ziel^[i]:=quelle^[i]
  End;  

  Procedure ReadVec (VAR fp: TEXT; n, ic: INTEGER; VAR x:pVECT);
 {====================================================================*
 *                                                                    *
 * Read a vector of length x from input file fp. Index starts at ZERO.*
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *                                                                    *
 *      n        lenght of vector of INTEGER type.                    *
 *      ic       number of items per line (INTEGER).                  *
 *      x        pointer to vector of pVECT type.                     *
 *                                                                    *
 *====================================================================}
  Var
    i : INTEGER;
  Begin
    for i := 0 to n-1 do
    begin
      read(fp,x^[i]);
      if (i+1 MOD ic) = 0 then readln(fp)
    end
  End;

  Procedure ReadVec1 (VAR fp: TEXT; n, ic: INTEGER; VAR x:pVECT);
 {====================================================================*
 *                                                                    *
 * Read a vector of length x from input file fp. Index starts at ONE. *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *                                                                    *
 *      n        lenght of vector of INTEGER type.                    *
 *      ic       number of items per line (INTEGER).                  *
 *      x        pointer to vector of pVECT type.                     *
 *                                                                    *
 *====================================================================}
  Var
    i : INTEGER;
  Begin
    for i := 1 to n do
    begin
      read(fp,x^[i]);
      if (i+1 MOD ic) = 0 then readln(fp)
    end
  End;
   
  Procedure SetVec;
 {====================================================================*
 *                                                                    *
 *  Initialize a vector of length n with a constant value val.        *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *      n        INTEGER ( n > 0 )                                    *
 *               length of vector                                     *
 *      a        pointer to vector of pVECT type.                     *
 *      val      constant value                                       *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      a        vector with constant value val in every position     *
 *                                                                    *
 *====================================================================}
  Var
    i:INTEGER;
  Begin
    for i := 0 to n-1 do
      x^[i] := val
  End;     
     
Procedure WriteVec;
 {====================================================================*
 *                                                                    *
 * Put out vector of length x to  output file. Index starts at ZERO.  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *                                                                    *
 *      n        lenght of vector of INTEGER type.                    *
 *      ic       number of items per line (INTEGER).                  *
 *      x        pointer to vector of pVECT type.                     *
 *                                                                    *
 *====================================================================}
  Var
    i,compte : INTEGER;
  Begin
    compte:=0;
    for i := 0 to n-1 do
    begin
      Inc(compte);
      f_aff_reel(fp,x^[i]);
      if compte = ic then
      begin
        writeln(fp);
        compte:=0
      end
    end
  End;

  Procedure WriteVec1;
 {====================================================================*
 *                                                                    *
 * Put out vector of length x to  output file. Index starts at ONE.   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *                                                                    *
 *      n        lenght of vector of INTEGER type.                    *
 *      ic       number of items per line (INTEGER).                  *
 *      x        pointer to vector of pVECT type.                     *
 *                                                                    *
 *=====================================================================
    ic = number of items per line }
  Var
    i,compte : INTEGER;
  Begin
    compte:=0;
    for i := 1 to n do
    begin
      Inc(compte);
      f_aff_reel(fp,x^[i]);
      if compte = ic then
      begin
        writeln(fp);
        compte:=0
      end
    end
  End;
     
Procedure ReadMat (VAR fp:TEXT; m, n, ic:INTEGER; VAR a:pVECT);
 {====================================================================*
 *                                                                    *
 *  Read an m x n matrix from input file. Index starts at ZERO.       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters :                                                *
 *  ==================                                                *
 *      m        int m; ( m > 0 )                                     *
 *               number of rows of matrix                             *
 *      n        int n; ( n > 0 )                                     *
 *               column number of  matrix                             *
 *      a        pointer of pVECT type.                               *
 *               matrix in vector of VECT type                        *
 *      ic       number of items per line                             *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      a        matrix                                               *
 *                                                                    *
 *   ATTENTION : WE do not allocate storage for a here.               *
 *                                                                    *
 *====================================================================}
Begin
  ReadVec (fp, m*n, ic, a);
End;

Procedure ReadMat1 (VAR fp:TEXT; m, n, ic:INTEGER; VAR a:pVECT);
 {====================================================================*
 *                                                                    *
 *  Read an m x n matrix from input file. Index starts at ONE.        *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters :                                                *
 *  ==================                                                *
 *      m        int m; ( m > 0 )                                     *
 *               number of rows of matrix                             *
 *      n        int n; ( n > 0 )                                     *
 *               column number of  matrix                             *
 *      a        pointer of pVECT type.                               *
 *               matrix in vector of VECT type                        *
 *      ic       number of items per line                             *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      a        pointer to matrix in vector of VECT type             *
 *                                                                    *
 *   ATTENTION : WE do not allocate storage for a here.               *
 *                                                                    *
 *====================================================================}
Var tmp : pVECT;
    i,j : INTEGER;
Begin
  New(tmp);
  for i:=1 to m do
  begin
    ReadVec1(fp, n, ic, tmp);
    for j:=1 to n do
      a^[i*(n+1)+j]:=tmp^[j]
  end;
  Dispose(tmp)
End;

Procedure WriteMat;
 {====================================================================*
 *                                                                    *
 *  Put out an m x n matrix in output file. Index starts at ZERO.     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *      m        int m; ( m > 0 )                                     *
 *               row number of matrix                                 *
 *      n        int n; ( n > 0 )                                     *
 *               column number of matrix                              *
 *      a        pointer to vector of pVECT type.                     *
 *               matrix                                               *
 *====================================================================}
Var i:INTEGER;
Begin
  WriteVec (fp, m*n, ic, a);
End;

Procedure WriteMat1;
 {====================================================================*
 *                                                                    *
 *  Put out an m x n matrix in output file. Index starts at ONE.      *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *      m        int m; ( m > 0 )                                     *
 *               row number of matrix                                 *
 *      n        int n; ( n > 0 )                                     *
 *               column number of matrix                              *
 *      a        pointer to vector of pVECT type.                     *
 *               matrix                                               *
 *====================================================================}
Var i,j:INTEGER;
    tmp:pVECT;
Begin
  New(tmp);
  for i:=1 to m do
  begin
    for j:=1 to n do tmp^[j]:=a^[i*(n+1)+j]; 
    WriteVec1(fp, n, ic, tmp);
  end;
  Dispose(tmp)
End;
         
Procedure SetMat(m, n:INTEGER; VAR a:pVECT; val:REAL);
 {====================================================================*
 *                                                                    *
 *  Initialize an m x n matrix with a constant value val .            *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *      m        int m; ( m > 0 )                                     *
 *               row number of matrix                                 *
 *      n        int n; ( n > 0 )                                     *
 *               column number of matrix                              *
 *      a        pointer to vector of pVECT type.                     *
 *               matrix                                               *
 *      val      constant value                                       *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      a        matrix with constant value val in every position     *
 *                                                                    *
 *====================================================================}
Begin
  SetVec(m*n, a, val)
End;

Procedure WriteHead (VAR fp:TEXT; nom : string);
 {====================================================================*
 *                                                                    *
 *  Put out header with text in string in file fp.                    *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *      nom      of type STRING.                                      *
 *               text of headertext                                   *
 *                                                                    *
 *====================================================================}
Begin
  Writeln(fp,Separator);
  Writeln(fp,nom);
  Writeln(fp,Separator)
End;

Procedure WriteEnd  (VAR fp:TEXT);
 {====================================================================*
 *                                                                    *
 *  Put out a separation line onto output file fp.                    *
 *                                                                    *
 *====================================================================}
Begin
  Writeln(fp,Separator)
End;

Function min;
Begin
  if a < b then min:=a else min:=b
End;

Function max;
Begin
  if a > b then max:=a else max:=b
End;

Procedure Swap;
Var temp : REAL;
Begin
  temp:=a;
  a:=b;
  b:=temp
End;

  PROCEDURE f_aff_reel(VAR fp:TEXT; l2 : REAL);

  {This procedure displays in a correct way with 10 characters any real number. 
   The pascal language does not put big numbers automatically into scientific
   notation (like C language for instance) which may cause some unpleasantness.
   The real scope is to display any real number on 10 characters, beginning by
   an empty space and using the scientific, if necessary. This version writes 
   to a TEXT file. }                            }

  VAR  {l equals abs(l2) allowing to work only on positive numbers}

  l : REAL;

  {texte_res is the 10 characters string displayed to screen}

  texte_res : STRING[10];

  {nom is a temporary variable used during real conversions}
  
  nom : STRING;

  {n equals the number of figures -1 of integer part if >1, else equals
   the number of zéros +1 at begin of its decimal part; this allows to
   build the scientific notation}

  n : INTEGER;


  BEGIN

    {if l2 coefficient equals zero,  special teatment (log problem) }

    IF l2=0 THEN texte_res:='  0.000000'

    ELSE
    BEGIN

      {else we work on its absolute value, sign will be taken into account later}

      l:=abs(l2);

      {n and texte_res initialization}

      texte_res:='';
      n:=0;

      IF (l>=10000) OR (l<0.001) THEN

      {if l is very big or very small,  we use scientific notation by
       multiplying log(l) by  10 power -n to get a number between 1 and 10}

      BEGIN

        n:=evalue_log(l);
        l:=l*puissance10(-n);

      END

      {else, the number is rounded off to avoid results such as 1.99999}

      ELSE IF (l<>0) THEN IF (abs((l-round(l))/l)<0.0000001) THEN l:=round(l);


      IF n=0 THEN

      {In this case, the number is displayed normally, n=0 means that scientific
       notation is not in use}

      BEGIN

        {The string first contains the integer part of the number}

        str(trunc(l),nom);
        texte_res:=nom;

        {then the decimal part}

        texte_res:=texte_res+'.';
        str(round((l-int(l))*puissance10(8-length(texte_res))),nom);
        IF (8-length(texte_res))>length(nom) THEN
        WHILE NOT ((8-length(texte_res))=length(nom)) DO nom:='0'+nom;
        texte_res:=texte_res+nom;

      END

      ELSE
      BEGIN

        {We begin to code the integer part that will take 1 character}

        str(trunc(l),nom);
        texte_res:=nom;
        texte_res:=texte_res+'.';

        {then the decimal part on two characters}

        str(round((l-int(l))*100),nom);
        IF length(nom)<2 THEN WHILE length(nom)<2 DO nom:='0'+nom;
        texte_res:=texte_res+nom;

        {Finally, the scientific notation is built with e symbol}

        IF n<0 THEN texte_res:=texte_res+'e-' else texte_res:=texte_res+'e+';
        str(abs(n),nom);
        IF length(nom)=1 THEN nom:='0'+nom;

        {Twox characters are used for the power of 10}

        texte_res:=texte_res+nom;

      END;

      {if the number was negative, one adds up a minus sign padds the string
       up to 10 characters}

      IF l2<0 THEN texte_res:='-'+texte_res;

      WHILE NOT (length(texte_res)>=10) DO texte_res:=' '+texte_res;

    END;

    WRITE(fp,texte_res);

  END;  {f_aff_reel}

END.

{End of file basis.pas}