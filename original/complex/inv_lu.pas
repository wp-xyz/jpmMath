{*********************************************************************************
*          Inversion of a complex square matrix by LU decomposition              *
*          with dynamic allocations                                              *
*                                                                                *
*                                            Pascal version by J-P Moreau, Paris *
*                                                    (www.jpmoreau.fr)           *
* ------------------------------------------------------------------------------ *
* Uses:  units Wincrt, Basis.pas, Lu.pas                                         *
*                                                                                *
* SAMPLE RUN:                                                                    *
*                                                                                *  
* Input file (inv_clu.dat):                                                      *
*                                                                                *
* 4                                                                              *
*  8  0  2    0  3  0  12    1                                                   *
*  2  0  4    0  7  0   0.25 2                                                   *
*  3  0  7    0  3  0   5    3                                                   *
* 12  0  0.25 0  5  0   2    4                                                   *
*                                                                                *
* Output file (inv_clu.lst):                                                     *
*                                                                                *
* INVERSION OF A COMPLEX SQUARE MATRIX:                                          *
*                                                                                *
* N=4                                                                            *
*                                                                                *
* 8.000000  0.000000  2.000000  0.000000  3.000000  0.000000  12.00000  1.000000 *
* 2.000000  0.000000  4.000000  0.000000  7.000000  0.000000  0.250000  2.000000 *
* 3.000000  0.000000  7.000000  0.000000  3.000000  0.000000  5.000000  3.000000 *
* 12.00000  0.000000  0.250000  0.000000  5.000000  0.000000  2.000000  4.000000 *
*                                                                                *
* Inverted matrix Y:                                                             *
*                                                                                *
*-0.030217 -0.041620 -0.084258 -0.006409  0.053057  0.014655  0.104257  0.025151 *
*-0.074213 -0.048474 -0.060236 -0.007464  0.198125  0.017069  0.009984  0.029293 *
* 0.054432  0.008897  0.201874  0.001370 -0.129568 -0.003133 -0.037542 -0.005377 *
* 0.104314  0.024909  0.016062  0.003835 -0.036731 -0.008771 -0.063037 -0.015052 *
*                                                                                *
* Verification A*Y = I:                                                          *
*                                                                                *
* 1.000000 -4.55e-13 -6.82e-13  0.000000 -4.55e-13  0.000000  9.09e-13  2.27e-13 *
*-3.98e-13 -2.27e-13  1.000000  0.000000 -2.52e-12  0.000000 -3.98e-13  0.000000 *
*-1.36e-12 -4.55e-13 -5.68e-13  0.000000  1.000000  2.27e-13  4.55e-13  4.55e-13 *
* 0.000000 -4.55e-13 -8.53e-13  0.000000  3.41e-13  0.000000  1.000000  0.000000 *
*                                                                                *
* End of file inv_clu.lst                                                        *
*                                                                                *
*********************************************************************************}
  Program Inversion_LU;
  Uses WinCrt, CLu;

  Var
  A : pCVEC;     { matrix 0:n x 0:n stored in a vector (see NOTA2) }
  A1: pCVec;     { copy of matrix A }
  Y : pCVEC;     { matrix 0:n x 0:n stored in a vector }
  temp : pCVEC;  { vector 0:n }
  INDX : pIVEC;  { integer vector 0:n }

  {NOTA1: index zero is not used here.
   NOTA2: The element i,j of matrix A(n,n) is the element i*(n+1)+j
          of the vector A, n being the 2nd dimension of A(n,n) }

  d, i, j, n, rc : integer;
  input, output, s : STRING;
  F1, F2 : TEXT;


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

  PROCEDURE f_aff_reel(VAR fp:TEXT; l2 : REAL);

  {----------------------------------------------------------------------------
   This procedure displays in a correct way with 10 characters any real number. 
   The pascal language does not put big numbers automatically into scientific
   notation (like C language for instance) which may cause some unpleasantness.
   The real scope is to display any real number on 10 characters, beginning by
   an empty space and using the scientific, if necessary. This version writes 
   to a TEXT file.
   ----------------------------------------------------------------------------}

  VAR   {l equals abs(l2) allowing to work only on positive numbers}

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

  Procedure ReadCVec (VAR fp: TEXT; n, ic: INTEGER; VAR x:pCVEC);
 {====================================================================*
 *                                                                    *
 * Read a complex vector of length x from input file fp. Index starts *
 * from ONE.                                                          *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *                                                                    *
 *      n        lenght of vector of INTEGER type.                    *
 *      ic       number of items per line (INTEGER).                  *
 *      x        pointer to vector of pCVEC type (See Unit CLU).      *
 *                                                                    *
 *====================================================================}
  Var
    i : INTEGER;
  Begin
    for i := 1 to n do
    begin
      read(fp,x^[i].R, x^[i].I);
      if (i+1 MOD ic) = 0 then readln(fp)
    end
  End;

  Procedure WriteCVec(VAR fp:TEXT; n, ic:INTEGER; VAR x:pCVEC);
 {====================================================================*
 *                                                                    *
 * Put out a complex vector of length x to  output file. Index starts *
 * from ONE.                                                          *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *                                                                    *
 *      n        lenght of vector of INTEGER type.                    *
 *      ic       number of items per line (INTEGER).                  *
 *      x        pointer to vector of pCVEC type.                     *
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
      f_aff_reel(fp,x^[i].R);
      f_aff_reel(fp,x^[i].I);
      if compte = ic then
      begin
        writeln(fp);
        compte:=0
      end
    end
  End;

  {******************************************                                     
  *      MULTIPLY TWO COMPLEX MATRICES      *
  * --------------------------------------- *                                     
  * INPUTS:    A  MATRIX N*N                *                                     
  *            B  MATRIX N*N                *                                     
  *            N  INTEGER                   *                                     
  * --------------------------------------- *                                     
  * OUTPUT:    C  MATRIX N*N, PRODUCT A*B   *                                     
  *                                         *                                     
  ******************************************}
  Procedure MatCMul(A,B:pCVec; VAR C:pCVec; n:Integer); 
  VAR SUM, tmp: Complex;
      I,J,K: INTEGER;
  Begin                                           
    for I:=1 to N do
    begin                                                                  
      for J:=1 to N do
      begin                                                                
        SUM.R:=0.0; SUM.I:=0.0;                                                                
        for K:=1 to N do
        begin
          CMUL(A^[I*(n+1)+K],B^[K*(n+1)+J],tmp);
          SUM.R:=SUM.R + tmp.R;
          SUM.I:=SUM.I + tmp.I
        end;                                                 
        C^[I*(n+1)+J]:=SUM                                                            
      end                                                                   
    end                                                                     
  End;


{main program}
BEGIN

  Writeln;
  Write(' Data file name (without .dat): '); read(s);
  input := s + '.dat';
  output := s + '.lst';

  Assign(F1,input); Reset(F1);
  Assign(F2,output); Rewrite(F2);

  readln(F1,n);     {size of matrix}

  New(A); New(A1); New(Y); New(temp); New(INDX);

  Writeln(F2);
  Writeln(F2,' INVERSION OF A COMPLEX SQUARE MATRIX:');
  Writeln(F2);
  writeln(F2,'  N=',n);
  writeln(F2);

  for i:=1 to n do
  begin
    ReadCVec(F1,n,n,temp);  {read a line}
    for j:=1 to n do
    begin
      A^[i*(n+1)+j] := temp^[j];
      A1^[i*(n+1)+j] := temp^[j];
      Y^[i*(n+1)+j].R := 0.0;
      Y^[i*(n+1)+j].I := 0.0
    end;
    Y^[i*(n+1)+i].R := 1.0;
    Y^[i*(n+1)+i].I := 0.0;
    WriteCVec(F2,n,n,temp)   {write a line}
  end;
  close(F1);

{call LU decomposition routine (only once) }
  LUDCMP(A,n,INDX,D,rc);

{call solver if previous return code is ok
 to obtain inverse of A one column at a time}
  if rc=0 then
    for j:=1 to n do
    begin
      for i:=1 to n do temp^[i]:=Y^[i*(n+1)+j];
      LUBKSB(A,n,INDX,temp);
      for i:=1 to n do Y^[i*(n+1)+j]:=temp^[i]
    end;
  {inverse of matrix A is now in matrix Y,
   matrix A is destroyed. }

  if rc=1 then
    writeln(F2,' The matrix is singular !')
  else
  begin
    writeln(F2);
    writeln(F2,'  Inverted matrix Y:');
    writeln(F2);  
    for i:=1 to n do
    begin
      for j:=1 to n do
      begin
        f_aff_reel(F2,Y^[i*(n+1)+j].R);
        f_aff_reel(F2,Y^[i*(n+1)+j].I)
      end;
      writeln(F2)
    end
  end;
  {verify A1 x Y = I (result put in A) }
  MatCMul(A1,Y,A,n);
  {A should now contain identity matrix}
  writeln(F2);
  writeln(F2,'  Verification A*Y = I:');
  writeln(F2);
  for i:=1 to n do
  begin
    for j:=1 to n do
    begin
      f_aff_reel(F2,A^[i*(n+1)+j].R);
      f_aff_reel(F2,A^[i*(n+1)+j].I)
    end;
    writeln(F2)
  end;
  Writeln(F2);
  Writeln(F2,' End of file ', output);
  Close(F2);
  Dispose(A); Dispose(A1); Dispose(Y); Dispose(temp); Dispose(INDX);
  Writeln;
  Writeln(' Results in file ',output,'.');
  Writeln;
  Readkey;
  DoneWinCrt

END.

{End of file test_lu.pas}