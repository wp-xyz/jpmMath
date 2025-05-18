  {****************************************************************
  * Utility procedures to read/write real numbers from/to screen  *
  * or text files.                                                *
  ****************************************************************}
  UNIT Utilit;

  INTERFACE

  USES WinCrt1, WinTypes, WinProcs;

  TYPE      AR_REEL = Double;

  FUNCTION  Read_integer(min:INTEGER;max:longint):INTEGER;
  PROCEDURE Read_filename(VAR fich:TEXT);
  PROCEDURE disp_real(l2 : ar_reel);
  PROCEDURE f_disp_real(VAR f:TEXT; l2 : ar_reel);
  PROCEDURE WriteHead(VAR f:TEXT; s:STRING);
  PROCEDURE WriteEnd(VAR f:TEXT);

  IMPLEMENTATION
  
  CONST Separator =
  '-------------------------------------------------------------------------';

  FUNCTION evalue_log (l : ar_reel) : INTEGER;
  {This procedure allows to calculate n, the role of which is explained in
   procedure Disp_real}
  BEGIN
    IF l<1 THEN evalue_log:=trunc(ln(l)/ln(10))-1
    ELSE evalue_log:=trunc(ln(l)/ln(10)*1.0000001);
    {Multiplying by 1.000 .. 0001 allows to avoid that the Pascal language,
    calculating by default, give such a result as log(100)=1.9999999
    leading to 1 pour n, which would be erroneous}
  END;

  FUNCTION power10 (n : INTEGER) : ar_reel;
   {This simple procedure iteratively calculates a power of 10 to avoid
   the rounds up of a formula such as exp(ln(10)*n)                   }
   VAR i : INTEGER;
  temp : ar_reel;
   BEGIN
    temp:=1;
    IF n>=0 THEN FOR i:=1 TO n DO temp:=temp*10
    ELSE FOR i:=1 TO -n DO temp:=temp*0.1;
    power10:=temp;
  END;

  {read an integer number controlling errors}
  FUNCTION Read_integer;
  VAR markx,marky:BYTE; loop:INTEGER;
  num:LONGINT;
  BEGIN
    REPEAT
      markx:=whereX; marky:=whereY;
      {$I-}
      READ(num);
      {$I+}
      loop:=IOResult;
      IF (loop<>0) OR (num<min) OR (num>max) THEN
      BEGIN
        GotoXY(markx,marky); ClrEol; MessageBeep(0)
      END;
    UNTIL (loop=0) AND (num>=min) AND (num<=max);
    Read_integer:=num
  END;

  {read a data input file name controlling existence}
  PROCEDURE Read_filename;
  VAR nomfich:STRING; ok:BOOLEAN;
  BEGIN
    REPEAT
      ok:=TRUE;
      gotoxy(3,2); Clreol;
      WRITE('Data filename: '); READLN(nomfich);
      {$I-}
      Assign(fich,nomfich); Reset(fich);
      {$I+}
      IF IOResult<>0 THEN
      BEGIN
        ok:=FALSE;
        MessageBeep(0);
        Writeln(' Error: File not found !')
      END
    UNTIL ok=TRUE
  END;

  PROCEDURE disp_real(l2 : ar_reel);
 {--------------------------------------------------------------------------
  This procedure allows to correctly display to screen a real number with 10 
  {characters. The pascal language does not provide an automatic display of
  big (or small) numbers in scientific notation (with exponant) which can
  cause some discomfort.
  The scope is to display a real number entirely with 10 characters beginning
  with a space and by usind if necessary the scientific notation.
  --------------------------------------------------------------------------}

  {l equals abs(l2) that allows to have positive numbers}

  VAR  l: ar_reel;

  {texte_res is the 10 characters string that will be displayed to screen}
    texte_res : STRING[10];

  {nom is a temporary variable used when  converting reals -> strings}
    nom : STRING;

  {n is used for the scientific notation}
    n : INTEGER;

  BEGIN  {disp_real}

    {if the l2 coefficient equals zero, special log. treatment}

    IF l2=0 THEN texte_res:='  0.000000'

    ELSE
    BEGIN

      {else we take the absolute value, sign will be examined later.}
      l:=abs(l2);

      {initialize n and texte_res}
      texte_res:='';
      n:=0;

      IF (l>=10000) OR (l<0.001) THEN

      {if l is very big or very small, use the scientific notation by
      calculating and multiplying l by power10(-n), this allows to 
      get a nuber between 1 and 10}
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
        str(round((l-int(l))*power10(8-length(texte_res))),nom);
        IF (8-length(texte_res))>length(nom) THEN
        WHILE NOT ((8-length(texte_res))=length(nom)) DO nom:='0'+nom;
        texte_res:=texte_res+nom;
      END
      ELSE
      BEGIN
        {First code the integer part with 1 character}
        str(trunc(l),nom);
        texte_res:=nom;
        texte_res:=texte_res+'.';
        {Then code with 2 characters the decimal part avoiding: 0.01 -> 0.1}
        str(round((l-int(l))*100),nom);
        IF length(nom)<2 THEN WHILE length(nom)<2 DO nom:='0'+nom;
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
      to 10 characters}

      IF l2<0 THEN texte_res:='-'+texte_res;
      WHILE NOT (length(texte_res)>=10) DO texte_res:=' '+texte_res;

    END;

    WRITE(texte_res);

  END;  {disp_real}

  {same as disp_real but to a text file}
  PROCEDURE f_disp_real(VAR f:TEXT; l2 : ar_reel);
  VAR
  l: ar_reel;
  texte_res: STRING[10];
  nom: STRING;
  n  : INTEGER;

  BEGIN  {f_disp_real}

    IF l2=0 THEN texte_res:='  0.000000'

    ELSE
    BEGIN

      l:=abs(l2);
      texte_res:='';
      n:=0;

      IF (l>=10000) OR (l<0.001) THEN
      BEGIN
        n:=evalue_log(l);
        l:=l*power10(-n);
      END
      ELSE IF (l<>0) THEN IF (abs((l-round(l))/l)<0.0000001) THEN l:=round(l);

      IF n=0 THEN
      BEGIN
        str(trunc(l),nom);
        texte_res:=nom;
        texte_res:=texte_res+'.';
        str(round((l-int(l))*power10(8-length(texte_res))),nom);
        IF (8-length(texte_res))>length(nom) THEN
        WHILE NOT ((8-length(texte_res))=length(nom)) DO nom:='0'+nom;
        texte_res:=texte_res+nom;
      END
      ELSE
      BEGIN
        str(trunc(l),nom);
        texte_res:=nom;
        texte_res:=texte_res+'.';
        str(round((l-int(l))*100),nom);
        IF length(nom)<2 THEN WHILE length(nom)<2 DO nom:='0'+nom;
        texte_res:=texte_res+nom;
        IF n<0 THEN texte_res:=texte_res+'e-' else texte_res:=texte_res+'e+';
        str(abs(n),nom);
        IF length(nom)=1 THEN nom:='0'+nom;
        texte_res:=texte_res+nom;
      END;

      IF l2<0 THEN texte_res:='-'+texte_res;
      WHILE NOT (length(texte_res)>=10) DO texte_res:=' '+texte_res;

    END;

    WRITE(f,texte_res);

  END;  {f_disp_real }

  PROCEDURE WriteHead;  {write a caption to output file}
  BEGIN
    Writeln(f,Separator);
    Writeln(f,s);
    Writeln(f,Separator)
  END;

  PROCEDURE WriteEnd;  {write a separation line in output file}
  BEGIN
    Writeln(f,Separator)
  END;

  END.

{end of file utilit.pas (english edition) }
