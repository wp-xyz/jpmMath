{***************************************************************
   * Recursive descent parser accepting variables beginning with  *
   * an upcase letter with possibility to save formulas.          *
   * ------------------------------------------------------------ *
   * Ref.: "Advanced Turbo C By Herbert Schildt, Borland-Osborne/ *
   *        McGraw-Hill, Berkeley, USA, 1987".                    *
   *                                                              * 
   *              TPW version 5, 01/07/1994 By J-P Moreau, Paris. *
   *                             (www.jpmoreau.fr)                *
   * ------------------------------------------------------------ *
   *         Examples of valid expressions:                       *
   *               A=10/4                                         *
   *               B=25.478                                       *
   *               C=-1                                           *
   *	           X=0.72                                         *
   *	           F=A*X^2+B*X+C                                  *
   *               F                                              *
   *               25.75*2+1.784                                  *
   *               2*pi*R                                         *
   * ------------------------------------------------------------ *
   *  Version 1, 12/10/1993 (basis version)                       * 
   *  Example 1:   F=A/B  )                                       *
   *               A=500  )  in any order                         *
   *               B=2    )                                       *
   *               F                                              *
   *  Result:      250.000000                                     *
   * ------------------------------------------------------------ *
   * Version 2, 13/12/1993: possibility to replace right hand va- *
   * riables by their formula, if already defined.                *
   * ------------------------------------------------------------ *
   *  Example 2 :  S=L*M                                          *
   *               V=S*H                                          *
   *               P=D*V                                          *
   *               L=100                                          *
   *               M=50                                           *
   *               H=5                                            *
   *               D=1.6                                          *
   *               P                                              *
   *  Result:      40000.000000                                   *
   * ------------------------------------------------------------ *
   * Version 3: added functions atan, cos, exp, log, root, sin,   *
   *            tan and constant pi.                              *
   * ------------------------------------------------------------ *
   *  Example 3:   exp(1) or e(1) or e1                           *
   *  Result:      2.718282                                       *
   *               sinus(2) or sin(2) or s(2) or s2               *
   *  Result:      0.909297                                       *
   *               root(225) or r(225) ou r225                    *
   *  Result:      15.000000                                      *
   * ------------------------------------------------------------ *
   * Version 4: Variables can have seveal characters,             *
   *            Intermediate values displayed,                    *
   *            Error messages in functions.                      *
   * ------------------------------------------------------------ *
   * Version 5: Possibility to read/save a session,               *
   *            added a main menu.                                *
   ***************************************************************}
PROGRAM Parser;
USES WinCrt,WinTypes,WinProcs,WObjects,Strings;

CONST
  DELIMITER=1;
  VARIABLE=2;
  NUMBER=3;
  NB_VAR=40;  {maximum number of variables}

TYPE
  tab_char = SET OF '!'..'z';
  variabl  = SET OF 'A'..'Z';
  S40 = STRING[40];
  p40 = ^S40;
  S80 = STRING[80];
  p80 = ^S80;

CONST
  FNC : tab_char = ['a','c','e','l','p','r','s','t'];
  FONC: tab_char = ['-','+','a','c','e','l','p','r','s','t'];
  FNCE: tab_char = ['a','c','e','g','i','l','o','n','p','r','s','t','u','x'];
  DELIM:tab_char = ['+','-','*','/','%','^','=','(',')'];
  DELI: tab_char = ['+','-','*','/','%','^','=','(',')','a','c','e','l','p','r','s','t'];
  VARI: variabl  = ['A'..'Z'];
  DIGIT:tab_char = ['-','+','.','0','1','2','3','4','5','6','7','8','9','e','E'];

VAR
  prog:STRING;    {the expression to evaluate}
  bak :STRING;    {store previous expression}
  token:S40;      {current element (operator, variable, number) }
  ii,k,tok_TYPE:integer;
  exist:boolean;

{table of formulas of kind A=... }
  vrbl   : ARRAY[0..NB_VAR] OF p40; {table of pointers to a variable}
  formule: ARRAY[0..NB_VAR] OF p80; {table of pointers to a formula}
  indice : integer;                 {formula identification number}
  reponse: double;                  {current result of calculation}
  flg    : ARRAY[0..NB_VAR] OF boolean;  { false=no  true=yes }
  ok     : boolean;
  num    : integer;
  choix  : char;
  f_IN, f_out : TEXT;
  nom_IN, nom_out : S40;

  {open application window with caption}
  PROCEDURE WinCrtInit(Nom:PChar);
  BEGIN
    WindowOrg.X:=25;
    WindowOrg.Y:=40;
    WindowSize.X:=590;
    WindowSize.Y:=400;
    StrCopy(WindowTitle,Nom);
    InitWinCrt;
  END;

  PROCEDURE get_token; forward;
  PROCEDURE serror(error:integer); forward;
  PROCEDURE level2(VAR resultat:double); forward;
  PROCEDURE level3(VAR resultat:double); forward;
  PROCEDURE level4(VAR resultat:double); forward;
  PROCEDURE level5(VAR resultat:double); forward;
  PROCEDURE level6(VAR resultat:double); forward;
  PROCEDURE arith(o:char;VAR r:double; h:double); forward;
  PROCEDURE unary(o:char;VAR r:double); forward;
  PROCEDURE primitive(VAR resultat:double); forward;
  FUNCTION  find_VAR: double; forward;
  FUNCTION  get_formula: boolean; forward;

  {entry point into parser}
  PROCEDURE calcule_exp(VAR resultat:double);
  BEGIN
    get_token;
    IF token='' THEN
    BEGIN
      serror(2);
      exit;
    END;
    level2(resultat);
  END;

  {add or substract two terms}
  PROCEDURE level2(VAR resultat:double);
  VAR op: char;
      hold:double;
  BEGIN
    level3(resultat);
    REPEAT
      op:=token[1];
      IF (op<>'+') AND (op<>'-') THEN exit;
      {write('  op=',op); readln;   }
      get_token;
      level3(hold);
      arith(op,resultat,hold);
    UNTIL (op<>'+') AND (op<>'-')
  END;

  {multiply or divise two terms}
  PROCEDURE level3(VAR resultat:double);
  VAR op: char;
      hold: double;
  BEGIN
    level4(resultat);
    REPEAT
      op:=token[1];
      IF (op<>'*') AND (op<>'/') AND (op<>'%') THEN exit;
      {write('  op=',op); readln; }
      get_token;
      level4(hold);
      arith(op,resultat,hold);
    UNTIL (op<>'*') AND (op<>'/') AND (op<>'%')
  END;

  {power integer or not}
  PROCEDURE level4(VAR resultat:double);
  VAR hold: double;
  BEGIN
    level5(resultat);
    IF(token='^') THEN
    BEGIN
      get_token;
      level4(hold);
      arith('^',resultat,hold)
    END
  END;

  {unary operator + , - or function}
  PROCEDURE level5(VAR resultat:double);
  VAR op: char;
  BEGIN
    IF ((tok_TYPE=DELIMITER) AND (token[1] IN FONC)) THEN
    BEGIN
      op:=token[1];
      get_token;
    END;
    { pas d'argument pour la constante pi }
    IF (op<>'p') THEN level6(resultat);
    IF (op IN FONC)  THEN unary(op,resultat);
  END;

  {expression with parentheses }
  PROCEDURE level6(VAR resultat:double);
  BEGIN
    IF ((token[1]='(') AND (tok_TYPE=DELIMITER)) THEN
    BEGIN
      get_token;
      level2(resultat);
      IF (token[1] <> ')') AND (ii < length(prog)+1) THEN serror(1);
      get_token;
    END
    ELSE primitive(resultat);
  END;

  {effective value of a number or a variable}
  PROCEDURE primitive(VAR resultat:double);
  VAR erreur: integer;
  BEGIN
    CASE (tok_TYPE) OF
      VARIABLE :
      BEGIN
	resultat:=find_VAR;
	get_token
      END;
      NUMBER :
      BEGIN
	Val(token,resultat,erreur);
	IF erreur <> 0 THEN resultat:=0.0;  
	get_token;
      END;
      ELSE serror(0);  {syntax error message}
    END
  END;

  FUNCTION Power(x,y:double): double;
  BEGIN
    IF x<0 THEN EXIT;
    Power:=Exp(x*Ln(y));
  END;

  FUNCTION Tan(x: double) : double;
  VAR Cx : double;
  BEGIN
    Cx:=Cos(x);
    IF (Cx=0) THEN
		WriteLn(' Argument of Tan is an odd multiple of Pi/2')
	      ELSE
	        Tan:=Sin(x)/Cx;
  END;

  {arithmetic  operations}
  PROCEDURE arith(o:char;VAR r:double; h:double);
  VAR t,ex: double;
  BEGIN
    CASE(o) OF
      '-': r:=r-h;
      '+': r:=r+h;
      '*': r:=r*h; 
      '/': IF h<>0 THEN r:=r/h
		   ELSE BEGIN
			  writeln('  Divide zero error!');
			  ok:=false
			END;
      '%': BEGIN t:=r/h; r:=r-(t*h); END;
      '^': r:=Power(h,r);
    END
  END;

  {process a unary operation}
  PROCEDURE unary(o:char;VAR r:double);
  BEGIN
    IF(o='-') THEN r:=-(r)
    ELSE IF (o='a') THEN r:=arctan(r)
    ELSE IF (o='c') THEN r:=cos(r)
    ELSE IF (o='e') THEN r:=exp(r)
    ELSE IF (o='l') THEN IF r>0 THEN r:=Ln(r)
		    ELSE BEGIN
			   writeln('  Argument of log negative or null !');
			   ok:=false
			END
    ELSE IF (o='p') THEN r:=PI
    ELSE IF (o='r') THEN IF (r>=0) THEN r:=sqrt(r)
  		                   ELSE BEGIN
		             	     writeln('  Argument of sqrt negative !');
				     ok:=false
				   END
    ELSE IF (o='s') THEN r:=sin(r)
    ELSE IF (o='t') THEN IF abs(r-PI/2)>1E-10 THEN r:=sin(r)/cos(r)
					      ELSE
                                              BEGIN
						writeln('  INFINITY');
						ok:=false
					      END
  END;                                        

{find value of a variable by using its formula
 if it exists, else return 0.   }
FUNCTION find_VAR: double;
VAR i,jj,n:integer;
    q:STRING;
    valeur:double;
BEGIN
  IF NOT(token[1] IN VARI) THEN
  BEGIN
    serror(1);
    find_VAR:=0;
    exit
  END;    
  {seek the identification number of the variable}
  n:=-1;
  FOR i:=0 TO num DO
    IF ((vrbl[i]<>NIL) AND (vrbl[i]^=token)) THEN n:=i;
  {case a new variable}
  IF n=-1 THEN
  BEGIN
    Inc(num);
    New(vrbl[num]);
    vrbl[num]^:=token;
    n:=num
  END;
  IF formule[n]<>NIL THEN   { such a formula exists}
  BEGIN
    q:=prog;                {current calculation saved in q}
    jj:=ii;                 {current position saved in j}
    prog:=formule[n]^;
    ii:=1;
    calcule_exp(valeur);    {calculate the variable with the formula}
    prog:=q;                {restore current calculation}
    ii:=jj;                 {and current position in expression}
    IF ((NOT flg[n]) AND (n>0)) THEN
      writeln('  ',vrbl[n]^,'=',valeur:16:8);
    flg[n]:=TRUE;
    find_VAR:=valeur;       {current value of the variable}
  END
  ELSE
  BEGIN
    IF ((NOT flg[n]) AND (n>0)) THEN
      writeln('  No formula found for ',vrbl[n]^);
    flg[n]:=TRUE;
    find_VAR:=0;            {no assigned value}
  END
END;

  {syntax error or other}
  PROCEDURE serror(error:integer);
  VAR e:ARRAY[0..2] OF STRING;
  BEGIN
    e[0]:='syntax error.';
    e[1]:='Odd number of parentheses.';
    e[2]:='No expression found.';
    writeln('  ',e[error]);
    ok:=FALSE
  END;

  {***********************************************
  *         Analyze a string's elements          *
  *      (delimiter, variable and number)        *
  * -------------------------------------------- *
  * Example : analyse the expression:            *
  *           Masse+100-(Base*S)/2               *
  * -------------------------------------------- *
  * The successive calls to get_token will give  *
  * the following result:                        *
  *        token         tok_TYPE                *
  *        -----         --------                *
  *        Masse         VARIABLE                *
  *          +           DELIMITER               *
  *         100           NUMBER                 *
  *          -           DELIMITER               *
  *          (           DELIMITER               *
  *         Base         VARIABLE                *
  *          *           DELIMITER               *
  *          S           VARIABLE                *
  *          )           DELIMITER               *
  *          /           DELIMITER               *
  *          2            NUMBER                 *
  *         NULL           NULL                  *
  ***********************************************}
  PROCEDURE get_token;
  VAR temp: S40;
  BEGIN
    tok_TYPE:=0;
    temp:='';
    IF (prog[ii] IN DELI) AND (ii<length(prog)+1) THEN
    BEGIN
      tok_TYPE:=DELIMITER;
      temp:=temp+prog[ii];
      {for functions we take only the first letter}
      IF prog[ii] IN FNC THEN
	WHILE prog[ii] IN FNCE DO inc(ii)
      ELSE inc(ii)
    END
    ELSE IF prog[ii] IN VARI THEN
    BEGIN
      WHILE (NOT (prog[ii] IN DELIM)) AND (ii<length(prog)+1) DO
      BEGIN
	temp:=temp+prog[ii];
	inc(ii)
      END;
      tok_TYPE:=VARIABLE;
    END
    ELSE IF prog[ii] IN DIGIT THEN 
    BEGIN
      WHILE (NOT (prog[ii] IN DELI)) AND (ii<length(prog)+1) DO
      BEGIN
	temp:=temp+prog[ii];
	inc(ii);
      END;
      tok_TYPE:=NUMBER;
    END;
    token:=temp;
    {writeln('  token=',token);   }
  END;
  
  FUNCTION get_formula: Boolean;
  {store a formula A...= à Z...= and return a flag FALSE
    (no value assigned) or TRUE (assignation OK) to the
    main program which prints no result if assignation OK}
  VAR i,r:integer;
      egal: boolean;
      nom : STRING;
  BEGIN
    egal:=FALSE;
    {seek a sign =} 
    FOR i:=2 TO length(prog) DO
      IF prog[i]='=' THEN
      BEGIN
	egal:=TRUE;
	r:=i;
	nom:=Copy(prog,1,r-1);
      END;
    IF egal THEN
    BEGIN       {sign = has been found}
      indice:=-1;
      {seek an existing variable}
      FOR i:=0 TO num DO
	IF vrbl[i]^=nom THEN indice:=i;
      {case a new variable}
      IF indice=-1 THEN
      BEGIN
	Inc(num);
	New(vrbl[num]);
	vrbl[num]^:=nom;
	indice:=num
      END;
      {copy right hand member in formule[indice]^ }
      IF formule[indice]=NIL THEN New(formule[indice]);
      IF formule[indice]<>NIL THEN
      BEGIN
        formule[indice]^:=Copy(prog,r+1,length(prog)); 
	get_formula:=TRUE;  { assignation Ok }
      END
      ELSE
	writeln('  Memory full !')
    END
    ELSE get_formula:=FALSE;   { pas d'assignation }
  END; 

  {menu option calculate}
  PROCEDURE exec_calcul;
  VAR  k: integer;
       fin,trouve: boolean;
  BEGIN
    { la variable REP contient le dernier résultat }
    New(vrbl[0]); New(formule[0]);
    vrbl[0]^:='REP'; formule[0]^:='';
    {take the expressions and  calculate result
     so long as a blank line is not found.}
    Clrscr;
    Writeln;
    Writeln('    *** CALCULATE ***');
    fin:=FALSE;
    REPEAT
      ii:=1;
      ok:=TRUE;
      writeln;
      write('  ');
      {read from keyboard the calculation to make}
      readln(prog);
      {normal exit with a blank line}
      IF length(prog)=0 THEN fin:=TRUE;
      IF NOT fin THEN
      BEGIN
        {case several calculations in a row}
        IF prog[1] IN DELIM THEN
          prog:=bak+prog
        ELSE
          bak:=prog;
        exist:=get_formula;           {treat assignation}
        trouve:=FALSE;
        {detect instruction with a variable name alone}
        FOR k:=0 TO num DO
	  IF vrbl[k]^=prog THEN
	  BEGIN
            trouve:=TRUE;
	    indice:=k
	  END;
        IF trouve THEN
	  IF formule[indice]=NIL THEN  {no formula defined}
          BEGIN
	    serror(3);
	    exist:=true;  {inhibit the calculation}
	  END
	  ELSE
	  BEGIN
	    FOR k:=0 TO num DO flg[k]:=FALSE; 
	    prog:=formule[indice]^;
	    writeln('  formula: ',prog)
          END;
        IF (exist=false) THEN  {no assignation}
        BEGIN
	  calcule_exp(reponse);
	  IF ok THEN
          BEGIN
	    writeln(reponse:16:8);
	    Str(reponse,formule[0]^);
	    bak:=formule[0]^
	  END
	  ELSE bak:=''
	END
      END
    UNTIL length(prog)=0;
  END;

  {List of current variables}
  PROCEDURE Affichage;
  VAR i:integer;
  BEGIN
    clrscr;
    writeln;
    writeln('              ***  LIST OF CURRENT VARIABLES  ***');
    writeln;
    writeln('            Name                       Formula or value');
    writeln('            ----                       ----------------');
    FOR i:=1 TO num DO
    BEGIN
      write(vrbl[i]^:15,formule[i]^:40);
      writeln
    END;
    readln
  END;

  PROCEDURE lecture_session;
  VAR i,j,n: integer;
      VAR1,form: S40;
      ligne    : STRING[65];
  BEGIN
    clrscr;
    writeln;
    writeln('     *** READ A SESSION ***');
    writeln;
    IF length(nom_IN)=0 THEN nom_IN:='c:\tpw\progr\calcul.asc';
    Assign(f_IN,nom_IN);
    Reset(f_IN);
    readln(f_IN,n);
    FOR i:=1 TO n+1 DO
    BEGIN
      readln(f_IN,ligne);
      j:=1; WHILE ligne[j]=' ' DO Inc(j);
      VAR1:=Copy(ligne,j,15-j+1);
      j:=16; WHILE ligne[j]=' ' DO Inc(j);
      form:=Copy(ligne,j,55-j+1);
      IF length(VAR1)>0 THEN
      BEGIN
	New(vrbl[i]);
	vrbl[i]^:=VAR1;
	New(formule[i]);
	formule[i]^:=form;
        num:=i
      END
    END;
    close(f_IN);
    writeln('  ',num,' items read.');
    readln
  END;

  PROCEDURE sauvegarde_session;
  VAR i: integer;
  BEGIN
    clrscr;
    writeln;
    writeln('     *** SAVE A SESSION ***');
    writeln;
    IF length(nom_out)=0 THEN nom_out:='c:\tpw\progr\calcul.asc';
    Assign(f_out,nom_out);
    Rewrite(f_out);
    writeln(f_out,num);
    FOR i:=1 TO num DO
      writeln(f_out,vrbl[i]^:15,formule[i]^:40);
    close(f_out);
    writeln('  ',num,' items written in ',nom_out);
    readln
  END;

  PROCEDURE menu;
  BEGIN
    clrscr;
    writeln;
    writeln('    ********************************************');
    writeln('    *              PARSER                      *');
    writeln('    *                                          *');
    writeln('    *  List variables and formulas       :  l  *');
    writeln('    *  Make calculations                 :  c  *');
    writeln('    *  Read a session from disk          :  r  *');
    writeln('    *  Save a session to disk            :  s  *');
    writeln('    *  Suppress a variable               :  x  *');
    writeln('    *  --------------------------------------  *');
    writeln('    *  Quit parser                       :  q  *');
    writeln('    *                                          *');
    writeln('    ********************************************');
    write('        Your choice ( l, c, r, s, x or q ) : ');
    choix:=readkey
  END;      

  PROCEDURE question;
  VAR c: char;
  BEGIN
    writeln;
    write('       Do you want to save your session (y/n): ');
    c:=readkey;
    IF c='y' THEN sauvegarde_session
  END;

  PROCEDURE suprim_variable;
  VAR nom: S40;
      i, ind: integer;
  BEGIN
    clrscr;
    writeln;
    writeln('       *** SUPPRESS A VARIABLE ***');
    writeln;
    write('  Name of the variable to suppress: ');
    readln(nom);
    FOR i:=1 TO num DO
      IF (vrbl[i]<>NIL) AND (vrbl[i]^=nom) THEN ind:=i;
    FOR i:=ind TO num-1 DO
    BEGIN
      vrbl[i]^:=vrbl[i+1]^;
      formule[i]^:=formule[i+1]^
    END;
    writeln('  Suppress done.');
    Dispose(vrbl[num]);
    Dispose(formule[num]);
    Dec(num)
  END;

  BEGIN   { programme principal }
    WinCrtInit(' PARSER');   {open main windows with caption}
    nom_IN:=ParamStr(1);
    nom_out:=ParamStr(2);
    REPEAT
      menu;
      CASE choix OF
      'l' : affichage;
      'c' : exec_calcul;
      'r' : lecture_session;
      's' : sauvegarde_session;
      'x' : suprim_variable;
      'q' : BEGIN
	      question;
	      FOR k:=0 TO num DO
	      BEGIN
		IF vrbl[k]<>NIL THEN Dispose(vrbl[k]);
		IF formule[k]<>NIL THEN Dispose(formule[k])
	      END;
	      DoneWinCrt
	    END
      END
    UNTIL choix='q'
  END.    	      	    

{end of file parser.pas}