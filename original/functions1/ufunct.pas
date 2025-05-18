{========================================================================*
*                                                                        *
*       PROCEDURES TO COMPILE AND EVALUATE USER DEFINED FUNCTIONS        *
*                  DEFINED BY ITS EQUATION F(x) OR F(t)                  *
*                                                                        *
*(This unit must be included in any program using mathematical functions)*
*                                                                        *
*                                 By Alain Reverchon, Paris [BIBLI 03].  *
* ---------------------------------------------------------------------- *
*  English reduced and autonomous version of FONCTION.INC for Turbo      *
*  Pascal for Windows.  Visible Functions:   Evaluate, CompileFunction,  *
*  EnterFunction, RealValue and DestroyFunction.                         *
*                                                                        *
*                                                  By J-P Moreau Paris.  *
*                                                   (www.jpmoreau.fr)    *
*========================================================================}
Unit UFunct;

Interface
Uses Type_def;  { for type REAL_AR : real or double

*========================================================================*
*                   SPECIFIC TYPES AND CONSTANTS                         *
*========================================================================}
TYPE
  TypeNoeud      = (CTE, VARIABLE, OPUNAIRE, OPBINAIRE);
  UneOperUnaire  =
    (FMOINS, FPLUS, FSIN, FCOS, FTAN, FARCSIN, FARCCOS,
     FARCTAN, FSH,  FCH, FTH, FARGSH, FARGCH, FARGTH,
     FABS, FINT, FFRAC, FFACT, FEXP, FLOG, FSQR);
  UneOperBinaire = (FSUM, FSUB, FMUL, FDIV, FPUI);

  PtrNoeud = ^noeud;
  noeud =
    RECORD
      left, right  : PtrNoeud;
      CASE typ : TypeNoeud OF
	CTE        : (value : Real_ar);
        VARIABLE   : (rang   : integer);
        OPUNAIRE   : (oper1  : UneOperUnaire);
        OPBINAIRE  : (oper2  : UneOperBinaire);
    END;

  FONC    = noeud;
  TWOFONC = ARRAY [1..2] OF FONC;
  CHAINE  = string[80];
  MAXSTRING = string[255];

CONST
  OperationsUnaires : ARRAY [UneOperUnaire] OF STRING [6] =
     ('-', '+',
      'SIN','COS','TAN','ARCSIN','ARCCOS','ARCTAN',
      'SH', 'CH', 'TH', 'ARGSH', 'ARGCH', 'ARGTH',
      'ABS','INT','FRAC', 'FACT', 'EXP','LN','SQR');

  lastfonc : CHAINE = '';

VAR f1, f2, f3, f4 : FONC;
    err : boolean;
    param : array[0..99] of Real_ar;

  {Visible functions for calling program}
  Function  Evaluate (f : FONC; x : Real_ar) : Real_ar;
  Function  CompileFunction (ch:MAXSTRING; VAR f:FONC) : boolean;
  Function  EnterFunction (ch:CHAINE; VAR f:FONC) : boolean;
  Function  RealValue (ch:CHAINE): Real_ar;
  Procedure DestroyFunction (VAR f:FONC);


Implementation

{**********************************************************
*              EVALUATE A FUNCTION AT POINT X             *
**********************************************************}
Function Evaluate (f : FONC; x : Real_ar) : Real_ar;
VAR  i, j : integer;
     r, s : Real_ar;
BEGIN
  err    := FALSE;
  evaluate := 0;
  CASE f.typ OF
    CTE       : evaluate := f.value;
    VARIABLE  : IF (f.rang = -1)
                   THEN evaluate := x
                   ELSE evaluate := param [f.rang];
    OPUNAIRE  : BEGIN
                   r := evaluate (f.left^, x);
                   IF err THEN exit;
                   CASE f.oper1 OF
                     FMOINS : evaluate := - r;
                     FPLUS  : evaluate := r;
                     FABS   : evaluate := abs  (r);
                     FINT   : evaluate := int  (r);
                     FFRAC  : evaluate := frac (r);
                     FSQR   : IF (r >= 0)
                                 THEN evaluate := sqrt (r)
                                 ELSE err    := TRUE;
                     FEXP   : evaluate := exp (r);
                     FLOG   : IF (r > 0)
                                 THEN evaluate := ln (r)
                                 ELSE err    := TRUE;
                     FSIN : evaluate := sin (r);
                     FCOS : evaluate := cos (r);
                     FTAN : BEGIN
                              s := cos (r);
                              IF (s <> 0)
                                 THEN evaluate := sin (r) / s
                                 ELSE err    := TRUE;
                            END;
                     FSH  : evaluate := (exp (r) - exp (-r)) / 2;
                     FCH  : evaluate := (exp (r) + exp (-r)) / 2;
                     FTH  : evaluate := (exp (2 * r) - 1) / (exp (2 * r) + 1);
                  FARCSIN : IF (r = -1)
                             THEN evaluate := -PI/2
                             ELSE
                              IF (r = 1)
                               THEN evaluate := PI/2
                               ELSE
                                IF (r > -1) AND (r < 1)
                                 THEN evaluate := ArcTan (r / Sqrt (1 - Sqr(r)))
                                 ELSE err    := TRUE;
                  FARCCOS : IF (r = -1)
                             THEN evaluate := PI
                             ELSE
                              IF (r = 1)
                               THEN evaluate := 0
                               ELSE
                                IF (r > -1) AND (r < 1)
                                 THEN evaluate := PI/2 - ArcTan (r / Sqrt (1 - Sqr(r)))
                                 ELSE err    := TRUE;
                  FARCTAN : evaluate := Arctan (r);
                   FARGSH : evaluate := ln (r + sqr (1 + r * r));
                   FARGCH : IF (r >= 1)
                                 THEN evaluate := ln (r + sqr (-1 + r * r))
                                 ELSE err    := TRUE;
                   FARGTH : IF (r > -1) AND (r < 1)
                                 THEN evaluate := ln ((1 + r) / (1 - r)) / 2
                                 ELSE err    := TRUE;
                    FFACT : IF (abs (r) < 32)
                              THEN
                                BEGIN
                                s := 1;
                                FOR i:=1 TO ROUND (r) DO s := s * i;
                                evaluate := s;
                                END
                              ELSE err    := TRUE;
                 END; {CASE}
                 END;
    OPBINAIRE : BEGIN
                   r := evaluate (f.left^, x);
                   IF err THEN exit;
                   s := evaluate (f.right^, x);
                   IF err THEN exit;
                   CASE f.oper2 OF
                     FSUM : evaluate := r + s;
                     FSUB : evaluate := r - s;
                     FMUL : evaluate := r * s;
                     FDIV : IF (s <> 0)
                               THEN evaluate := r / s
                               ELSE err    := TRUE;
                     FPUI : IF (f.right^.typ IN [CTE, VARIABLE])
                              AND (int (s) = s) THEN
                               BEGIN
                                 j := ROUND (s);
                                 IF (j < 10)
                                  THEN
                                  BEGIN
                                    s := 1;
                                    FOR i:=1 TO j DO s := s * r;
                                    evaluate := s;
                                  END
                                  ELSE
                                  IF (r = 0)
                                   THEN evaluate := 0
                                   ELSE
                                   IF odd (j)
                                     THEN
                                      IF (r > 0)
                                        THEN evaluate := exp (ln (r) * s)
                                        ELSE evaluate := - exp (ln (-r) * s)
                                     ELSE  evaluate := exp (ln (abs (r)) * s);
                               END
                               ELSE
                               IF (r = 0) THEN evaluate := 0
                               ELSE
                                 IF (r > 0) THEN evaluate := exp (s * ln (r))
                                            ELSE err    := TRUE;
                 END; {CASE oper2}
          END; {cas des operations binaires}
    END; {CASE}
END;


{************************************************************************
*                         COPY A FUNCTION                               *
************************************************************************}
Procedure CopyFunction (f:FONC; VAR g:FONC);
BEGIN
  g := f;
  CASE f.typ OF
    OPUNAIRE  : BEGIN
                   new (g.left);
                   CopyFunction (f.left^, g.left^);
                 END;
    OPBINAIRE : BEGIN
                   new (g.left);
                   CopyFunction (f.left^, g.left^);
                   new (g.right);
                   CopyFunction (f.right^, g.right^);
                 END;
  END;
END;

{************************************************************************
*                        ERASE A FUNCTION                               *
************************************************************************}
Procedure DestroyFunction (VAR f:FONC);
BEGIN
  WITH f DO
  CASE typ OF
    OPBINAIRE : BEGIN
                   DestroyFunction (right^);
                   dispose (right);
                   DestroyFunction (left^);
                   dispose (left);
                 END;
    OPUNAIRE  : BEGIN
                   DestroyFunction (left^);
                   dispose (left);
                 END;
  END;
END;

{************************************************************************
*                    CREATE THE FUNCTION TREE                           *
************************************************************************}
Function CreateTree (ch:MAXSTRING; VAR n:noeud) : boolean;

VAR  i, j, l     : integer;
     r           : Real_ar;
     k           : UneOperUnaire;
     p1, p2      : Ptrnoeud;

Function SkipParenthesis (VAR i: integer; sens:integer) : boolean;
VAR  np : integer;
BEGIN
  np  := 0;
  REPEAT
    CASE ch [i] OF
      '(' : np := np + 1;
      ')' : np := np - 1;
    END;
    i := i + sens;
  UNTIL (np = 0) OR (i > l) OR (i < 1);
  i := i - sens;
  SkipParenthesis := (np = 0);
END;

Procedure Cut (c:CHAINE; x:UneOperUnaire);
BEGIN
  n.typ   := OPUNAIRE;
  n.oper1 := x;
  new (p1);
  n.left := p1;
  CreateTree := CreateTree (c, n.left^);
END;

Procedure CutBi (x : UneOperBinaire);
BEGIN
  n.typ   := OPBINAIRE;
  n.oper2 := x;
  new (p1); new (p2);
  n.left := p1; n.right := p2;
  IF NOT CreateTree (copy (ch, 1,   i-1), n.left^) THEN exit;
  IF NOT CreateTree (copy (ch, i+1, l-i), n.right^) THEN exit;
  CreateTree := TRUE;
END;

Function Seek (c1,c2:char; x1,x2:UneOperBinaire) : boolean;
BEGIN
  Seek := FALSE;
  i := l;
  REPEAT
     IF (ch [i] = ')') THEN
       IF NOT SkipParenthesis (i, -1) THEN
         BEGIN
         Seek := TRUE;
         exit;
         END;
     IF (ch [i] = c1) THEN
       IF (i > 1) AND (i < l) AND NOT (ch [i-1] IN ['(', 'E', '*', '/', '^']) THEN
          BEGIN
          CutBi (x1);
          Seek := TRUE;
          exit;
          END;
     IF (ch [i] = c2) THEN
       IF (i > 1) AND (i < l) AND NOT (ch [i-1] IN ['(', 'E', '*', '/', '^']) THEN
          BEGIN
          CutBi (x2);
          Seek := TRUE;
          exit;
          END;
     i := i - 1;
   UNTIL (i <= 1);
END;

BEGIN  {CreateTree}
   {------------------------ 1: Initializations }
   CreateTree := FALSE;
   n.typ     := CTE;
   l := length (ch);
   IF (l = 0) THEN exit;

   {------------------------ 2: Seek binary operations }
   IF Seek ('-', '+', FSUB, FSUM) THEN exit;
   IF Seek ('/', '*', FDIV, FMUL) THEN exit;
  {------------------------ 4:  Seek unary - and + }
  IF (ch [1] = '-') THEN
    BEGIN
    Cut (copy (ch, 2, l-1), FMOINS);
    exit;
    END;
  IF (ch [1] = '+') THEN
    BEGIN
    Cut (copy (ch, 2, l-1), FPLUS);
    exit;
    END;

   IF Seek ('^', '^', FPUI, FPUI) THEN exit;

  {------------------------ 3: Seek parentheses }
   i := 1;
   IF (ch [1] = '(') THEN
     BEGIN
       IF SkipParenthesis (i, 1) THEN
          IF (i = l) THEN CreateTree := CreateTree (copy (ch, 2, l-2), n);
       exit;
     END;

  {------------------------ 5:  Seek sign ! }
  IF (ch [l] = '!') THEN
    BEGIN
    Cut (copy (ch, 1, l-1), FFACT);
    exit;
    END;

  {------------------------ 6:  Seek unary functions }
   FOR k := FSIN TO FSQR DO
   BEGIN
     j := length (OperationsUnaires [k]);
     IF (OperationsUnaires [k] = copy (ch, 1, j)) THEN
       BEGIN
       IF (ch [j + 1] = '(') THEN Cut(copy (ch, j+1, l-j), k);
       exit;
       END;
   END;

  {------------------------ 7:  Seek sub-functions }
   IF (ch = 'U') OR (ch = 'U''') OR (ch = 'V') OR (ch = 'V''') THEN
     BEGIN
     IF (ch = 'U')   THEN CopyFunction (f1, n);
     IF (ch = 'U''') THEN CopyFunction (f2, n);
     IF (ch = 'V')   THEN CopyFunction (f3, n);
     IF (ch = 'V''') THEN CopyFunction (f4, n);
     CreateTree := TRUE;
     exit;
     END;

  {------------------------ 8:  Seek unknownss }
   IF (ch = 'X') OR (ch = 'T') OR (ch = 'I') OR (ch = 'N') THEN
    BEGIN
     n.typ      := VARIABLE;
     n.rang     := -1;
     CreateTree := TRUE;
     exit;
    END;
   IF (ch [1] = 'Y') THEN
    BEGIN
      IF (l = 1)
       THEN n.rang := 0
       ELSE
        BEGIN
         IF (ch [2] = '''')
          THEN
           BEGIN
             n.rang := l - 1;
             FOR j := 2 TO l DO
              IF (ch [j] <> '''') THEN exit;
           END
          ELSE
           BEGIN
            val (copy (ch, 2, l-1), i, j);
            IF (j <> 0) THEN exit;
            n.rang := i;
           END;
        END;
      n.typ      := VARIABLE;
      CreateTree := TRUE;
      exit;
    END;

  {------------------------ 9:   Seek predefined constants }
   IF (ch = 'PI') THEN
   BEGIN
     n.typ      := CTE;
     n.value   := PI;
     CreateTree := TRUE;
     exit;
   END;

  {------------------------ 10:  Seek numerical constants }
   val (ch, r, i);
   IF (i = 0) THEN
   BEGIN
     n.typ      := CTE;
     n.value   := r;
     CreateTree := TRUE;
   END;
END; {Function CreateTree}


{************************************************************************
*                 COMPILE A USER DEFINED FONCTION                       *
************************************************************************}
Function CompileFunction (ch:MAXSTRING; VAR f:FONC) : boolean;
VAR  i, j : integer;
BEGIN
   CompileFunction := FALSE;
   j := 0;
   FOR i:=1 TO length (ch) DO
    IF (ch [i] <> ' ') THEN
     BEGIN
       j:= j + 1;
       ch [j] := upcase (ch [i]);
     END;
   ch [0] := chr (j);
   IF CreateTree (ch, f)
      THEN CompileFunction := TRUE
      ELSE CompileFunction := FALSE
END;

{************************************************************
* Transform the string ch into a compiled function pointed  *
* to by f of type FONC. Return TRUE if the compilation      *
* succeeded or FALSE if the compilation failed (syntax      *
* error or unknown operator). Accepted variables are x, t.  *
************************************************************}   
FUNCTION EnterFunction (ch:CHAINE; VAR f:FONC) : boolean;
VAR  i  : integer;
BEGIN
  EnterFunction := TRUE;
  IF (ch = '') THEN exit
               ELSE lastfonc := ch;
  IF CompileFunction (lastfonc, f) THEN exit;
  EnterFunction := FALSE;
END;

{convert a string (ex: -PI/2 ) into a real number}
FUNCTION  RealValue (ch:CHAINE): Real_ar;
VAR  local_f:FONC; 
BEGIN
  IF CompileFunction(ch,local_f) THEN
  BEGIN
    RealValue:=Evaluate(local_f,0);
    IF NOT err THEN exit;
  END;
  RealValue:=0
END;


END.   { UFunct.pas }