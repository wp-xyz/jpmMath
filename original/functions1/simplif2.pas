{****************************************************************************
*                                                                           *
*               PROCEDURE TO SIMPLIFY A USER-DEFINED FUNCTION               *
*                         after a formal derivation                         *
*                                                                           *
* ------------------------------------------------------------------------- *
*               @Alain Reverchon et Marc Ducamp, 1987, France               *
*                                                                           *
*      OTHER USED FILE:                                                     *
*                                                                           *
*       UFUNCT.PAS                                                          *
*                                                                           *
****************************************************************************}
UNIT SIMPLIF2;      {TPW English Version By J-P Moreau, Paris}

INTERFACE

  USES UFunct;

  PROCEDURE SimplifyFunction (VAR f:FONC);


IMPLEMENTATION

PROCEDURE SimplifyFunction (VAR f:FONC);

VAR  i, j : integer;
     p1   : PtrNoeud;

BEGIN
  CASE f.typ OF
    OPUNAIRE  :
     BEGIN
       SimplifyFunction (f.left^);
       CASE f.oper1 OF
         FMOINS :
          BEGIN
      {----------------------------------------------- -(K) = (-K) }
            IF (f.left^.typ = CTE) THEN
               BEGIN
                 p1       := f.left;
                 f.typ    := CTE;
                 f.value := - p1^.value;
                 dispose (p1);
                 exit;
               END;
      {------------------------------------------------ -(-K) = K }
            IF (f.left^.typ = OPUNAIRE) AND (f.left^.oper1 = FMOINS) THEN
               BEGIN
                 p1       := f.left;
                 f        := p1^.left^;
                 dispose (p1^.left);
                 dispose (p1);
                 exit;
               END;
          END;
       END;  {CASE oper1}
     END; {OPUNAIRE}

    OPBINAIRE :
     BEGIN
       SimplifyFunction (f.left^);
       SimplifyFunction (f.right^);
       CASE f.oper2 OF
         FSUM :
          BEGIN
      {---------------------------------------------- CTE + CTE }
            IF (f.left^.typ = CTE) AND (f.right^.typ = CTE)
               THEN
               BEGIN
                 f.typ := CTE;
                 f.value := f.left^.value + f.right^.value;
                 dispose (f.left);
                 dispose (f.right);
                 exit;
               END;
      {---------------------------------------------- 0 + u }
            IF (f.left^.typ = CTE) AND (f.left^.value = 0)
               THEN
               BEGIN
                 dispose (f.left);
                 p1 := f.right;
                 f  := f.right^;
                 dispose (p1);
                 exit;
               END;
      {---------------------------------------------- u + 0 }
            IF (f.right^.typ = CTE) AND (f.right^.value = 0)
               THEN
               BEGIN
                 dispose (f.right);
                 p1 := f.left;
                 f  := f.left^;
                 dispose (p1);
                 exit;
               END;
      {---------------------------------------------- CTE + '+' }
            IF (f.left^.typ = CTE) AND (f.right^.typ = OPBINAIRE)
               AND (f.right^.oper2 = FSUM)
               THEN
               IF (f.right^.left^.typ = CTE) THEN
               BEGIN
                 f.left^.value := f.left^.value + f.right^.left^.value;
                 p1       := f.right;
                 f.right := p1^.right;
                 dispose (p1^.left);
                 dispose (p1);
                 exit;
               END
               ELSE
               IF (f.right^.right^.typ = CTE) THEN
               BEGIN
                 f.left^.value := f.left^.value + f.right^.right^.value;
                 p1       := f.right;
                 f.right := p1^.left;
                 dispose (p1^.right);
                 dispose (p1);
                 exit;
               END;
      {---------------------------------------------- '+' + CTE }
            IF (f.right^.typ = CTE) AND (f.left^.typ = OPBINAIRE)
               AND (f.left^.oper2 = FSUM)
               THEN
               IF (f.left^.left^.typ = CTE) THEN
               BEGIN
                 f.right^.value := f.right^.value + f.left^.left^.value;
                 p1       := f.left;
                 f.left := p1^.right;
                 dispose (p1^.left);
                 dispose (p1);
                 exit;
               END
               ELSE
               IF (f.left^.right^.typ = CTE) THEN
               BEGIN
                 f.right^.value := f.right^.value + f.left^.right^.value;
                 p1       := f.left;
                 f.left := p1^.left;
                 dispose (p1^.right);
                 dispose (p1);
                 exit;
               END;
          END;

         FSUB : BEGIN
      {---------------------------------------------- CTE - CTE }
                  IF (f.left^.typ = CTE) AND (f.right^.typ = CTE)
                     THEN
                     BEGIN
                       f.typ := CTE;
                       f.value := f.left^.value - f.right^.value;
                       dispose (f.left);
                       dispose (f.right);
                       exit;
                     END;
      {---------------------------------------------- u - 0 }
                  IF (f.left^.typ = CTE) AND (f.left^.value = 0)
                     THEN
                     BEGIN
                       dispose (f.left);
                       f.left := f.right;
                       f.typ := OPUNAIRE;
                       f.oper1 := FMOINS;
                       exit;
                     END;
      {---------------------------------------------- 0 - u }
                  IF (f.right^.typ = CTE) AND (f.right^.value = 0)
                     THEN
                     BEGIN
                       dispose (f.right);
                       p1 := f.left;
                       f  := f.left^;
                       dispose (p1);
                       exit;
                     END;
          END;

         FMUL : BEGIN
      {---------------------------------------------- CTE * CTE }
                  IF (f.left^.typ = CTE) AND (f.right^.typ = CTE)
                     THEN
                     BEGIN
                       f.typ := CTE;
                       f.value := f.left^.value * f.right^.value;
                       dispose (f.left);
                       dispose (f.right);
                       exit;
                     END;
      {---------------------------------------------- u * 0 }
                  IF (f.left^.typ = CTE) AND (f.left^.value = 0)
                     THEN
                     BEGIN
                       f.typ := CTE;
                       f.value := 0;
                       dispose (f.left);
                       dispose (f.right);
                       exit;
                     END;
      {---------------------------------------------- 0 * u }
                  IF (f.right^.typ = CTE) AND (f.right^.value = 0)
                     THEN
                     BEGIN
                       f.typ := CTE;
                       f.value := 0;
                       dispose (f.left);
                       dispose (f.right);
                       exit;
                     END;
      {---------------------------------------------- 1 * u }
                  IF (f.left^.typ = CTE) AND (f.left^.value = 1)
                     THEN
                     BEGIN
                       dispose (f.left);
                       p1 := f.right;
                       f := f.right^;
                       dispose (p1);
                       exit;
                     END;
      {---------------------------------------------- u * 1 }
                  IF (f.right^.typ = CTE) AND (f.right^.value = 1)
                     THEN
                     BEGIN
                       dispose (f.right);
                       p1 := f.left;
                       f  := f.left^;
                       dispose (p1);
                       exit;
                     END;
      {---------------------------------------------- -1 * u }
                  IF (f.left^.typ = CTE) AND (f.left^.value = -1)
                     THEN
                     BEGIN
                       f.typ    := OPUNAIRE;
                       f.oper1  := FMOINS;
                       p1       := f.left;
                       f.left := f.right;
                       dispose (p1);
                       exit;
                     END;
      {---------------------------------------------- u * -1 }
                  IF (f.right^.typ = CTE) AND (f.right^.value = -1)
                     THEN
                     BEGIN
                       f.typ    := OPUNAIRE;
                       f.oper1  := FMOINS;
                       dispose (f.right);
                       exit;
                     END;
      {---------------------------------------------- CTE + '*' }
            IF (f.left^.typ = CTE) AND (f.right^.typ = OPBINAIRE)
               AND (f.right^.oper2 = FMUL)
               THEN
               IF (f.right^.left^.typ = CTE) THEN
               BEGIN
                 f.left^.value := f.left^.value * f.right^.left^.value;
                 p1       := f.right;
                 f.right := p1^.right;
                 dispose (p1^.left);
                 dispose (p1);
                 exit;
               END
               ELSE
               IF (f.right^.right^.typ = CTE) THEN
               BEGIN
                 f.left^.value := f.left^.value * f.right^.right^.value;
                 p1       := f.right;
                 f.right := p1^.left;
                 dispose (p1^.right);
                 dispose (p1);
                 exit;
               END;
      {---------------------------------------------- '*' + CTE }
            IF (f.right^.typ = CTE) AND (f.left^.typ = OPBINAIRE)
               AND (f.left^.oper2 = FMUL)
               THEN
               IF (f.left^.left^.typ = CTE) THEN
               BEGIN
                 f.right^.value := f.right^.value * f.left^.left^.value;
                 p1       := f.left;
                 f.left := p1^.right;
                 dispose (p1^.left);
                 dispose (p1);
                 exit;
               END
               ELSE
               IF (f.left^.right^.typ = CTE) THEN
               BEGIN
                 f.right^.value := f.right^.value * f.left^.right^.value;
                 p1       := f.left;
                 f.left := p1^.left;
                 dispose (p1^.right);
                 dispose (p1);
                 exit;
               END;
          END;

         FDIV : BEGIN
      {---------------------------------------------- CTE / CTE }
                  IF (f.left^.typ = CTE) AND (f.right^.typ = CTE)
                     THEN
                     BEGIN
                       f.typ := CTE;
                       f.value := f.left^.value /f.right^.value;
                       dispose (f.left);
                       dispose (f.right);
                       exit;
                     END;
      {---------------------------------------------- 0 / u }
                  IF (f.right^.typ = CTE) AND (f.right^.value = 0)
                     THEN
                     BEGIN
                       f.typ := CTE;
                       f.value := 0;
                       dispose (f.left);
                       dispose (f.right);
                       exit;
                     END;
      {---------------------------------------------- u / 1 }
                  IF (f.right^.typ = CTE) AND (f.right^.value = 1)
                     THEN
                     BEGIN
                       dispose (f.right);
                       p1 := f.left;
                       f  := f.left^;
                       dispose (p1);
                       exit;
                     END;
      {---------------------------------------------- u / -1 }
                  IF (f.right^.typ = CTE) AND (f.right^.value = -1)
                     THEN
                     BEGIN
                       f.typ    := OPUNAIRE;
                       f.oper1  := FMOINS;
                       dispose (f.right);
                       exit;
                     END;
                END;

         FPUI : BEGIN
      {---------------------------------------------- CTE ^ CTE }
                  IF (f.left^.typ = CTE) AND (f.right^.typ = CTE)
                     THEN
                     BEGIN
                       f.typ := CTE;
                       f.value := exp(ln(f.left^.value)*f.right^.value);
                       dispose (f.left);
                       dispose (f.right);
                       exit;
                     END;
      {---------------------------------------------- u ^ 0 }
                  IF (f.left^.typ = CTE) AND (f.left^.value = 0)
                     THEN
                     BEGIN
                       f.typ := CTE;
                       f.value := 1;
                       dispose (f.left);
                       dispose (f.right);
                       exit;
                     END;
      {---------------------------------------------- 0 ^ u }
                  IF (f.right^.typ = CTE) AND (f.right^.value = 0)
                     THEN
                     BEGIN
                       f.typ := CTE;
                       f.value := 0;
                       dispose (f.left);
                       dispose (f.right);
                       exit;
                     END;
      {---------------------------------------------- 1 ^ u }
                  IF (f.left^.typ = CTE) AND (f.left^.value = 1)
                     THEN
                     BEGIN
                       f.typ := CTE;
                       f.value := 1;
                       dispose (f.left);
                       dispose (f.right);
                       exit;
                     END;
      {---------------------------------------------- u ^ 1 }
                  IF (f.right^.typ = CTE) AND (f.right^.value = 1)
                     THEN
                     BEGIN
                       dispose (f.right);
                       p1 := f.left;
                       f  := f.left^;
                       dispose (p1);
                       exit;
                     END;
      {---------------------------------------------- u ^ -1 }
                  IF (f.right^.typ = CTE) AND (f.right^.value = -1)
                     THEN
                     BEGIN
                       f.typ     := OPBINAIRE;
                       f.oper2   := FDIV;
                       f.right^ := f.left^;
                       f.left^.typ    := CTE;
                       f.left^.value := 1;
                       exit;
                     END;
                END;
          END; {CASE oper2}
     END;
  END; {CASE}
END;

END.