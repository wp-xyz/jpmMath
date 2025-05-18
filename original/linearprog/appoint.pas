{********************************************************
*                 APPOINTMENT METHOD                    *
*                 ------------------                    *
*                                                       *
* LIST OF MAIN VARIABLES:                               *
*                                                       *
*  NP:         NUMBER OF JOBS                           *
*  C(NP,NP):   APPOINTMENT COST MATRIX                  *
*  MP(NP,NP):  JOB/APPLICANT MATRIX                     *
*  IF1:        FLAG=1 IF OPTIMAL APPOINTMENT IS FOUND   *
*  IZ:         FLAG TO MARK A ZERO                      *
*  XMIN:       MINIMUM VALUE                            *
*  XCASE:      APPOINTMENT OF A ZERO IN RELATION WITH   *
*              THE ZEROES' MINIMUM                      *
*  NZ:         NUMBER OF ZEROES                         *
*  ZI,ZJ:      COORDINATES OF CURRENT CASE              *
*  M,N:        COORDINATES OF A MARKED ZERO             *
*  A,B:        AUXILIARY VARIABLES TO CHANGE            *
*              APPOINTMENTS.                            *
* ----------------------------------------------------- *
* PROBLEM DESCRIPTION                                   *
* An employer receives four applicants for four avail-  *
* able jobs. They must give a mark from 1 to 100 to     *
* each job (1=very poor preference, 100=very good).     *
* See table below:                                      *
*                    J O B S                            *
*                 J1  J2  J3  J4                        *
*              A1 90  80  20  40                        *
*  APPLICANTS  A2 90  70  30  80                        *
*              A3 40  70  20  80                        *
*              A4 50  40  20  60                        *
* The employer's purpose is to appoint the four jobs so *
* as to maximize the applicants' satisfaction.          *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
* (appoint four jobs with the following regrets matrix: *
*                    J O B S                            *
*                 J1  J2  J3  J4                        *
*              A1 10  20  80  60                        *
*  APPLICANTS  A2 10  30  70  20                        *
*              A3 60  30  80  20                        *
*              A4 50  60  80  40   )                    *
*                                                       *
* LINEAR PROGRAMMING                                    *
*                                                       *
* APPOINTMENT METHOD                                    *
*                                                       *
* NUMBER OF JOBS ? 4                                    *
*                                                       *
* INPUT APPOINTMENT COSTS/REGRETS OF APPLICANTS:        *
*                                                       *
* APPLICANT #1:                                         *
*    JOB #1 ? 10                                        *
*    JOB #2 ? 20                                        *
*    JOB #3 ? 80                                        *
*    JOB #4 ? 60                                        *
*                                                       *
* APPLICANT #2:                                         *
*    JOB #1 ? 10                                        *
*    JOB #2 ? 30                                        *
*    JOB #3 ? 70                                        *
*    JOB #4 ? 20                                        *
*                                                       *
* APPLICANT #3:                                         *
*    JOB #1 ? 60                                        *
*    JOB #2 ? 30                                        *
*    JOB #3 ? 80                                        *
*    JOB #4 ? 20                                        *
*                                                       *
* APPLICANT #4:                                         *
*    JOB #1 ? 50                                        *
*    JOB #2 ? 60                                        *
*    JOB #3 ? 80                                        *
*    JOB #4 ? 40                                        *
*                                                       *
* APPOINTMENTS:                                         *
*                                                       *
*      APPLICANT #1 => JOB #2                           *
*      APPLICANT #2 => JOB #1                           *
*      APPLICANT #3 => JOB #4                           *
*      APPLICANT #4 => JOB #3                           *
*                                                       *
* ----------------------------------------------------- *
* REFERENCE:                                            *
* Modèles pratiques de décision Tome 2, By Jean-Pierre  *
* Blanger, PSI Editions, France, 1982.                  *
*                                                       *
*            Pascal Release 1.0 By J-P Moreau, Paris.   *
*                       (www.jpmoreau.fr)               *
********************************************************}
PROGRAM APPOINTMENT;

Uses WinCrt;

Const
      NMAX = 10;

Var
    C: Array[0..NMAX,0..NMAX] of REAL;
    MP: Array[0..NMAX,0..NMAX] of Integer;
    IF1,IZ,M,N,NP,NZ,ZI,ZJ: Integer;
    A,B,XCASE,XMIN: REAL;


Procedure Data;  {INPUT COST/APPOINT}
Var I,J: Integer;
Begin
  writeln;
  writeln(' LINEAR PROGRAMMING');
  writeln;
  writeln(' APPOINTMENT METHOD');
  writeln;
  write(' NUMBER OF JOBS ? '); readln(NP);
  writeln;
  writeln(' INPUT APPOINTMENT COSTS/REGRETS OF APPLICANTS:');
  FOR I := 1 TO NP DO
  begin
    writeln;
    writeln(' APPLICANT #', I,':');
    FOR J := 1 TO NP DO
    begin
      write('    JOB #',J,' ? ');
      readln(C[I, J])
    end
  end;
  writeln;
  writeln(' APPOINTMENTS:');
  writeln
End;

    Procedure Zeroes; Forward;
    Procedure Appoint; Forward;
    Procedure Mark; Forward;
    Procedure Sub_add; Forward;

Procedure Main;  {HUNGARIAN ALGORITHM}
Label 10;
Begin
10: Zeroes;
    Appoint;
    Mark;
    Sub_add;
    Appoint;
    IF IF1 <> NP THEN GOTO 10
End;

Procedure Zeroes;
Label 100,200;
Var I,J: Integer;
    R:REAL;
Begin
     IZ := 0; R := 0;
     FOR I := 1 TO NP DO
     begin
       XMIN := C[I, 1];
       FOR J := 1 TO NP DO
       begin
	 IF C[I, J] = 0 THEN IZ := 1;
	 IF C[I, J] < XMIN THEN XMIN := C[I, J]
       end;
       IF IZ = 1 THEN
       begin
         IZ := 0; GOTO 100
       end;
       FOR J := 1 TO NP DO
	 C[I, J] := C[I, J] - XMIN;
100: end;
     FOR J := 1 TO NP DO
     begin
       XMIN := C[1, J];
       FOR I := 1 TO NP DO
       begin
	 IF C[I, J] = 0 THEN IZ := 1;
	 IF C[I, J] < XMIN THEN XMIN := C[I, J]
       end;
       IF IZ = 1 THEN
       begin
	 IZ := 0; GOTO 200
       end;
       FOR I := 1 TO NP DO
	 C[I, J] := C[I, J] - XMIN;
200: end
End;

Procedure Appoint;
Var I,J,K: Integer;
Label 10;
Begin
     FOR I := 1 TO NP DO
     begin
       FOR J := 1 TO NP DO
       begin
	 MP[I, J] := 0; C[0, J] := 0
       end;
       C[I, 0] := 0
     end;
     FOR I := 1 TO NP DO
     begin
       XCASE := 999999;
       FOR J := 1 TO NP DO
       begin
	 IF (C[I, J] <> 0) OR (MP[I, J] <> 0) THEN GOTO 10;
	 NZ := 0;
	 FOR K := 1 TO NP DO
	   IF C[K, J] = 0 THEN Inc(NZ);
	 IF NZ < XCASE THEN
         begin
	   XCASE := 1.0*NZ; ZI := I; ZJ := J
	 end;
10:    end;
       MP[ZI, ZJ] := 1;
       FOR K := 1 TO NP DO
	 IF (C[K, ZJ] = 0) AND (MP[K, ZJ] = 0) THEN MP[K, ZJ] := -1;
       FOR K := 1 TO NP DO
	 IF (C[ZI, K] = 0) AND (MP[ZI, K] = 0) THEN MP[ZI, K] := -1
     end;
     IF1 := 0;
     FOR I := 1 TO NP DO
       FOR J := 1 TO NP DO
	 IF MP[I, J] = 1 THEN Inc(IF1)
End;

Procedure Mark;
Label 10;
Var I,J: Integer;
Begin
10:  FOR I := 1 TO NP DO
     begin
       N := 0;
       FOR J := 1 TO NP DO
	 IF MP[I, J] = 1 THEN N := 1;
       IF (N = 0) AND (C[I, 0] = 0) THEN
       begin
	 C[I, 0] := 1; M := 1
       end
     end;
     FOR J := 1 TO NP DO
       FOR I := 1 TO NP DO
	 IF (MP[I, J] = -1) AND (C[I, 0] = 1) AND (C[0, J] = 0) THEN
         begin
	   C[0, J] := 1; M := 1
	 end;
     FOR I := 1 TO NP DO
       FOR J := 1 TO NP DO
	 IF (MP[I, J] = 1) AND (C[0, J] = 1) AND (C[I, 0] = 0) THEN
         begin
	   C[I, 0] := 1; M := 1
	 end;
     IF M = 1 THEN
     begin
       M := 0; GOTO 10
     end
End;

Procedure Sub_add;
Var I,J: Integer;
Begin
  XMIN := 999999;
  FOR I := 1 TO NP DO
    FOR J := 1 TO NP DO
    begin
      A := C[I, 0]; B := C[0, J];
      IF (A = 1) AND (B = 0) AND (C[I, J] < XMIN) THEN XMIN := C[I, J]
    end;
  FOR I := 1 TO NP DO
    FOR J := 1 TO NP DO
    begin
      A := C[I, 0]; B := C[0, J];
      IF (A = 1) AND (B = 0) THEN C[I, J] := C[I, J] - XMIN;
      IF (A = 0) AND (B = 1) THEN C[I, J] := C[I, J] + 2 * XMIN
    end
End;

Procedure Results;
Label 10;
Var I,J: Integer;
Begin
    FOR I := 1 TO NP DO
      FOR J := 1 TO NP DO
      begin
        IF MP[I, J] <> 1 THEN GOTO 10;
	writeln('      APPLICANT #',I,' => JOB #', J);
10:   end;
    writeln
End;

{main program}
BEGIN
  ClrScr;
  Data;      {INPUT COST/APPOINT}
  Main;      {HUNGARIAN ALGORITHM}
  Results;   {PRINT RESULTS}
  ReadKey;
  DoneWinCrt
END.

{END OF FILE APPOINT.PAS}