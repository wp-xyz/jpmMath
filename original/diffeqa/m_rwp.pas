{***********************************************************************
*  Solve a boundary value problem for a first order DE system of the   *
*  form:                                                               *
*                                                                      *
*                   Y1' = F1(X,Y1,Y2,...,YN)                           *
*                   Y2' = F2(X,Y1,Y2,...,YN)                           *
*                             ...                                      *
*                   YN' = FN(X,Y1,Y2,...,YN)                           *
*                                                                      *
*  via the shooting method by determining an approximation for the     *
*  initial value Y(A).                                                 *
* -------------------------------------------------------------------- *
* SAMPLE RUNS;                                                         *  
*                                                                      *
*   TEST EXAMPLE (WITH METHOD #1):                                     *
*   Method: Runge-Kutta embedding formula of 4/5th order.              *
*   =============                                                      *
*   ANALYZED SET OF DIFFERENTIAL EQUATIONS:                            *
*   ---------------------------------------                            *
*            Y'(1) = Y(2)                                              *
*            Y'(2) = -Y(1)**3                                          *
*   WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)                          *
*     AT POINT  A =      0.000:      0.000                             *
*     AT POINT  B =      1.000:      0.000                             *
*                                                                      *
*   AT LOCATION  A =      0.000 THE FOLLOWING VALUES ARE PROVIDED:     *
*                                                                      * 
*     Y(1) =  0.00000000000000E+0000                                   *
*     Y(2) =  1.20000000000000E+0001                                   *
*                                                                      *
*   REQUIRED PARAMETER:                                                *
*   ------------------:                                                *
*   -PRECISION FOR THE IVP =        0.0000001000                       * 
*   -PRECISION FOR THE BVP =        0.0000000100                       *
*   -STEP WIDTH FOR THE IVP=        0.0100000000                       *
*   -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  20000   *
*   -MAX. NUMBER OF NEWTON-ITERATION STEPS =   1000                    *
*                                                                      *
*   SOLUTION OF THE PROBLEM:                                           * 
*   ------------------------                                           *
*   THE START VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED           *
*   IN POINT A =      0.000 AS FOLLOWS:                                *
*                                                                      *
*     Y(1) =         0.0000000000                                      *
*     Y(2) =         9.7229810241                                      *
*                                                                      *
*   THIS REQUIRES 4 NEWTON-ITERATIONS!                                 *
*                                                                      *
*   DETERMINATION COMPLETED AS PLANNED!                                *
*                                                                      *  
*                                                                      *
*   TEST EXAMPLE (WITH METHOD #2):                                     *
*   Method: Predictor-corrector method of order 4 by Adams-Bashforth-  *
*           Moulton.                                                   *
*   =============                                                      *
*   ANALYZED SET OF DIFFERENTIAL EQUATIONS:                            *
*   ---------------------------------------                            *
*            Y'(1) = Y(2)                                              *
*            Y'(2) = -Y(1)**3                                          *
*   WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)                          *
*     AT POINT  A =      0.000:      0.000                             *
*     AT POINT  B =      1.000:      0.000                             *
*                                                                      * 
*   AT LOCATION  A =      0.000 THE FOLLOWING VALUES ARE PROVIDED:     *
*                                                                      *
*     Y(1) =  0.00000000000000E+0000                                   *
*     Y(2) =  1.20000000000000E+0001                                   *
*                                                                      * 
*   REQUIRED PARAMETER:                                                *
*   ------------------:                                                *
*   -PRECISION FOR THE IVP =    0.0000001000                           *
*   -PRECISION FOR THE BVP =    0.0000000100                           * 
*   -STEP WIDTH FOR THE IVP=    0.0100000000                           *
*   -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  20000   *
*   -MAX. NUMBER OF NEWTON-ITERATION STEPS =   1000                    *
*                                                                      *
*   SOLUTION OF THE PROBLEM:                                           *
*   ------------------------                                           *
*   THE START VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED           *
*   IN POINT A =      0.000 AS FOLLOWS:                                *
*                                                                      *
*     Y(1) =         0.0000000000                                      *
*     Y(2) =         9.7229809524                                      *
*                                                                      *
*   THIS REQUIRES 4 NEWTON-ITERATIONS!                                 * 
*                                                                      *
*   DETERMINATION COMPLETED AS PLANNED!                                *
*                                                                      *  
*                                                                      *
*   TEST EXAMPLE (WITH METHOD #3):                                     *
*   Method: Extrapolation method of  Bulirsch-Stoer.                   *
*   =============                                                      *
*   ANALYZED SET OF DIFFERENTIAL EQUATIONS:                            *
*   ---------------------------------------                            *
*            Y'(1) = Y(2)                                              *
*            Y'(2) = -Y(1)**3                                          *
*   WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)                          *
*     AT POINT  A =      0.000:      0.000                             *
*     AT POINT  B =      1.000:      0.000                             * 
*                                                                      *
*   AT LOCATION  A =      0.000 THE FOLLOWING VALUES ARE PROVIDED:     *
*                                                                      *
*     Y(1) =  0.00000000000000E+0000                                   *
*     Y(2) =  1.20000000000000E+0001                                   *
*                                                                      *
*   REQUIRED PARAMETER:                                                *
*   ------------------:                                                *
*   -PRECISION FOR THE IVP =    0.0000001000                           *
*   -PRECISION FOR THE BVP =    0.0000000100                           * 
*   -STEP WIDTH FOR THE IVP=    0.0100000000                           *
*   -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  20000   *
*   -MAX. NUMBER OF NEWTON-ITERATION STEPS =   1000                    *
*                                                                      *
*   SOLUTION OF THE PROBLEM:                                           *
*   ------------------------                                           *
*   THE START VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED           *
*   IN POINT A =      0.000 AS FOLLOWS:                                *
*                                                                      *
*     Y(1) =     0.0000000000                                          *
*     Y(2) =     9.7229810195                                          *
*                                                                      *
*   THIS REQUIRES 4 NEWTON-ITERATIONS!                                 *
*                                                                      *
*   DETERMINATION COMPLETED AS PLANNED!                                *
*                                                                      *  
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C,  By Gisela Engeln-Muellges       *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
*                                                                      *
*                                  TPW Release By J-P Moreau, Paris.   *
*                                          (www.jpmoreau.fr)           *
************************************************************************
   Initialize data; test examples can be exchanged by adjusting
   the functions dgl() in unit Uawp, randbed() in unit Urwp and
   the variable N in main program.                                     }
Program Test_Rwp;

Uses WinCrt1, Fgauss, Urwp, T_dgls;

Var
     TMethod: Array[1..5] of String[80];
     yanf,yanf0: pVec;
     N1: Integer;

     a, b, h, wa, wb: Double;
     ifmax, itmax: integer;
     epsawp, epsrb: Double; 
     fehler, iter, method, j: integer;


{Init table of method names}
Procedure Init_TMethod;
Begin
  TMethod[1]:='Runge-Kutta embedding formula of 4/5th order.';
  TMethod[2]:='Predictor-corrector method of order 4 by Adams-Bashforth-Moulton.';
  TMethod[3]:='Extrapolation method of  Bulirsch-Stoer.';
End;

Procedure Copy_Vector(a,b:pVec;n:integer);
Var i:integer;
Begin
  For i:=0 to n-1 do a^[i]:=b^[i]
End;


{main program}
BEGIN
                  
  New(yanf0); New(yanf);
  Init_TMethod;

  N1:=2;                     {Number of equations}

  yanf0^[0] :=  0.0;
  yanf0^[1] := 12.0;

  bspnummer := 2;            {Example #2 - see t_dgls.pas}
  a := 0.0;                  {starting x}
  b := 1.0;;                 {ending x}
  h := 0.01;                 {starting step in x}
  ifmax := 20000;            {maximum calls of dgl() }
  itmax := 1000;             {maximum iterations}
  epsawp := 1e-7;            {absolute error}
  epsrb := 1e-8;             {relative error}
  wa := 0.0;                 {boundary condition in a}
  wb := 0.0;                 {boundary condition in b}

  Copy_vector(yanf, yanf0, N1);                   {yanf = yanf0}
  dgl(a, yanf, yanf0);  {dummy call of dgl() to initialize dgltxt in t_dgls}
  Copy_vector(yanf0, yanf, N1);   {restore yanf0}

  {output of the test example #2 for three varying methods}
  writeln;
  for method := 1 to 3 do
  begin
    if method > 1 then Copy_vector(yanf, yanf0, N1);     {yanf = yanf0}
    ClrScr;
    writeln('*   TEST EXAMPLE (WITH METHOD #',method,')');
    writeln('*   Method: ', TMethod[method]); 
    writeln('=============');
    writeln('ANALYZED SET OF DIFFERENTIAL EQUATIONS:');
    writeln('---------------------------------------');
    for j := 0 to N1-1 do
      writeln('*    ', dgltxt[j]);
    writeln('WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)');
    writeln('*     AT POINT  A = ',a:10:3,': ',wa:10:3);
    writeln('*     AT POINT  B = ',b:10:3,': ',wb:10:3);
    writeln('*   AT LOCATION  A = ',a:10:3,', THE FOLLOWING VALUES ARE PROVIDED:');
    for j := 0 to N1-1 do
      writeln('*     Y(',j+1,') = ', yanf^[j]);
    writeln('REQUIRED PARAMETER:');
    writeln('------------------:');
    writeln('*   -PRECISION FOR THE IVP = ', epsawp:20:10);
    writeln('*   -PRECISION FOR THE BVP = ', epsrb:20:10);
    writeln('*   -STEP WIDTH FOR THE IVP= ', h:20:10);
    writeln('*   -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP''S = ', ifmax:6);
    writeln('*   -MAX. NUMBER OF NEWTON-ITERATION STEPS = ', itmax:6);

    fehler := rwp(a, b, h, yanf, N1, method, epsawp, epsrb, ifmax, itmax, iter);
 
    {print results}
     if fehler <> 0 then
    begin
      writeln('*');
      Case fehler of
        1:  writeln('ERRORBOUND(S) TOO SMALL');
        2:  writeln('ERROR: B <= A');
        3:  writeln('ERROR: STEP SIZE H <= 0');
        4:  writeln('ERROR: NUMBER N OF DEQ''S NOT CORRECT');
        5:  writeln('ERROR: WRONG IVP-PROGRAM CHOICE');
        6:  writeln('IFMAX TO SMALL FOR SOLUTION OF THE IVP-PROGRAM');
        7:  writeln('AFTER ITMAX STEPS PRECISION WAS NOT REACHED!');
        8:  writeln('ERROR: JACOBI-MATRIX IS SINGULAR!');
        9:  writeln('LACK of MEMORY!')
      end
    end
    else                                             {no error? }
    begin
      {Output of solution}
      writeln('*');
      writeln('SOLUTION OF THE PROBLEM:');
      writeln('------------------------');
      writeln('THE START VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED');
      writeln('*   IN POINT A = ', a:10:3,' AS FOLLOWS:');
      for j := 0 to N1-1 do
        writeln('*     Y(',j+1,') = ', yanf^[j]:20:10);
      writeln('*');
      writeln('*   THIS REQUIRES ',iter,' NEWTON-ITERATIONS!');
      write('DETERMINATION COMPLETED AS PLANNED! ')
    end;
    ReadKey
  end;
  Dispose(yanf); Dispose(yanf0);
  DoneWinCrt
END.

{ --------------------------- END m_rwp.pas --------------------------