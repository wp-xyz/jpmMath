{ ------------------------ UNIT fseidel.pas ------------------------- }
UNIT FSEIDEL;

{BIBLI 11}

INTERFACE

Const
       SIZE = 25;                       {Maximum size of linear system}
       ITERMAX = 300;                    {Maximal number of iterations}

       ONE  = 1.0;
       TWO  = 2.0;
       ZERO = 0.0;

Type
       pMAT = ^MATR;
       MATR = Array[1..SIZE,1..SIZE] of Double;
       pVEC = ^VEC;
       VEC =  Array[1..SIZE] of Double;

       Procedure seidel           {Gauss Seidel Method with relaxation}
           (
            crit: integer;         {crit = 0, 1, 2, 3 ................}
            n: integer;            {size of matrix ...................}
            mat: pMAT;             {pointer to matrix ................}
            b: pVEC;               {pointer to right hand side .......}
            omega: double;         {relaxaktion coefficient ..........}
            var x: pVEC;           {solution vector...................}
            var residu: pVEC;      {residuum vector ..................}
            var iter: integer;     {# of iterations ..................}
            var rc: integer        {return code.......................}
           );

IMPLEMENTATION

 {====================================================================*
 *                                                                    *
 *  seidel solves the linear system  mat * x = b  iteratively.        *
 *  Here  mat  is a nonsingular  n x n  matrix, b is the right hand   *
 *  side for the linear system and x is the solution.                 *
 *                                                                    *
 *  seidel uses the Gauss Seidel Method with relaxation for a given   *
 *  relaxation coefficient 0 < omega < 2.                             *
 *  If  omega = 1, the standard Gauss Seidel method (without          *
 *  relaxation) is performed.                                         *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Applications:                                                    *
 *   =============                                                    *
 *      Solve linear systems with nonsingular system matrices that    *
 *      satisfy one of the following criteria: row sum criterion,     *
 *      column sum criterion or the criterion of Schmidt and v. Mises.*
 *      Only if at least one of these criteria is satisfied for mat,  *
 *      convergence of the scheme is guaranteed.                      *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      crit     int crit;                                            *
 *               select criterion                                     *
 *               =1 : row sum criterion                               *
 *               =2 : column sum criterion                            *
 *               =3 : criterion of Schmidt-v.Mises                    *
 *               other : no check                                     *
 *      n        int n;  ( n > 0 )                                    *
 *               size of mat, b and x                                 *
 *      mat      REAL    mat[n];                                      *
 *               Matrix of the liear system                           *
 *      b        REAL    b[n];                                        *
 *               Right hand side                                      *
 *      omega    REAL   omega; ( 0.0 < omega < 2.0 )                  *
 *               Relaxation coefficient.                              *
 *      x        REAL   x[n];                                         *
 *               Starting vector for iteration                        *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      x        REAL   x[n];                                         *
 *               solution vector                                      *
 *      residu   REAL   residu[n];                                    *
 *               residual vector  b - mat * x; close to zero vector   *
 *      iter     int *iter;                                           *
 *               Number of iterations performed                       *
 *                                                                    *
 *   Return code rc:                                                  *
 *   ==============                                                   *
 *      =  0     solution has been found                              *
 *      =  1     n < 1  or omega <= 0 or omega >= 2                   *
 *      =  2     improper mat or b or x (not used here)               *
 *      =  3     one diagonal element of mat vanishes                 *
 *      =  4     Iteration number exceeded                            *
 *      = 11     column sum criterion violated                        *
 *      = 12     row sum criterion violated                           *
 *      = 13     Schmidt-v.Mises criterion violated                   *
 *                                                                    *
 *====================================================================}
Procedure Seidel;
Label 10, 20, return;
Var
  i, j: integer;
  tmp, eps: double;

Begin

  rc:=0;
  iter := 0;                             {Initialize iteration counter}

  if (n<1) or (omega<=ZERO) or (omega>=TWO) then          {Check omega}
  begin
    rc:=1;
    goto return
  end;

  eps := 1e-8;

  for i := 1 to n do                        {transform mat so that all}
  begin                                     {diagonals equal 1.       } 
    if mat^[i,i] = ZERO then
    begin
      rc := 3;
      goto return
    end;
    tmp := ONE / mat^[i,i];
    for j := 1 to n do mat^[i,j] := mat^[i,j] * tmp;
    b^[i] := b^[i] * tmp;                    {adjust right hand side b}
  end;

  Case crit of                             {check convergence criteria}
  
       1: for i := 1 to n do                        {row sum criterion}
          begin
            tmp:=ZERO;
            for j := 1 to n do
              tmp := tmp + ABS(mat^[i,j]);
            if tmp >= TWO then
            begin
              rc:=11;
              goto return
            end
          end;

       2: for j:=1 to n do                       {column sum criterion}
          begin
            tmp:=ZERO;
            for i:=1 to n do
              tmp := tmp + ABS(mat^[i,j]);
            if tmp >= TWO then
            begin
              rc:=12;
              goto return
            end
          end;

       3: begin
          tmp:=ZERO;
          for i:=1 to n do
            for j:=1 to n do                   {criterion of Schmidt, }
              tmp := tmp  + SQR(mat^[i,j]);    {Von Mises.            }
            tmp := SQRT(tmp - ONE);
            if tmp >= ONE then
            begin
              rc:=13;
              goto return
            end
          end
  end;

  for i:=1 to n do residu^[i] := x^[i];             {store x in residu}

  while iter <= ITERMAX do                            {Begin iteration}
  begin

    Inc(iter);

    for i:=1 to n do
    begin
      tmp:=b^[i];
      for j:=1 to n do
        tmp := tmp - mat^[i,j] * residu^[j];
      residu^[i] := residu^[i] + omega * tmp
    end;

    for i:=1 to n do                        {check break-off criterion}
    begin
      tmp := x^[i] - residu^[i];

      if ABS(tmp) <= eps then
      begin
        x^[i] := residu^[i];                {if rc = 0 at end of loop }
        rc := 0;                            { -> stop iteration.      }
      end
      else
      begin
        for j:=1 to n do x^[j] := residu^[j];
        rc := 4;
        goto 10
      end
    end;
    if rc = 0 then goto 20;                           {solution found }
10:end;                                                {End iteration }

20:for i := 1 to n do                           {find residual vector }
  begin
    tmp := b^[i];
    for j:=1 to n do
      tmp := tmp - mat^[i,j] * x^[j];
    residu^[i] := tmp
  end;

return: end;

END.

{ ------------------------ END fseidel.pas --------------------------- }