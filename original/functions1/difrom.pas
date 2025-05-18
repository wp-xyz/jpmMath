{ ------------------------ UNIT difrom.pas ----------------------------
!     "Numerical Algorithms with C, By Gisela Engeln-Muellges         !
!        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11]           !
!                                                                     !
!           Pascal Version By Jean-Pierre Moreau, Paris.              !
!                          (www.jpmoreau.fr)
----------------------------------------------------------------------}
UNIT DIFROM;

Interface

  Const NMAX = 10;  MACH_EPS = 1.2e-16;

 Procedure difrom1(x0, eps: Double;
                   n: Integer;
                   h: Double;
               Var res, er_app: Double;
               Var nend: Integer;
               Var hend: Double;
               Var error:Integer);

Implementation

{test function to derivate}
Function func(x: Double): Double;
begin
  func := 1 / x
end;

Procedure difrom1(x0, eps: Double;
                  n: Integer;
                  h: Double;
              Var res, er_app: Double;
              Var nend: Integer;
              Var hend: Double;
              Var error:Integer);
{************************************************************************
*  Computes an approximation for the first derivative of func at x0     *
*  using the ROMBERG method.                                            *
* --------------------------------------------------------------------- *
*  Input parameters :                                                   *
*                                                                       *
*    Double func(Double)  (External) Name of function to be             *
*                         differentiated.                               *
*    Double x0            value of abscissa at which derivative is to   *
*                         be found.                                     *
*    Double eps           desired accuracy.                             *
*    int    n             max. number of columns in the Romberg scheme  *
*                         (n > 1).                                      *
*    Double h             initial step size.                            *
*                                                                       *
*  Output parameters:                                                   *
*                                                                       *
*    Double res           approximate derivative                        *
*    Double er_app        error estimate for  res                       *
*    int    nend          number of columns actually used in scheme     *
*    Double hend          final step size                               *
*                                                                       *
*  Return value :                                                       *
*                                                                       *
*    0:   no error: er_app < eps                                        *
*    1:   n < 1  or  eps <= 0  or  h < MACH_EPS                         *
*    2:   desired accuracy not reached after n steps                    *
*    3:   step size drooped below  MACH_EPS                             *
*    4:   lack of sufficient memory (not used here)                     *
*                                                                       *
************************************************************************}
Label  10, return;
Var
  i, j, m   :Integer;
  h2, d1, d2: Double;
  d: Array[0..NMAX] of Double;

Begin

  if (n <= 1) or (eps <= 0.0) or (h < MACH_EPS) then    {check input}
  begin
    error:=1;
    goto return
  end;

  h2 := 2.0 * h;
  d[0] := (func(x0 + h) - func(x0 - h)) / h2;

{***********************************************************************
* This loop runs until the maximum of Romberg rows is filled or until  *
* the step size use drops below the machine constant or if the desired *
* accuracy is reached.                                                 *
***********************************************************************}
  error:=2;
  for j := 1 to n - 1 do
  begin
    d[j] := 0.0;
    d1   := d[0];
    h2   := h;
    h    := h * 0.5;
    nend := j;

    if h < MACH_EPS then
    begin
      error := 3;                 {step size less than machine constant}
      goto 10
    end;

    d[0] := (func(x0 + h) - func(x0 - h)) / h2;

    m:=4;
    for i := 1 to j do
    begin
      d2 := d[i];
      d [i] := (m * d[i-1] - d1) / (m-1);
      d1 := d2;
      m := m * 4
    end;

    er_app := ABS (d[j] - d[j-1]);

    if er_app < eps then
    begin
      error := 0;                         {desired accuracy reached}
      goto 10
    end
  end;

10:res := d[nend];                        {save final values}
   hend := h;

Return: End;

END.

{ ------------------------- END difrom.pas ---------------------------- }