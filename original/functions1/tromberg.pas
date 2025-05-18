{*************************************************************
* Program to demonstrate integration of a real function F(x) *
*                    by Romberg's method                     *
* ---------------------------------------------------------- *
* REFERENCE: "Mathematiques en Turbo-Pascal (Part 1) by      *
*             Marc Ducamp and Alain Reverchon, Eyrolles,     *
*             Paris, 1991" [BIBLI 03].                       *
* ---------------------------------------------------------- *
* SAMPLE RUN:                                                *
* (Integrate sin(x) from x=0 to x=1)                         *
*                                                            *
* Integral of a function F(X) by Romberg's method            *
*                                                            *
* Input begin and end values for x variable:                 *
*                                                            *
*         X0 = 0                                             *
*         X1 = 1                                             *
*                                                            *
* Input desired precision: 1e-10                             *
*                                                            *
*                                                            *
* Value of integral   :  4.59697694131851E-0001              *
*                                                            *
* Obtained precision  :  9.59908263986620E-0011              *
*                                                            *
* Number of iterations:  4                                   *
*                                                            *
*************************************************************}
Program Test_Romberg;
Uses WinCrt;

Var
        x0      : double;     {begin x value}
        x1      : double;     {end   x value}
        prec    : double;     {desired precision}
        integral: double;     {result of integral}
        obtprec : double;     {obtained precision}

        niter   : integer;    {number of actual iterations}


  {Given function to integrate}
  Function FUNC(x:double): double;
  Begin
    FUNC := SIN(x)
  End;


  {******************************************************
  * Integral of a function FUNC(X) by Romberg's method *
  * --------------------------------------------------- *
  * INPUTS:                                             *
  *          a      begin value of x variable           *
  *          b      end value of x variable             *
  *       prec      desired precision                   *
  *    itermin      minimum number of ietrations        *
  *    itermax      maximum number of iterations        *
  *                                                     *
  * OUTPUTS:                                            *
  *    obtprec      obtained precision for integral     *
  *          n      number of iterations done           *
  *                                                     *
  * RETURNED VALUE  the integral of FUNC(X) from a to b *
  *                                                     *
  ******************************************************} 
  Function RombergIntegral(a,b,prec:double; VAR obtprec:double;
                               VAR n:integer; itermin,itermax:integer) : double;
  Const MAXITER = 15;
  Var   i,j        : integer;
        pas,r,s,ta : double;
        t          : array[0..MAXITER,0..MAXITER] of double;
  Begin
    if itermax>MAXITER then itermax:=MAXITER;
    r := FUNC(a);
    ta := (r + FUNC(b)) / 2;
    n:=0;
    pas:=b-a;
    t[0,0]:=ta*pas;
    Repeat
      Inc(n);
      pas:=pas/2;
      s:=ta;
      for i:=1 to pred(1 SHL n) do   {2^n-1}
        s:=s + FUNC(a+pas*i);
      t[0,n]:=s*pas;
      r:=1;
      for i:=1 to n do
      begin
        r:=r*4;
        j:=n-i;
        t[i,j]:=(r*t[i-1,j+1] - t[i-1,j])/(r-1)
      end;
      obtprec := ABS(t[n,0] - t[n-1,0])
    Until (n>=itermax) OR ((obtprec<prec) AND (n>=itermin));
    RombergIntegral := t[n,0]
  End;

          
{main program}
BEGIN
  writeln;
  writeln(' Integral of a function F(X) by Romberg''s method');
  writeln;
  writeln(' Input begin and end values for x variable:');
  writeln;
  write('         X0 = '); read(x0);
  write('         X1 = '); read(x1);
  writeln;
  write(' Input desired precision: '); read(prec);

  integral := RombergIntegral(x0,x1,prec,obtprec,niter,1,50);

  writeln;
  writeln;
  writeln(' Value of integral   : ', integral);
  writeln;
  writeln(' Obtained precision  : ', obtprec);
  writeln;
  writeln(' Number of iterations:  ', niter);
  ReadKey; DoneWinCrt

END.

{End of file tromberg.pas}