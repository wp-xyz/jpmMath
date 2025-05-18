{utility functions for Fourier series}
Unit Fourier;

Interface

CONST NSIZE = 201;

TYPE
     pTab = ^Table;    {used by discreet functions}
     Table = Array[1..NSIZE] of Double;

VAR  flag: BOOLEAN;    {FALSE for an, TRUE for bn}
     om,T: Double;     {om=2*pi*n/T}


  Procedure DiscreetFourierHn(ndata:Integer; VAR X,Y:pTab; n:Integer; VAR a,b: Double);
  Procedure AnalyticFourierHn(t1,t2:Double; n:Integer; VAR a,b: Double);


Implementation

  {analytical function to study (1 period from x=-pi to x=pi):
   ---------------------------------------------------------------
  Note: Exact Fourier coefficients for this periodic function are:
    an=0 if n is even, -1/2n if n=4p+1, 1/2n if n=4p+3
    bn=0 if n=4p, -1/2n if n=4p+1, 1/n if n=4p+2, -1/2n if n=4p+3
------------------------------------------------------------------
  Function F(x:Double):Double;
  begin
    if x<-PI then f:=0.0
    else if x<-PI/2 then f:=PI/4
    else if x<=PI then f:=-PI/4
    else f:=0.0
  End;           }

{ Note: Exact Fourier coefficients for this periodic function are:
    a0=1/pi, an = 0 if n is odd, an = -2/(pi*(n2-1)) if n is even;
    b1=1/2,  bn = 0 for n>=2.
------------------------------------------------------------------}
  Function F1(x:Double):Double;
  begin
    if x < 0 then F1:=0.0
    else if x < PI then F1:=sin(x)
    else F1:=0
  end;

  {Function to integrate by Romberg method}
  Function FUNC(x:Double):Double;
  Begin
    if Not flag then FUNC:=F1(x)*cos(om*x)  {for an}
                else FUNC:=F1(x)*sin(om*x)  {for bn}
  End;

  {***************************************************************
  *   Calculate the Fourier harmonic #n of a periodic discreet   *
  *   function F(x) defined by ndata points.                     *
  * ------------------------------------------------------------ *
  * Inputs:                                                      *
  *            ndata: number of points of discreet function.     *
  *            X    : pointer to table storing xi abscissas.     *
  *            Y    : pointer to table storing yi ordinates.     *
  *                                                              *
  * Outputs:                                                     *
  *            a    : coefficient an of the Fourier series.      *
  *            b:   : coefficient bn of the Fourier series.      *
  * ------------------------------------------------------------ *
  * Reference: "Mathematiques en Turbo-Pascal By Marc Ducamp and *
  *             Alain Reverchon, 1. Analyse, Editions Eyrolles,  *
  *             Paris, 1991" [BIBLI 03].                         *
  *                                                              *
  *                           TPW version by J-P Moreau, Paris.  *
  ****************************************************************
  Note: The Fourier series  of a periodic discreet function F(x)
        can be written under the form:
                    n=inf.
        F(x) = a0 + Sum ( an cos(2 n pi/T x) + bn sin(2 n pi/T x)  
                    n=1
        (b0 is always zero).
  ---------------------------------------------------------------}  
  Procedure DiscreetFourierHn(ndata:Integer; VAR X,Y:pTab; n:Integer; VAR a,b: Double);
  Var i: Integer;
      xi,wa,wb,wc,wd,wg,wh,wi,wl,wm,wn,wp: Double;
  Begin
    T:=X^[ndata]-X^[1]; xi:=(X^[ndata]+X^[1])/2;
    om := 2*PI*n/T; a:=0; b:=0;
    for i:=1 to ndata-1 do
    begin
      wa:=X^[i]; wb:=X^[i+1];
      wc:=Y^[i]; wd:=Y^[i+1];
      if wa<>wb then
      begin
        wg := (wd-wc)/(wb-wa);
        wh := om*(wa-xi); wi:=om*(wb-xi);
        if n=0 then
          a := a+(wb-wa)*(wc+wg/2*(wb-wa))
        else
        begin
          wl := cos(wh); wm := sin(wh);
          wn := cos(wi); wp := sin(wi);
          a := a + wg/om*(wn-wl) + wd*wp - wc*wm;
          b := b + wg/om*(wp-wm) - wd*wn + wc*wl
        end
      end
    end;
    a := a/T; b := b/T;
    if n<>0  then
    begin
      a := a * 2/om; b := b * 2/om
    end
  End;

  {******************************************************
  * Integral of a function FUNC(X) by Romberg's method  *
  * --------------------------------------------------- *
  * INPUTS:                                             *
  *          a      begin value of x variable           *
  *          b      end value of x variable             *
  *       prec      desired precision                   *
  *    itermin      minimum number of iterations        *
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

  {***************************************************************
  * Calculate the Fourier harmonic #n of a periodic function F(x)*
  * analytically defined.                                        *
  * ------------------------------------------------------------ *
  * Inputs:                                                      *
  *            t1   : -period/2                                  *
  *            t2   : period/2                                   *
  *            n    : order of harmonic                          *          
  *                                                              *
  * Outputs:                                                     *
  *            a    : coefficient an of the Fourier series.      *
  *            b:   : coefficient bn of the Fourier series.      *
  * ------------------------------------------------------------ *
  * Reference: "Mathematiques en Turbo-Pascal By Marc Ducamp and *
  *             Alain Reverchon, 1. Analyse, Editions Eyrolles,  *
  *             Paris, 1991" [BIBLI 03].                         *
  *                                                              *
  *                           TPW version By J-P Moreau, Paris.  *
  ****************************************************************
  Note:  When a periodic function f(x), of period T, can be developped
         in a Fourier series, this one is unique:
                     n=inf.
         F(x) = a0 + Sum ( an cos(2 n pi/T x) + bn sin(2 n pi/T x)  
                     n=1
         (b0 is always zero).

         The coefficients an, bn are given by:
                      T/2
             a0 = 1/T Sum (f(x) dx)
                      -T/2
                      T/2
             an = 2/T Sum (f(x) cos(2npi/T x) dx)
                      -T/2
                      T/2
             an = 2/T Sum (f(x) sin(2npi/T x) dx)
                      -T/2
         Here, the integrals are calculated by the Romberg method.
  ---------------------------------------------------------------------}
  Procedure AnalyticFourierHn(t1,t2:Double; n:Integer; VAR a,b: Double);
  Var precision,temp: Double;
      itmax,niter:Integer;
  Begin
    T:=ABS(t2-t1);
    om:=2.0*PI*n/T;
    precision:=1e-10; itmax:=15;
    flag:=FALSE; {calculate an}
    a:=RombergIntegral(t1,t2,precision,temp,niter,5,itmax)*2/T;
    flag:=TRUE; {calculate bn}
    b:=RombergIntegral(t1,t2,precision,temp,niter,5,itmax)*2/T;
    if n=0 then a:=a/2;
  end;

END.

{end of file fourier.pas (last modified 06/11/2002) }