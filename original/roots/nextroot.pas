{****************************************************
*    Program to demonstrate NextRoot subroutine     *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                   *
*             Pascal version By J-P Moreau, Paris.  *
*                      (www.jpmoreau.fr)            *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
* (Example: Find another root of sin(x/2))          *
*                                                   *
* How many rootd have been determined: 3            *
*                                                   *
* Input the roots as prompted:                      *
*                                                   *
*     A( 1) = 0                                     *
*     A( 2) = 119.38052                             *
*     A( 3) = 75.398223                             *
*                                                   *
* Input the initial guess:                          *
*                                                   *
*     X0 = 75                                       *
*                                                   *
* Convergence criterion: .001                       *
*                                                   *
* Maximum number of iterations: 20                  *
*                                                   *
*                                                   *
* The calculated zero is X = 50.265482              *
*                                                   *
* The associated Y value is Y =  0.00021190         *
*                                                   *
* The number of steps was: 4                        *
*                                                   *
****************************************************}
PROGRAM NextRoot;
Uses WinCrt;

VAR
        A : Array[0..10] of DOUBLE;
        i,l,m,n : integer;
        e,x,x0,y1 : DOUBLE;



{********************************************
  Function subroutine                       }
  Function Y(x:DOUBLE;var y1:DOUBLE):DOUBLE;
  Begin
    Y := SIN(x/2.0);
    {derivative}
    y1 := 0.5*COS(x/2.0)
  End;
{*******************************************}


{***********************************************
*            NextRoot  subroutine              *
* -------------------------------------------- *
* This routine determines additional roots of  *
* a function Y(x), given a set of already      *
* established roots. Method applied is Newton- *
* Raphson iteration. The l established roots   *
* are in A(i). The routine requires an initial *
* guess, x0, and an accuracy criteria, e. It   *
* also requires a maximum number of iterations.*
***********************************************}
PROCEDURE Next_Root;
Label 100,200;
Var b,x1,yy:DOUBLE;
    i:integer; 
Begin
  n:=0;
  {Given x0, find Y/Y' }
100: yy:=Y(x0,y1);
  b:=y1/yy;
  for i:=1 to l do
    b:=b-1.0/(x0-A[i]);
  {Newton-Raphson iteration}
  x1:=x0-1.0/b;
  n:=n+1;
  {Test for convergence}
  if n>=m then goto 200;
  if ABS(x1-x0)<e then goto 200;
  x0:=x1;
  goto 100;
200: x:=x1
End;


{main program}
BEGIN
  clrscr;
  writeln;
  write(' How many roots have been determined: '); read(l);
  writeln;
  writeln(' Input the roots as prompted:');
  writeln;
  for i:=1 to l do
  begin
    write('   A(',i,') = '); read(A[i])
  end;
  writeln;
  writeln(' Input the initial guess:');
  writeln;
  write('    X0 = '); read(x0);
  writeln;
  write(' Convergence criterion: '); read(e);
  writeln;
  write(' Maximum number of iterations: '); read(m);

  Next_Root;  {Call NextRoot subroutine}

  writeln;
  writeln;
  writeln(' The calculated zero is X = ',x:9:6);
  writeln;                                  
  writeln(' The associated Y value is Y = ',Y(x,y1):9:6);
  writeln;
  writeln(' The number of steps was: ',n);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file nextroot.pas}