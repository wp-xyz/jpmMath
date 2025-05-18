{**************************************************************************
*                      FRACTALS:   SET OF JULIA                           *
* ----------------------------------------------------------------------- *
* TUTORIAL:                                                               *
*                                                                         *
* PRELIMINARY: THE SET OF MANDELBROT (see program mandel.pas)             *
*                                                                         *
* The set of Mandelbrot is a figure in the complex plane, that is (for no *
* mathematicians) the plane allowing to represent complex numbers having  *
* the form a + i b, a is the real part (ox componant) and b is the        *
* imaginary part (oy componant). The imaginary number i is the unity of   *
* axis oy, defined as i x i  = -1. In the following, we will use the      *
* multiplication of two complex numbers and the module value of a complex *
* number before defining the set of Mandelbrot itself.                    *
*                                                                         *
* Multiply z1 by z2                                                       *
*                                                                         *
* Let be z1 = a1 + i b1 and z2 = a2 + i b2, two complex numbers. To cal-  *
* culate Z = z1 x z2, one just has to algebraically develop the expression*
* (a1 + i b1)*(a2 + i b2), replacing i² by -1 and sepatating constant     *
* terms from i terms. One easily obtains:                                 *
*                                                                         *
*          Z = (a² - b²) + i (2ab)                                        *
*                                                                         *
* Module of a complex number                                              *
*                                                                         *
* Every complex number a + i b can be represented as the point of coordi- *
* nates (a, b) in the complex plane. Its module is nothing other than the *
* distance of this point from origin (0,0), that is to say:               *
*                                                                         *
*        | Z | = square root of (a² + b²)                                 *
*                                                                         *
* We can now define the iterative process that will lead us to the magni- *
* ficent "Set of Mandelbrot":                                             *
*                                                                         *
* We take, as a starting point, a fixed complex number c, and we calculate*
* the complex expression z² + c, z being a variable complex number.       *
* Let us take z = 0, then z² + c equals c. Let us replace z by c in the   *
* formula z² + c, then we obtain c² + c. Let us again remplace z by this  *
* new value, then we obtain (c²+c)² + c. Such a process is indefinitly    *
* resumed, each new value z² + c is taken as new value z. This provides   *
* an unlimited series of z numbers.                                       *
*                                                                         *
* The mathematician Mandelbrot was the first to notice that, according to *
* the chosen value c, the series of complex numbers z so obtained could   *
* have a module, either rapidly tending towards infinity, or tending      *
* towards a finite value, whatever the number of iterations may be.       * 
*                                                                         *
* Two examples:                                                           *
*                                                                         *
* 1)  c = 1 + i                                                           *
*                                                                         *
*	iteration	new z	        module of z                       *
*	____________________________________________                      *
*	1		1 + 3i		3,16227...                        *
*	2		-7 + 7i		9,89949...                        *
*	3		1-97i		97,00515...                       *
*       4       	-9407-193i	9408,97965                        *
*                                                                         *
* The module of z increases very rapidly.                                 *
*                                                                         *
* 2)  c = -1 + 0,25 i                                                     *
*                                                                         *
*	iteration	new z	        module de z                       *
*	____________________________________________                      *
*	1		-0,5 - 0i	0,5                               *
*	2		-0,25 + 0,5i	0,55902...                        *
*	3		-0,687 + 0,25i	0,73154...                        *
*                                                                         *
* At 80th iteration, z = 0.40868 + 0,27513 i, the module equals 0,49266...*
* The module remains with a finite value, whatever the number of itera-   *
* tions may be. Practically, we will consider that the limit is obtained  *
* after 100 iterations, if the module of z is < 4. The set of Mandelbrot, *
* still called the Mandelbrot man, because of its shape, is constituted   *
* by all the points for which  the expression z² + c has a finite value   *
* whatever the number of iterations may be.                               *
*                                                                         *
*                                                                         *
* SET OF JULIA                                                            *
*                                                                         *
* If you take as a rule to fix the value of complex number c and to make  *
* vary the  starting value of complex number z, you define a new set,     *
* different from the previous set of Mandelbrot, called Set of Julia,     *
* after the name of a french mathematician who studied it forst. Each new *
* c value injected into the iterative formula z² + c generates a new set  *
* of Julia which  are in unlimited amount. These sets of Julia can have   *
* amazingly different forms according to the starting value.              *
*                                                                         *
* If the starting point c is chosen inside the main body of the Mandelbrot*
* man, then the corresponding set of Julia is in one single part.         * 
* According to other possible starting points, the set of Julia can have  *
* the form of an infinity of so called fractal circles, or a "dendritic"  *
* line (with many ramifications), or can be composed of several parts.    *
* At the limit, if c is far away from the Mandelbrot man's body, the set  *
* of Julia is reduced to a "dust" of points, called "Fatou's dust" after  *
* the name of another french mathematicians.                              *
*                                                                         *
* The program asks for the following parameters:                          *
*                                                                         *
*     x1, x2: interval for ox axis                                        *
*     y1, y2: interval for oy axis                                        *
*     Lim   : maximum number of iterations                                *
*     p, q  : coordinates of starting point c                             *
*                                                                         *
* Some interesting examples:                                              *
*                                                                         *
*       (Dragon shape)                                                    *
* 	x1 = -1.5		x2 = 1.5                                              *
*	y1 = -1.5		y2 = 1.5		Lim = 200                             *
*	p  = 0.32		q  = 0.043                                            *
*                                                                         *
*       (Jewel shape)                                                     * 
* 	x1 = -1.8		x2 = 1.8                                              *
*	y1 = -1.8		y2 = 1.8		Lim = 200                             *
*	p  = -0.74543		q  = 0.11301                                      *
*                                                                         *
*       (cactus shape)                                                    *
* 	x1 = -1.4		x2 = 1.4                                              *
*	y1 = -1.4		y2 = 1.4		Lim = 200                             *
*       p  = -0.12		q  = 0.74                                         *
*                                                                         *
* Up to you to discover other ones !                                      *
* ----------------------------------------------------------------------- *
* REFERENCE:                                                              *
*     "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0      *
*      By R. Dony - MASSON, Paris 1990" [BIBLI 12].                       *
*                                                                         *
*                                          TPW version by J-P Moreau      *
*                                              (www.jpmoreau.fr)          *
**************************************************************************}
Program Set_of_julia;
Uses WinCrtMy,WinTypes,WinProcs,WObjects,Strings,Type_Def,CrtGr2D;

 Var x1,x2,y1,y2,incX,incY,p,q: real_ar;
     indice,limite:byte;
     CrtPen : HPen; {blue pen}

    Procedure Data;
    begin
      WinCrtInit('JULIA');
      clrscr;
      writeln;
      writeln('       SET OF JULIA');
      writeln('       ------------');
      writeln; writeln;
      write('  Interval X1, X2: '); readln(x1,x2);
      writeln;
      write('  Interval Y1, Y2: '); readln(y1,y2);
      writeln;
      write('  Number of iterations: '); readln(limite);
      writeln;
      write('  Starting values (Xc, Yc): '); readln(p,q);
      incX:=(x2-x1)/(MaxX*1.0);
      incY:=(y2-y1)/(MaxY*1.0);
    end;


  PROCEDURE Draw_Julia;    {draw a set of Julia}
  VAR compt: INTEGER;
  module,x,y,aux,xx,yy: REAL_AR;
  BEGIN
    CrtPen:=CreatePen(ps_Solid,1,RGB(0,0,255));
    SelectObject(CrtDC,CrtPen);
    Clrscr;
    Fenetre(x1,x2,y1,y2);
    Cloture(10,MaxX-25,80,MaxY-10);
    Bordure(CrtDc);
    xx:=x1;
    WHILE xx <= x2 DO
    BEGIN
      yy:=y1;
      WHILE yy <= 0 DO
      BEGIN
        compt:=0;
        module:=0;
        x:=xx; y:=yy;
        WHILE (compt <= limite) AND (module <= 4) DO
        BEGIN
          aux:=x;
          x:=x*x-y*y+p;
          y:=2*y*aux+q;
          module:=x*x+y*y;
          inc(compt)
        END;
        IF module <= 4 THEN
        BEGIN
          MoveXY(CrtDC,-xx,-yy);
          LineXY(CrtDC,-xx+incX,-yy);
          MoveXY(CrtDC,xx,yy);
          LineXY(CrtDC,xx+incX,yy);
        END;
        yy:=yy+incY
      END;
      xx:=xx+incX
    END;
    TextOut(CrtDC,25,15,'SET OF JULIA',12);
    SortieGraphique
  END;

  {main program}
  Begin
    rep:='o';
    Repeat
      Data;
      Draw_Julia
    Until rep='n';
    DeleteObject(CrtPen);
    DoneWinCrt
  End.

{end of file julia.pas}