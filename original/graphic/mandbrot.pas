{**************************************************************************
*                 FRACTALS:   ZOOM ON THE SET OF MANDELBROT               *
* ----------------------------------------------------------------------- *
* TUTORIAL:                                                               *
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
*	iteration	new z	module of z                                       *
*	____________________________________________                          *
*	1		1 + 3i		3,16227...                                        *
*	2		-7 + 7i		9,89949...                                        *
*	3		1-97i		97,00515...                                       *
*   4      	-9407-193i	9408,97965                                        *
*                                                                         *
* The module of z increases very rapidly.                                 *
*                                                                         *
* 2)  c = -1 + 0,25 i                                                     *
*                                                                         *
*	iteration	new z	    module de z                                   *
*	____________________________________________                          *
*	1		-0,5   - 0i	    0,5                                           *
*	2		-0,25  + 0,5i	0,55902...                                    *
*	3		-0,687 + 0,25i	0,73154...                                    *
*                                                                         *
* At 80th iteration, z = 0.40868 + 0,27513 i, the module equals 0,49266...*
* The module remains with a finite value, whatever the number of itera-   *
* tions may be. Practically, we will consider that the limit is obtained  *
* after 100 iterations, if the module of z is < 4. The set of Mandelbrot, *
* still called the Mandelbrot man, because of its shape, is constituted   *
* by all the points for which  the expression z² + c has a finite value   *
* whatever the number of iterations may be.                               *
*                                                                         *
* This program allows to make a zoom on a particular zone of the Mandel-  *
* brot domain by fixing the following parameters:                         *
*                                                                         *
* 	x1 et x2		limits in ox                                          *
*	y1 et y2		limits in oy                                          *
*	Lim		        number of iterations (usually 100)                    *
*                                                                         *
* Some interesting examples:                                              *
*                                                                         *
*	x1 = -0.67166		x2 = -0.44953                                     *
*	y1 = 0.49216		y2 = 0.71429		Lim = 250                     *
*                                                                         *
*	 x1 = -0.19050		x2 = -0.13824                                     *
*	y1 = 1.01480		y2 = 1.06707		Lim = 320                     *
*                                                                         *
*       x1 = -0.74591		x2 = -0.74448                                 *
*	y1 = 0.11196		y2 = 0.11339		Lim = 250                     *
*                                                                         *
* These three examples are quite acceptable with Lim = 100.               *
*                                                                         *
* Up to you to discover other ones !                                      *
*                                                                         *
*                                                                         *
* See also programs mandel.pas and julia.pas.                             *
*                                                                         *
* ----------------------------------------------------------------------- *
* REFERENCE:                                                              *
*     "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0      *
*      By R. Dony - MASSON, Paris 1990" [BIBLI 12].                       *
*                                                                         *
*                                          TPW version by J-P Moreau      *
*                                              (www.jpmoreau.fr)          *
* ----------------------------------------------------------------------- *
* Note: the print option is not implemented here.                         *
**************************************************************************}
Program mandelbrot;   {WINCRTMY color version without printing}
Uses WinCrtMy,WinTypes,WinProcs,WObjects,Strings,Type_Def,CrtGr2D;

Var x1,x2,y1,y2,incX,incY,p,q: REAL_AR;
    indice,limite:INTEGER;
    compt: integer;
    aborted: boolean;
    p0,q0,module,x,y,aux,xx,yy: REAL_AR;
    CrtPen,CrtPen1,CrtPen2,CrtPen3,OldPen: HPen;
    

    Procedure Donnees;
    begin
      clrscr;
      writeln;
      writeln;
      write('  Interval X1, X2: '); readln(x1,x2);
      writeln;
      write('  Interval Y1, Y2: '); readln(y1,y2);
      writeln;
      write('  Number of iterations: '); readln(limite);
      if (X2<=X1) OR (Y2<=Y1) OR (limite<10) then
      begin
        MessageBox(CrtWindow,'Wrong data !',
                   'ERROR',mb_Ok);
        aborted:=TRUE
      end;
      IncX:=(x2-x1)/MaxX;
      IncY:=(y2-y1)/MaxY;
      clrscr;
      gotoxy(10,1);
      write('X1=',x1:7:4,'  X2=',x2:7:4,'  Y1=',y1:7:4,'  Y2=',y2:7:4,'  Lim=',limite:3);
    end;

    
 BEGIN
   WinCrtInit('ZOOM on Mandelbrot');
   Repeat
      Donnees;
      aborted:=false;
      CrtPen:=CreatePen(ps_Solid,1,RGB(0,0,255));
      CrtPen1:=CreatePen(ps_Solid,1,RGB(0,255,0));
      CrtPen2:=CreatePen(ps_Solid,1,RGB(255,0,0));
      CrtPen3:=CreatePen(ps_Solid,1,RGB(255,255,255));
      OldPen:=SelectObject(CrtDC,CrtPen);
      Fenetre(x1,x2,y1,y2);
      Cloture(15,MaxX-25,80,MaxY-15);
      Bordure(CrtDC);
      xx:=x1;
      while (xx <= x2) and (not aborted) do
      begin
        yy:=y1;
        while yy <= y2 do
        begin
          compt:=0;
          module:=0;
          p:=xx; q:=yy; x:=0; y:=0;
          while (compt <= limite) and (module <= 4) do
          begin
            aux:=x;
            x:=x*x-y*y+p;
            y:=2*y*aux+q;
            module:=x*x+y*y;
            inc(compt)
          end;
          if module > 4 then
          begin
            indice:=compt mod 5;
            case indice of
              0: SelectObject(CrtDC,OldPen);
              1: SelectObject(CrtDC,CrtPen);
              2: SelectObject(CrtDC,CrtPen1);
              3: SelectObject(CrtDC,CrtPen2);
              4: SelectObject(CrtDC,CrtPen3);
            end;  
            MoveXY(CrtDc,-xx,-yy);
            LineXY(CrtDc,-xx+incX,-yy);
            MoveXY(CrtDc,xx,yy);
            LineXY(CrtDc,xx+incX,yy);
          end;
          yy:=yy+incY
        end;
        xx:=xx+incX;
        if Keypressed then aborted:=true
      end;
      DeleteObject(CrtPen);
      DeleteObject(CrtPen1);
      DeleteObject(CrtPen2);
      DeleteObject(CrtPen3);
      SortieGraphique;
      if rep='i' then MessageBox(CrtWindow,'Print option not implemented !',
                   'WARNING',mb_Ok);
   Until rep='n';
   DoneWinCrt
 END.

{end of file mandbrot.pas}