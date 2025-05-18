{**************************************************************************
*                    FRACTALS:   SET OF MANDELBROT                        *
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
* See also programs mandbrot.pas and julia.pas.                           *
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
  PROGRAM Mandelbrot;
  Uses WinCrtMy,Wintypes,Winprocs,Strings,Type_def,CrtGr2D;

  Const
        xmargin = 50;
        ymargin = 50;

  Var
        incX,incY,x1,x2,y1,y2 : real_ar;
        Lim,ligne,colonne,compt,sx : integer;
        p0,q0,module,x,y,aux: real_ar;
        sx1,sx2,sy1,sy2,slim : string[12];
        s: array[0..80] of char;
        CrtPen : HPen;  {pen identificator}

  BEGIN

    WinCrtInit('MANDELBROT');
    {define and select blue pen}
    CrtPen:=CreatePen(ps_Solid,1,RGB(0,0,255));
    SelectObject(CrtDC,CrtPen);
    {draw outside frame}
    Rectangle(CrtDC,25,25,MaxX-25,MaxY-80);  {all commands here are in screen pixels}
    x1:=-1.6; x2:=0.55; y1:=-1.15; y2:=-y1;
    Lim:=250;  {this parameter commands the finesse of the drawing}
    incX:=(x2-x1)/(MaxX-2*xmargin);
    incY:=(y2-y1)/(MaxY-2*ymargin);
    {prepare title}
    Str(x1:5:2,sx1); str(x2:5:2,sx2); str(y1:5:2,sy1); str(y2:5:2,sy2); Str(lim:3,slim);
    StrPCopy(s,'x1='+sx1+'  x2='+sx2+'  y1='+sy1+'  y2='+sy2+'  Lim='+slim);
    sx:=(MaxX-xmargin-6*strlen(s)) div 2;  {center title horizontaly}
    TextOut(CrtDC,sx,5,s,strlen(s));       {display title}

    {main loop to display Mandelbrot man line by line
     with a symmetry with respect to 0x axis }
    colonne:=1;
    While colonne <= MaxX-xmargin do
    begin
      p0:=x1+colonne*incX;
      ligne:=1;
      while ligne <= ((MaxY-ymargin) Div 2) do
      begin
        q0:=y1+ligne*incY;
        x:=0; y:=0;
        compt:=1; module:=0;
        while (compt <= lim) and (module < 4) do
        begin
          aux:=x;
          x:=sqr(x)-sqr(y)+p0;
          y:=2*y*aux+q0;
          module:=sqr(x)+sqr(y);
          Inc(compt)
        end;
        if module < 4 then
	begin
          MoveTo(CrtDC,xmargin+colonne,ymargin+ligne-25);
          LineTo(CrtDC,xmargin+colonne+1,ymargin+ligne-25);
          MoveTo(CrtDC,xmargin+colonne,MaxY-ymargin-ligne-25);
	  LineTo(CrtDC,xmargin+colonne+1,MaxY-ymargin-ligne-25)
	end;
        Ligne:=ligne+1
      end;
      Colonne:=colonne+1
    end;
    TextOut(CrtDC,40,40,'DOMAIN OF MANDELBROT',20);
    SortieGraphique;
    if rep='n' then DoneWinCrt

  END.

{end of file mandel.pas}