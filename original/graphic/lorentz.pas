{**************************************************************************
*                  L O R E N T Z   A T T R A C T O R                      * 
* ----------------------------------------------------------------------- *
* TUTORIAL:                                                               *
*                                                                         *
* Until recently, the only known attractors were the fixed point, the     *
* limit cycle and the torus. In 1963, Edwards Lorenz, meteorologist at the*
* M.I.T.,discovered a practical example of a simple dynamic system presen-*
* ting a complex behaviour. To adapt them to computers available at that  *
* time, he began with simplifying the equations of metéorology to obtain  *
* finally a model composed of three differential equations with three     *
* unknown variables x, y, z and three parameters a, b, c:                 *
*                                                                         *
*       dx / dt = - a x + a y                                             *
*		dy / dt = b x - y - x z                                           *
*		dz / dt = -c z + x y                                              *
*                                                                         *
* During very long simulations on a computer, Lorentz decided, to check a *
* result, to restart the same calculation halfway in order to spare time. *
* For that, he reinjected into the computer the intermediate date obtained*
* earlier. He was very surprised to see that the new results were comple- *
* tely different from the first ones. After he suspected some computer    *
* failure, Lorenz understood at last that the big difference between both *
* solutions came from very small differences in data. These small pertur- *
* bations exponentially amplified themselves, being doubled every four    *
* days in simulated time, so after two months the results became entirely *
* different !                                                             *
*                                                                         *
* Lorenz then realized that it would be very difficult to make meteorolo- *
* gical foresights in the long term, the slightest change in initial con- *
* ditions leading to a radically different evolution of atmosphere.       *
*                                                                         *
* This is still the case today with atmospheric models much more sophis-  *
* ticated.                                                                *
*                                                                         *
* None of the three attractors known at the time could predict the beha-  *
* viour of such a dynamic system. Lorenz had just discovered a strange or *
* "chaotic" attractor to which his name was given.                        *
* ----------------------------------------------------------------------- *
* REFERENCE:                                                              *
*     "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0      *
*      By R. Dony - MASSON, Paris 1990" [BIBLI 12].                       *
*                                                                         *
*                                          TPW version by J-P Moreau      *
*                                              (www.jpmoreau.fr)          *
**************************************************************************}
Program Lorentz_Attractor;
Uses WinCrtMy,WinTypes,WinProcs,Type_def,CrtGr2D;

VAR  A,B,C,x,y,z : REAL_AR;
     CrtPen : HPEN;

  Procedure Lorentz;

    Procedure f;
    Const delta = 0.01;
    Var   dx,dy,dz: REAL_AR;
    begin
      dx:=A*(y-x);
      dy:=x*(B-z)-y;
      dz:=x*y-C*z;
      x:=x+delta*dx;
      y:=y+delta*dy;
      z:=z+delta*dz
    end;

  Begin
    f;
    MoveXY(CrtDc,1+x,z);
    Repeat
      f;
      LineXY(CrtDc,1+x,z)
    Until KeyPressed;
    gotoxy(50,26); write('LORENTZ ATTRACTOR');
  End;

  {main program}
  BEGIN
    WinCrtInit('LORENTZ ATTRACTOR');
    CrtPen:=CreatePen(ps_Solid,1,RGB(0,0,255));
    SelectObject(CrtDC,CrtPen);
    A:=10; B:=30; C:=2.6666;
    x:=1; y:=1; z:=1;
    Fenetre(-22,25,-8,60);
    Cloture(15,MaxX-25,80,MaxY-15);
    Repeat
      clrscr;
      Bordure(CrtDc);
      Lorentz;
      SortieGraphique;
    Until rep='n';
    DoneWinCrt
  END.

{end of file Lorentz.pas}