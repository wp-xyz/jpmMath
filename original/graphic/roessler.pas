{****************************************************************************
*                   R O E S S L E R   A T T R A C T O R                     * 
* ------------------------------------------------------------------------- *
* TUTORIAL:                                                                 *
*                                                                           *
*  Until recently, the only known attractors were the fixed point, the      *
*  limit cycle and the torus. In 1963, Edwards Lorenz, meteorologist at the *
*  M.I.T.,discovered a practical example of a simple dynamic system presen- *
*  ting a complex behaviour. To adapt them to computers available at that   *
*  time, he began with simplifying the equations of metéorology to obtain   *
*  finally a model composed of three differential equations with three      *
*  unknown variables x, y, z and three parameters a, b, c:                  *
*                                                                           *
*         dx / dt = - a x + a y                                             *
*		  dy / dt = b x - y - x z                                           *
*		  dz / dt = -c z + x y                                              *
*                                                                           *
*  During very long simulations on a computer, Lorentz decided, to check a  *
*  result, to restart the same calculation halfway in order to spare time.  *
*  For that, he reinjected into the computer the intermediate date obtaine  *
*  earlier. He was very surprised to see that the new results were comple-  *
*  tely different from the first ones. After he suspected some computer     *
*  failure, Lorenz understood at last that the big difference between both  *
*  solutions came from very small differences in data. These small pertur-  *
*  bations exponentially amplified themselves, being doubled every four     *
*  days in simulated time, so after two months the results became entirely  *
*  different !                                                              *
*                                                                           *
*  Lorenz then realized that it would be very difficult to make meteorolo-  *
*  gical foresights in the long term, the slightest change in initial con-  *
*  ditions leading to a radically different evolution of atmosphere.        *
*                                                                           *
*  This is still the case today with atmospheric models much more sophis-   *
*  ticated.                                                                 *
*                                                                           *
*  None of the three attractors known at the time could predict the beha-   *
*  viour of such a dynamic system. Lorenz had just discovered a strange or  *
*  "chaotic" attractor to which his name was given, see program lorentz.pas *
*                                                                           *
*  The Roessler attractor is similar to the Lorenz attractor, by taking the *
*  following differential system of equations:                              *
*                                                                           *
*		  dx / dt = -( y + z )                                              *
*		  dy / dt = x + ( y / 5 )                                           *
*		  dz / dt = 1/5 + z ( x - 5,7 )                                     *
*                                                                           *
* ------------------------------------------------------------------------- *
*  REFERENCE:                                                               *
*     "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0        *
*      By R. Dony - MASSON, Paris 1990" [BIBLI 12].                         *
*                                                                           *
*                                          TPW Version By J-P Moreau        *
*                                              (www.jpmoreau.fr)            *
****************************************************************************}
Program Roessler_Attractor;
Uses WinCrtMy,WinProcs,Type_def,CrtGr2D;

VAR  A,B,C,x,y,z : REAL_AR;
     flag : boolean;
     i    : word;

  Procedure Roessler;

    Procedure f;
    Const delta = 0.005;
    Var   dx,dy,dz: REAL_AR;
    begin
      dx:=-(y+z);
      dy:=x+y*A;
      dz:=B+z*(x-C);
      x:=x+delta*dx;
      y:=y+delta*dy;
      z:=z+delta*dz
    end;

  Begin
    
    f; 
    if Not flag then
      Repeat
        f;
        Inc(i)
      Until i=1000;
    MoveXY(CrtDc,x,y+z+z);
    Repeat
      f;
      LineXY(CrtDc,x,y+z+z)
    Until KeyPressed
  End;


  {main program}
  BEGIN
    WinCrtInit('ROESSLER ATTRACTOR');
      A:=0.2; B:=0.2; C:=5.7;
      x:=-10; y:=-1; z:=-1;
      flag:=FALSE;  i:=0;
      Fenetre(-10,13,-17,55);
      Cloture(15,MaxX-25,80,MaxY-15);
      Bordure(CrtDc);
      Repeat
        Roessler;
        TextOut(CrtDC,25,20,'ROESSLER  ATTRACTOR',19);
        SortieGraphique;
        flag:=TRUE
      Until rep='n';
    DoneWinCrt
  END.

{end of file roessler.pas}