{**************************************************************************
*                   FRACTALS:  FEIGENBAUM  DIAGRAM                        *
* ----------------------------------------------------------------------- *
* TUTORIAL:                                                               *
*                                                                         *
* The dynamic process that allows to very simply simulate "chaos" has the *
* following mathematical form:                                            *
*		  x    = f (x , C ).                                      *
*                  n+1       n                                            *
* The only condition is that the relation between input x n and output    *
* x n+1 be non linear. If this iterative process starts with an arbitrary * 
* value x0, it will give a series of values:  x1, x2,...,x n... The long  *
* term behaviour of this series is interesting.                           *
*                                                                         *  
* Let us consider the classical example of the growth of a population over*
* several years. Let be an initial population x0 which after n years has  *
* become x n. The growth ratio is:                                        *
*                 r = ( x n+1 - x n ) / x n.                              *
* If this ratio remains constant, the dynamic law is:                     *
*                                                                         *
*                 x n+1 = ( 1 + r )  x n                                  *
*                                                                         *
* After n years, the population will be:    x n = ( 1 + r )  x0.          *
*                                                                         *
* To clarify the situation, let us consider the case where x0 = 0,001.    *
* If the growth ratio equals 5%, the population will roughly double       *
* every 15 years. As a matter of fact:                                    *
*                                                                         *
*   x0 = 0,001	       x1 = 0,00105	 x2 = 0,0011025                   *
*   x3 = 0,00115763    x4 = 0,00121551	 x5 = 0,00127628                  *
*   x6 = 0,00134010    x7 = 0,00140710	 x8 = 0,00147746                  *
*   x9 = 0,00155133    x10 = 0,00162889	 x11 = 0,00171034                 *
*   x12 = 0,00179586   x13 = 0,00188565	 x14 = 0,00197993                 *
*   x15 = 0,00207893   ...                                                *
*                                                                         *
* This kind of growth is exponential. But this dynamic law is linear,     *
* which is not judicious. Actually, the real growth of a poulation is     *
* necessarily limited by the ressources of its habitat, which are not in  *
* infinite quantity, such as food, for example. The belgian Verhulst was  *
* one of the first to study this phenomenon in 1845.                      *
* He advised to consider a variable growth ratio r, taking the form       *
* r = r - C x n.  The dynamic law of growth then takes the form:          * 
*                                                                         *
*		  x n+1 = ( 1 + r ) x n - C x² n                          * 
*                                                                         *
* By having C = r / X, the population increases up to the value X, then   *
* stabilizes at this value. At least, this remains true so long as the    *
* growth ratio is < 200 %. A human population has never such a high growth*  
* ratio, but in the case of insects, for example, this can be quite       *
* possible. For a growth ratio even higher, one can observe surprising    *
* results (see verhulst.pas program).                                     *
*                                                                         *
* The calculation begins at x0 = 0,1 X.                                   *
*                                                                         *
* Case r = 1.8                                                            *
*                                                                         *
* The response curve climbs up rapidly and after some oscillations reaches*
* a limit that is an "attractor".                                         *
*                                                                         *
* Case r = 2.3                                                            *
*                                                                         *
* The curve oscilates rapidly between two levels that frame the value X.  *
* The suite has two attractors.                                           *  
*                                                                         *
* Case r = 2.5                                                            *
*                                                                         *
* The suite has four attractors.                                          *
*                                                                         *
* Case r = 3.0                                                            *
*                                                                         *
* The numbers x jump from one value to another one, without having twice  *
* the same value. Such a behaviour can be qualified as "chaotic".         *
*                                                                         *
* The Feigenbaum diagram:                                                 *
*                                                                         *
* To better observe the behaviour of such suites x when r varies, we now  *
* only consider the "attractors" for each r value. The first 100 transient*
* values are skipped then at each r value 300 points are displayed.       *
* For r < 2,57, the behaviour is non chaotic: the attractors are in       *
* limited number. When r > 2,57, the attractors become queer and the dia- *
* gram has more and more ramifications until being fully inextricable:    *
* now we have a chaotic bahaviour!                                        *
*                                                                         *
*  The obtained picture is called the bifurcation diagram or Feigenbaum   *
* tree. By an accurate analysis of the bifurcation points, the mathema-   *
* tician Feigenbaum discovered a new universal constant. The length of the*
* r intervals for which a stable period is obtained,  is shortened, when  *
* the period is doubled, by a factor that tends toward the universal      *
* constant k = 4,669201660910...                                          *
*                                                                         *
* This constant, called the Feigenbaum Constant, can be found in other    *
* chaotic phenomenons, such as fluidic turbulences, chemical reactions,   *
* and even in human heart!                                                *
* ----------------------------------------------------------------------- *
* From "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0     *
* By R. Dony - MASSON, Paris 1990, pages 189-192" [BIBLI 12].             *
*                                                                         *
*                                      TPW version by J-P Moreau, Paris   *
*                                             (www.jpmoreau.fr)           *
**************************************************************************}
 Program CHAOS;   {WINCRT version with printing capability}
 Uses WinCrtMy,WinTypes,WinProcs,WObjects,Strings,Type_def,CrtGr2D,Winprint;

 Label 100;

 Var  dr,r,x : real_ar;
      i,n    : integer;
      ch     : char;
      OldPen,CrtPen : HPen;
      Pinfo  : PPrinterInfo;


   Procedure Draw_chaos(P:HDC; flag:boolean);
   var margex : integer;
   begin
     if Not flag then
     begin
       CrtPen:=CreatePen(ps_Solid,1,RGB(127,0,255));
       OldPen:=SelectObject(P,CrtPen);    { plume bleue }
     end;
     Fenetre(1.9,3.0,0.0,1.6);
     if MaxX > 1000 then margex:=200 else margex:=45;
     Cloture(margex,MaxX-25,95,MaxY-5);
     Axes(P);
     Grille(P,0.1,0.1);
     Gradue(P,0.2,0.3);
     Bordure(P);
     r:=1.9; dr:=0.0025;
     {write titles}
     if MaxX < 1000 then   {screen}
     begin
       TextOut(P,75,30,'FEIGENBAUM DIAGRAM',18);
       TextOut(P,50,10,'X',1);
       TextOut(P,490,295,'R',1)
     end
     else                  {HP laser printer}
     begin                 
       TextOut(P,300,80,'FEIGENBAUM DIAGRAM',18);
       TextOut(P,100,27,'X',1);
       TextOut(P,2950,2050,'R',1)
     end;
     {main r loop}
     repeat
       x:=0.3;
       for i:=1 to 200 do x:=(1+r)*x-r*x*x;
       for i:=1 to 300 do
       begin
         x:=(1+r)*x-r*x*x;
         if r>=1.95 then
         begin
           MoveXY(P,r,x);
	   LineXY(P,r+dr,x);
         end
       end;
       r:=r+dr
     until r > 3;
     if Not flag then
     begin
       SelectObject(P,OldPen);
       DeleteObject(CrtPen)
     end
   end;


 {main program}
 Begin
   WinCrtInit('CHAOS');
   Draw_chaos(CrtDc,FALSE);
   Sortiegraphique;
   if rep='i' then   {send to printer}
   begin
     New(Pinfo,Init);
     With Pinfo^ do
     BEGIN
       MaxX:=3200; MaxY:=2200; {HP laser landscape 300 dpi}
       StartDoc('chaos');
       Draw_chaos(PrintDc,TRUE);
       NewFrame;
       EndDoc
     END;
     Dispose(Pinfo,Done)
   end;
   DoneWinCrt  
 End.

{end of file chaos.pas}