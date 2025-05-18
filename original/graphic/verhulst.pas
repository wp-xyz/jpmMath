{**************************************************************************
*                   FRACTALS:  THE VERHULST DIAGRAM                       *
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
*   x0 = 0,001	       x1 = 0,00105	     x2 = 0,0011025                   *
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
*		  x n+1 = ( 1 + r ) x n - C x² n                                  * 
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
* See also program chaos.pas (Feigenbaum diagram).                        *
*                                                                         *
* ----------------------------------------------------------------------- *
* From "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0     *
* By R. Dony - MASSON, Paris 1990, page 189" [BIBLI 12].                  *
*                                                                         *
*                                      TPW Version By J-P Moreau, Paris   *
*                                             (www.jpmoreau.fr)           *
**************************************************************************}
 Program verhulst;  {WINCRT version with printing capability}
 Uses WinCrtMy,WinProcs,Strings,Type_def,Crtgr2d,WinPrint;

 var r,x,xx : REAL_AR;
     ch     : string[6];
     ch1    : string[8];
     s      : array[0..20] of char;
     PInfo  : PPrinterInfo;


 BEGIN
   WinCrtInit('VERHULST');
   Repeat
     Clrscr;
     writeln;
     write('  Input value of R ( 0 to 3 ): '); readln(r);
     if r<0.0 then r:=0.0;
     if r>3.0 then r:=3.0;
     ClrScr;
     Fenetre(0.0,170,0.0,1.5);
     Cloture(50,MaxX-25,95,MaxY-5);
     gotoxy(30,2); write('VERHULST DIAGRAM');
     Grille(CrtDc,10.0,0.1);
     Gradue(CrtDc,20.0,0.2);
     Bordure(CrtDc);
     x:=1.0; xx:=0.1; 
     MoveXY(CrtDc,0.0,xx);
     Repeat
       x:=x+1.0;
       xx:=(1+r)*xx-r*xx*xx;
       LineXY(CrtDc,x,xx);
     Until x > 165.0;
     str(r,ch); ch1:=concat('R=',ch);
     gotoxy(60,21); write(ch1);
     Sortiegraphique;
     if rep='i' then
     begin
       New(PInfo,Init);
       With Pinfo^ do
       BEGIN
         MaxX:=3200; MaxY:=2200;
         Cloture(200,MaxX,62,MaxY-25);
         StartDoc('verhulst');
         TextOut(PrintDc,1400,75,'VERHULST DIAGRAM',16);
         Grille(PrintDc,10.0,0.1);
         Gradue(PrintDc,20.0,0.2);
         Bordure(PrintDc);
         x:=1.0; xx:=0.1; 
         MoveXY(PrintDc,0.0,xx);
         Repeat
           x:=x+1.0;
           xx:=(1+r)*xx-r*xx*xx;
           LineXY(PrintDc,x,xx);
         Until x > 165.0;
         str(r,ch); ch1:=concat('R=',ch);
         StrPCopy(s,ch1);
         TextOut(PrintDc,2750,1900,s,strlen(s));
         NewFrame;
         EndDoc;
         MaxX:=596; MaxY:=401
       END;
       Dispose(Pinfo,Done);
       Sortiegraphique
     end
   Until rep='n';
   DoneWinCrt
 END.

{End of file verhulst.pas}