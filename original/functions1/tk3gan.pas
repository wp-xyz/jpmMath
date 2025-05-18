{************************************************************************
* Test programm for cubature over triangles via summed Gaussian n-point *
* formula                                                               *
* --------------------------------------------------------------------- *
* SAMPLE RUN:                                                           *
* (Integrate function EXP(-(x*x + y*y)) over triangle defined by three  *
*  points (0,0), (10,0) and (0,10). )                                   *
*                                                                       *
* # Method     Number       Value      Error Code                       *
*           Subtriangles   Integral    (must be 0)                      *
* -------------------------------------------------                     *
*     1          1       0.0000000112      0                            *
*     1          4       0.0483240045      0                            *
*     1          9       0.4706075316      0                            *
*     1         16       0.7913533084      0                            *
*     1         25       0.8814598010      0                            *
*     2          1       0.0000000000      0                            *
*     2          4       0.0120496973      0                            *
*     2          9       0.2538493244      0                            *
*     2         16       0.5593099250      0                            *
*     2         25       0.7075499085      0                            * 
*     3          1       0.0644320023      0                            *
*     3          4       1.0390297430      0                            *
*     3          9       1.0188794444      0                            *
*     3         16       0.8590974900      0                            *
*     3         25       0.7974378946      0                            *
*     7          1       0.8091876107      0                            *
*     7          4       0.9654995491      0                            *
*     7          9       0.7922397459      0                            *
*     7         16       0.7849640612      0                            *
*     7         25       0.7862848602      0                            *
* -------------------------------------------------                     *
*                                                                       *
* --------------------------------------------------------------------- *
* Reference:                                                            *
* "Numerical Algorithms with C, By Gisela Engeln-Muellges               *
*  and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].                  *
*                                                                       *
*                                 TPW Release 1.0 By J-P Moreau, Paris. *
*                                          (www.jpmoreau.fr)            *
************************************************************************} 
Program TK3GAN;

Uses WinCrt, Kubgauss;

Var
    i, Verfahren, Kantenteilung: Integer;
    WertDerKubatur: Double;
    Px, Py, Qx, Qy, Rx, Ry: Double;

{Note: Function func(x,y) to integrate is defined in unit Kubgauss.}


BEGIN

  {define summits of triangle}
  Px :=  0.0;  Py :=  0.0;
  Qx := 10.0;  Qy :=  0.0;
  Rx :=  0.0;  Ry := 10.0;

  {print header}
  writeln(' # Method     Number       Value      Error Code  ');
  writeln('           Subtriangles   Integral    (must be 0) ');
  writeln(' -------------------------------------------------');
  {main integration loop}
  for Verfahren := 1 to 7 do
    {only methods 1,2,3,7 are implemented}
    if (Verfahren < 4) or (Verfahren > 6) then
    begin
      for Kantenteilung := 1 to 5 do
      begin
        write('     ', Verfahren);   
        write('         ', (Kantenteilung * Kantenteilung):2);

        {call appropriate Gaussian cubature method}
        i := Kub3GauN (Px, Py, Qx, Qy, Rx, Ry, Verfahren,
                       Kantenteilung,
                       WertDerKubatur);

        writeln(' ',WertDerKubatur:18:10,'      ',i)

      end
    end;
    writeln(' -------------------------------------------------');

    ReadKey;
    DoneWinCrt

END.

{ end of file tk3gan.pas }