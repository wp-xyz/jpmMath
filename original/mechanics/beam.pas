{*****************************************************************************
* Calculates beam deflection (in) for four different support/loading systems *
* -------------------------------------------------------------------------- *
* SAMPLE RUN:                                                                *
*                                                                            *
* Give elasticity (lb/in^2) and moment of inertia (in^4): 30e6 797           *
* Give the beam length in ft: 20                                             *
* Choose one of these support/loading systems:                               *
* 1 - supported at each end, concentrated load                               *
* 2 - supported at each end, uniformly distributed load                      *
* 3 - supported at one end, concentrated load at free end                    *
* 4 - supported at one end, distributed load                                 *
* Input your choice (1 to 4): 1                                              *
* Give the concentrated force: 50000                                         *
*                                                                            *
* Deflection =  -0.6023                                                      *
*                                                                            *
* More? (y/n): y                                                             *
*                                                                            *
* Give elasticity (lb/in^2) and moment of inertia (in^4): 30e6 797           *
* Give the beam length in ft: 20                                             *
* Choose one of these support/loading systems:                               *
* 1 - supported at each end, concentrated load                               *
* 2 - supported at each end, uniformly distributed load                      *
* 3 - supported at one end, concentrated load at free end                    *
* 4 - supported at one end, distributed load                                 *
* Input your choice (1 to 4): 3                                              *
* Give the concentrated force: 10000                                         *
*                                                                            *
* Deflection =  -1.9272                                                      *
*                                                                            *
* More? (y/n): n                                                             *
* -------------------------------------------------------------------------- *
* Reference: "Problem Solving with Fortran 90 By David R. Brooks, Springer-  *
*             Verlag New York, 1997".                                        *
*                                                                            *
*                                       Pascal Version By J-P Moreau, Paris. *
*                                               (www.jpmoreau.fr)            *
*****************************************************************************}
    PROGRAM Beam2;
    USES WINCRT;

    VAR
          elasticity,        {lb/in^2}
          moment_of_inertia, {in^4}
          length,            {ft}
          load,              {lb}
          deflection: REAL;  {in}

          systemID: INTEGER; {1 - supported at each end, concentrated load
                              2 - supported at each end, distributed load
                              3 - supported one end, concentrated at free end
                              4 - supported one end, distributed }
          YesNo: CHAR;


    BEGIN

      Repeat
        Clrscr;
        Writeln;
        Write(' Give elasticity (lb/in^2) and moment of inertia (in^4): ');
        Readln(elasticity, moment_of_inertia);
        Write(' Give the beam length in ft: ');
        Readln(length);
        Writeln(' Choose one of these support/loading systems:');
        Writeln(' 1 - supported at each end, concentrated load');
        Writeln(' 2 - supported at each end, uniformly distributed load');
        Writeln(' 3 - supported at one end, concentrated load at free end');
        Writeln(' 4 - supported at one end, distributed load');
        Write(' Input your choice (1 to 4): '); Readln(systemID);
        CASE systemID of
          1: write(' Give the concentrated force: ');
          2: write(' Give the distributed weight: ');
          3: write(' Give the concentrated force: ');
          4: write(' Give the distributed weight: ');
        END;
        Readln(load);

        length:=length*12.0;
        CASE systemID of
          1: deflection:=-load*length*length*length/(48.0*elasticity*moment_of_inertia);
          2: deflection:=-5.0*load*length*length*length/(384.0*elasticity*moment_of_inertia);
          3: deflection:=-load*length*length*length/(3.0*elasticity*moment_of_inertia);
          4: deflection:=-load*length*length*length/(8.0*elasticity*moment_of_inertia)
        END;

        Writeln;
        Writeln(' Deflection = ', deflection:8:4);
        Writeln;
        Write(' More? (y/n): '); YesNo:=Readkey;
      Until YesNo<>'y';

      DoneWinCrt

    END.

{end of file beam.pas}