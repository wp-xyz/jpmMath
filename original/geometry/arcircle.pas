{****************************************************
* This program allows calculating the 4 parameters  *
* of an arc of circle knowing two of them.          *
*  R:radius of circle S:curvilinear abscissa of arc *
*  C:length of chord  T:angle of arc of circle (rad)*
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
*   R=0                                             *
*   S=0                                             *
*   C=2                                             *
*   T=1                                             *
*                                                   *
* SUCCESS: TRUE                                     *
*                                                   *
* Radius              :   2.0858296429              *
* Curvilinear abscissa:   2.0858296429              *
* Chord               :   2.0000000000              *
* Angle (radians)     :   1.0000000000              *
*                                                   *
* Surface of angular sector     :   2.1753426497    *
* Surface between arc and chord :   0.3448549280    *
*                                                   *
****************************************************}
Program arc_of_circle;
Uses WinCrt;

LABEL fin;

VAR
     c,r,s,t,s1,s2 : REAL;
     a,b,u,x,y: REAL;
     ndata: INTEGER;
     done : BOOLEAN;

BEGIN
  writeln;
  writeln('             ARC OF CIRCLE');
  writeln(' -----------------------------------------');
  writeln('  Give two parameters between:');
  writeln('    R: radius of circle');
  writeln('    S: curvilinear abscissa');
  writeln('    C: chord');
  writeln('    T: angle of arc in radians');
  writeln('  The other two must be given as zero.');
  writeln('  The program then calculates the');
  writeln('  missing parameters and the two surfaces:');
  writeln('    S1: total surface of angular sector');
  writeln('    S2: surface between chord and arc s.');
  writeln(' -----------------------------------------');
  write('     R= '); read(r);
  write('     S= '); read(s);
  write('     C= '); read(c);
  write('     T= '); read(t);
  writeln;
  done:=FALSE; ndata:=0;
  {verify consistency of data r,s,c,t}
  if r<>0 then Inc(ndata);
  if s<>0 then Inc(ndata);
  if c<>0 then Inc(ndata);
  if t<>0 then Inc(ndata);
  if ndata<>2 then
  begin
    writeln(' Too many data <> 0 (maximum two).');
    goto fin
  end;
  if (t<0) or (t>PI) then
  begin
    writeln(' Angle t(radians) must be between 0 and PI.');
    goto fin
  end;
  if (s>0) and (c>=s) then
  begin
    writeln(' Chord c must be smaller than curvilinear abscissa s.');
    goto fin
  end;
  if (c>0) and (s>c*PI/2) then
  begin
    writeln(' Curvilinear abscissa s must be smaller than chord times PI/2.');
    goto fin
  end;
  if (r>0) and (s>PI*r) then
  begin
    writeln(' Curvilinear abscissa s must be smaller then radius r times PI.');
    goto fin
  end;
  {end of verify data}
  if t>0 then
    if r>0 then
      {---------------------- t and r are given ----}
      begin
        s:=r*t; c:=2*r*sin(t/2)
      end
    else
      if s>0 then
      {---------------------- t and s are given ----}
      begin
        r:=s/t; c:=2*r*sin(t/2)
      end
      else
      {---------------------- t and c are given ----}
      begin
        r:=c/2/sin(t/2); s:=r*t
      end
  else {t is not given}
    if c=0 then
      {---------------------- r and s are given ----}
      begin
        t:=s/r; c:=2*r*sin(t/2)
      end
    else
      if r>0 then
      {---------------------- r and c are given ----}
      begin
        if c<2*r then
          t:=2*arctan(c/sqrt(4*r*r-c*c))
        else
          t:=PI;
        s:=r*t
      end
      else
      {---------------------- s and c are given}
      begin
        y:=c;
        repeat
          x:=y;        u:=s/x;
          a:=s*cos(u); b:=x*sin(u)-a;
          if (c-a)*b<=0 then goto fin;
          y:=x*(c-a)/b
        until abs(x-y)<1E-8;
        r:=y/2; t:=s/r
      end;
  {calculate surfaces s1 and s2}
  s1:=s*r/2;
  s2:=(s-c*cos(t/2))*r/2;
  done:=TRUE;

  {print results}
  Clrscr;
  writeln;
  writeln(' SUCCESS: ',done);
  writeln;
  if done then
  begin
    writeln(' Radius              : ',r:15:10);
    writeln(' Curvilinear abscissa: ',s:15:10);
    writeln(' Chord               : ',c:15:10);
    writeln(' Angle (radians)     : ',t:15:10);
    writeln;
    writeln(' Surface of angular sector     : ',s1:15:10);
    writeln(' Surface between arc and chord : ',s2:15:10);
    writeln
  end;

fin: readkey; donewincrt
END.

{end of file arcircle.pas}