{****************************************************
*            RESOLUTION OF A TRIANGLE               *
* ------------------------------------------------- *
* Given three side or angle elements of a triangle  *
* out of six, this program will determine the mis-  *
* sing elements and calculate the surface. The not  *
* given elements must be put to zero.               *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
*  RESOLUTION OF A TRIANGLE                         *
*                                                   *
*    A = 18                                         *
*    B = 0                                          *
*    C = 0                                          *
*                                                   *
*    Angle unit: PI = 180                           *
*                                                   *
*    TA (angle opposite to side A) = 0              *
*    TB (angle opposite to side B) = 110            *
*    TC (angle opposite to side C) = 52.2           *
*                                                   *
*                                                   *
*  A  =       18.0000000000                         *
*  B  =       55.3311316842                         *
*  C  =       46.5260342304                         *
*                                                   *
*  TA =       17.8000000000                         *
*  TB =      110.0000000000                         *
*  TC =       52.2000000000                         *
*                                                   *
*  Surface:      393.4815393691                     *
*                                                   *
*                       TPW Version By J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
****************************************************}
Program Triangle;
Uses WinCrt;

VAR
     a,b,c,ta,tb,tc,one_pi: REAL;
     two: BOOLEAN; 
   

  Function Triangle1(var a,b,c,ta,tb,tc,s: REAL; var two:boolean): boolean;
  Var
       side : array[1..3] of REAL;
       angle: array[1..3] of REAL;
       index: array[1..3] of BYTE;
       nbside,nbangle: INTEGER;
       i             : INTEGER;
       r             : REAL;

       Procedure Exchange(i,j:INTEGER);
       var r:REAL; n:BYTE;
       begin
         n:=index[i]; index[i]:=index[j]; index[j]:=n;
         r:=angle[i]; angle[i]:=angle[j]; angle[j]:=r;
         r:=side [i]; side [i]:=side [j]; side [j]:=r
       end;

       Procedure CalculateAngle(i:INTEGER);
       var j,k:INTEGER; r:REAL;
       begin
         j:=1+i MOD 3;  k:=1+j MOD 3;
         r:=(SQR(side[j])+SQR(side[k])-SQR(side[i]))/2/side[j]/side[k];
         angle[i]:=PI/2-ArcTan(r/SQRT(1-r*r))
       end;

  Begin
    {------------ prepare index tables -------------}
    side[1]:=a; angle[1]:=ta;
    side[2]:=b; angle[2]:=tb;
    side[3]:=c; angle[3]:=tc;
    for i:=1 to 3 do index[i]:=i;
    {------------ verify consistency of data -------}
    Triangle1:=FALSE; nbside:=0; nbangle:=0;
    {angles must be in 0,PI and sides > 0}
    for i:=1 to 3 do
      if (angle[i]<0) or (angle[i]>PI) or (side[i]<0) then exit;
    for i:=1 to 3 do
      if side[i]>0 then Inc(nbside);
    for i:=1 to 3 do
      if angle[i]>0 then Inc(nbangle);
    if (nbside=0) or (nbside+nbangle<>3) or (angle[1]+angle[2]+angle[3]>PI) then exit;
    {end verify consistency of data}
    Case nbside of
      {-------------------- 1 side given -----------}
      1: begin
           for i:=1 to 3 do
             if angle[i]=0 then  {calculate missing angle}
               angle[i]:=PI-angle[1]-angle[2]-angle[3];
           if side[2]>0 then Exchange(1,2);  {put given side}
           if side[3]>0 then Exchange(1,3);  {in position 1 }
           side[2]:=side[1]*sin(angle[2])/sin(angle[1]);
           side[3]:=side[1]*sin(angle[3])/sin(angle[1])
         end;
      {-------------------- 2 sides given ----------}
      2: begin
           if side[3]>0 then  {put unknown side in position 3}
             if side[1]>0 then Exchange(2,3)
                          else Exchange(1,3);
             if angle[3]>0 then  {a, b, tc given} 
             begin
               side[3]:=SQRT(a*a+b*b-2*a*b*cos(angle[3]));
               for i:=1 to 2 do CalculateAngle(i)
             end
             else                {a, b, ta given}
             begin
               if angle[1]=0 then Exchange(1,2);
               r:=side[2]*sin(angle[1]);
               if r>side[1] then exit;
               r:=SQRT(SQR(side[1])-SQR(r));
               if (angle[1] >= PI/2) and (side[1] <= side[2]) then exit;
               if (angle[i] < PI/2) and (side[1] < side[2]) then
                 if two then
                   side[3]:=side[2]*cos(angle[1]) - r
                 else
                 begin
                   side[3]:=side[2]*cos(angle[1]) + r;
                   two:=TRUE
                 end
               else
                 side[3]:=side[2]*cos(angle[1]) + r;
               for i:=2 to 3 do CalculateAngle(i)
             end
           end;
      {-------------------- 3 sides given ----------}
      3: begin
           if (c <= ABS(a-b)) or (c >= a+b) then exit;
           for i:=1 to 3 do  CalculateAngle(i)
         end
    End; {case       
    ----------------------- desindex ---------------}
    a:=side[index[1]]; ta:=angle[index[1]];
    b:=side[index[2]]; tb:=angle[index[2]];
    c:=side[index[3]]; tc:=angle[index[3]];
    {---------------------- calculate surface ------}
    r:=(a+b+c)/2;
    s:=SQRT(r*(r-a)*(r-b)*(r-c));
    Triangle1:=TRUE  {success}

  End; {Function Triangle1}


  Procedure DisplayTriangle;
  Var xa,xb,xc,xta,xtb,xtc,s: REAL;
  Begin
    xa:=a; xb:=b; xc:=c;
    xta:=ta; xtb:=tb; xtc:=tc;
    if Triangle1(xa,xb,xc,xta,xtb,xtc,s,two) then
    begin
      writeln('  A  = ',xa:21:10);
      writeln('  B  = ',xb:21:10);
      writeln('  C  = ',xc:21:10);
      writeln;
      writeln('  TA = ',xta*(one_pi/PI):21:10);
      writeln('  TB = ',xtb*(one_pi/PI):21:10);
      writeln('  TC = ',xtc*(one_pi/PI):21:10);
      writeln;
      writeln('  Surface: ', s:21:10)
    end
    else
      writeln(' Wrong data or no solution found.');
      writeln
  End;


{main program}
BEGIN
  writeln;
  writeln('                 RESOLUTION OF A TRIANGLE');
  writeln;
  writeln(' Give three side or angle elements out of six. The other elements');
  writeln(' must be given as zero.  The program  will determine  the missing');
  writeln(' elements and display results.');
  writeln;
  write('      Side A = '); readln(a);
  write('      Side B = '); readln(b);
  write('      Side C = '); readln(c);
  writeln;
  write('      Angle unit: PI = '); readln(one_pi);
  writeln;
  write('      TA (angle opposite to side A) = '); readln(ta);
  write('      TB (angle opposite to side B) = '); readln(tb);
  write('      TC (angle opposite to side C) = '); readln(tc);
  ta:=ta*PI/one_pi; tb:=tb*PI/one_pi; tc:=tc*PI/one_pi;
  writeln;
  two:=FALSE;
  DisplayTriangle;
  if two then
  begin
    writeln(' There is a second triangle solution:');
    Readkey;
    writeln;
    DisplayTriangle
  end;
  Readkey;
  DoneWinCrt
END.

{end of file triangle.pas}