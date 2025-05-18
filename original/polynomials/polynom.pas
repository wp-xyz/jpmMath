{****************************************************
*   Elementary operations on Polynomials in Pascal  *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
*                                                   *
*                 TPW Version By J-P Moreau, Paris. *
*                        (www.jpmoreau.fr)          *
****************************************************}
UNIT POLYNOMS;

INTERFACE

  Uses WinCrt, WinProcs;

  Const     MAXPOL = 100;        {max degree of a polynomial}
            BIG    = 1E9;        {big number}
            SMALL  = 1E-8;

  Type
            NUMBER = Record
              is_real : BOOLEAN; {TRUE=real number,FALSE=fractional or integer number} 
              value   : REAL;    {value of real number}
              p,q     : LONGINT  {p/q if fractional number, q=1 for integer number}
            End;

            POLYNOM = Record
              degree: INTEGER;                     {degree of polynomial}
              coeff : Array[0..MAXPOL] of NUMBER   {coefficients of polynomial}
            End;


  Function  StrReal(r:REAL;long,deci:INTEGER):STRING;
  Function  GCD(a,b:REAL):REAL;
  Function  SetNumber(VAR n:NUMBER; ch:STRING): Boolean;
  Function  AddNumber(x,y:NUMBER; VAR z:NUMBER): Boolean;
  Function  MultNumber(x,y:NUMBER; VAR z:NUMBER): Boolean;
  Function  DivNumber(x,y:NUMBER; VAR z:NUMBER): Boolean;
  Procedure ReadNumber(tx:STRING; VAR n:NUMBER);
  Procedure WriteNumber(n:NUMBER);
  Function  EnterPolynom(tx:STRING; VAR P:POLYNOM): Boolean;
  Procedure DisplayPolynom(VAR P:POLYNOM);
  Function  EvaluatePolynom(VAR P:POLYNOM; x:NUMBER; VAR y:NUMBER):Boolean;


IMPLEMENTATION

  {return greater common divisor of two reals}
  Function GCD(a,b:REAL):REAL;
  Var  x,r: REAL;
  Begin
    a:=INT(ABS(a));
    b:=INT(ABS(b));
    gcd:=1;
    if (a>1E10) or (b>1E10) then exit;
    if (a=0) or (b=0) then exit;
    if a<b then
    begin
      x:=a; a:=b; b:=x
    end;
    Repeat
      r:=a-b*INT(a/b);
      a:=b; b:=r
    Until abs(r)<1E-10;
    gcd:=a
  End;

  {put a string into a number. Valid strings: 5.56874 or 15/187 or 255}
  Function SetNumber(VAR n:NUMBER; ch:STRING): Boolean;
  Var  i,j,k: INTEGER; r: REAL;
  Begin
    SetNumber:=FALSE; n.is_real:=(Pos('.',ch)>0);  {detect a point in string}
    if n.is_real then
    begin
      Val(ch,n.value,j); if j<>0 then exit  {error in number}
    end
    else
    begin
      j:=Pos('/',ch);  {detect / in string}
      if j>0 then
      begin            {case of a rational number}
        Val(Copy(ch,1,Pred(j)),n.p,k);
        if k<>0 then exit;  {error in numerator}
        Val(Copy(ch,j+1,length(ch)-j),n.q,k);
        if k<>0 then exit;  {error in denominator}
        {simplify p and q ?}
        r:=GCD(n.p,n.q);
        n.p:=n.p div Round(r); n.q:=n.q div Round(r); n.value:=n.p/n.q
      end
      else
      begin            {case of an integer}
        Val(ch,n.value,j); n.p:=Round(n.value); n.q:=1;
        if j<>0 then exit {error in integer}          
      end
    end;
    SetNumber:=TRUE
  End;

  {Add two numbers of type NUMBER (real, integer or fractional) }
  Function AddNumber(x,y:NUMBER; VAR z:NUMBER): Boolean;
  Var  s:REAL; negative:Boolean;
  Begin
    AddNumber:=TRUE;
    z.is_real:=(x.is_real) or (y.is_real);
    z.value:=x.value+y.value;
    if x.value=0 then z:=y
    else
      if y.value=0 then z:=x
      else
      begin
        if (abs(x.p) > BIG) or (abs(x.q) > BIG) or
          (abs(y.p) > BIG) or (abs(y.q) > BIG) then
            z.is_real := TRUE;
        if Not z.is_real then
        begin
          s:=GCD(x.q,y.q);
          z.p:=x.p*y.q+x.q*y.p; z.q:=x.q*y.q;
          negative:=(z.p<0) XOR (z.q<0);
          z.p:=Abs(z.p) div Round(s);
          z.q:=Abs(z.q) div Round(s);
          if negative then z.p:=-z.p
        end
      end
  End;

  {Multiply two numbers of type NUMBER (real, integer or fractional) }
  Function MultNumber(x,y:NUMBER; VAR z:NUMBER): Boolean;
  Var  s:REAL; negative:Boolean;
  Begin
    MultNumber:=TRUE;
    z.is_real:=(x.is_real) or (y.is_real);
    z.value:=x.value * y.value;
    if z.value=0 then
    begin
      z.p:=0; z.q:=1
    end
    else
    begin
      if (abs(x.p) > BIG) or (abs(x.q) > BIG) or
        (abs(y.p) > BIG) or (abs(y.q) > BIG) then
          z.is_real := TRUE;
      if Not z.is_real then
      begin
        z.p:=x.p*y.p; z.q:=x.q*y.q;
        s:=GCD(z.p,z.q);
        negative:=(z.p<0) XOR (z.q<0);
        z.p:=Abs(z.p) div Round(s);
        z.q:=Abs(z.q) div Round(s);
        if negative then z.p:=-z.p
      end
    end
  End;

  {Divide two numbers of type NUMBER (real, integer or fractional) }
  Function  DivNumber(x,y:NUMBER; VAR z:NUMBER): Boolean;
  Var s:REAL; negative:Boolean;
  Begin
    DivNumber:=FALSE; if y.value=0 then exit;
    z.is_real:=(x.is_real) or (y.is_real);
    if (abs(x.p)>BIG) or (abs(x.q)>BIG) or
      (abs(y.p)>BIG) or (abs(y.q)>BIG) then z.is_real:=TRUE;
    if Not z.is_real then
    begin
      z.p:=x.p*y.q; z.q:=x.q*y.p; s:=GCD(z.p,z.q);
      negative:=(z.p<0) XOR (z.q<0);
      z.p:=abs(z.p) div Round(s); z.q:=abs(z.q) div Round(s);
      if negative then z.p:=-z.p
    end;
    z.value:=x.value/y.value;
    DivNumber:=TRUE
  End;

  {return a real number as a string}
  Function StrReal(r:REAL;long,deci:INTEGER):STRING;
  var prov:STRING;
  Begin
    if long>0 then Str(r:long:deci,prov) else Str(r,prov);
    StrReal:=prov
  End;

  {read from screen a number of type NUMBER (real, integer or fractional) }
  Procedure ReadNumber(tx:STRING; VAR n:NUMBER);
  Var  ch:STRING; i:INTEGER;
  Begin
    for i:=1 to 3 do
    begin
      write(tx); readln(ch); if SetNumber(n,ch) then exit;
      MessageBeep(0)
    end;
    Fillchar(n,sizeof(n),0); writeln(tx,0)
  End;

  {write to screen a number of type NUMBER (real, integer or fractional) }
  Procedure WriteNumber(n:NUMBER);
  Begin
    if n.value>0 then write('+ ') else write('- ');
    if n.is_real then
      write(StrReal(abs(n.value),8,4))
    else
    begin
      write(abs(n.p));
      if abs(n.q)<>1 then write('/',abs(n.q))
    end
  End;

  {convert a valid string into a polynomial of type POLYNOM.
   Example of valid string:   X5 +3/5X4 -12X2 +8X -1/4
   Note that a blank space is required between monomes.    }    
  Function EnterPolynom(tx:STRING; VAR P:POLYNOM): Boolean;

  Var  i: INTEGER;
       ch: STRING;

    {internal function}
    Function ExtractMonom(VAR ch:STRING; VAR P:POLYNOM): Boolean;
    Var  i,degree,xpos,sppos,signe: INTEGER;
         coef: NUMBER; chn: STRING;
    Begin
      ExtractMonom:=FALSE;
      if ch='' then exit;
      While ch[1]=' ' do Delete(ch,1,1);  {take away blanks}
      if ch[1]='-' then signe:=-1 else signe:=1;
      if (ch[1] in ['+','-']) and (length(ch)>1) then
        ch:=Copy(ch,2,pred(length(ch)));
      While ch[1]=' ' do Delete(ch,1,1);  {take away blanks}
      xpos:=Pos('X',ch);
      if xpos=0 then
      begin
        degree:=0;
        if Not SetNumber(coef,ch) then exit;
        ch:=''
      end
      else
      begin
        if xpos=1 then chn:='1'
                  else chn:=Copy(ch,1,pred(xpos));
        if Not SetNumber(coef,chn) then exit;
        sppos:=Pos(' ',ch);
        if sppos=0 then sppos:=Succ(length(ch));
        if xpos=Pred(sppos) then
          degree:=1
        else
        begin
          Val(Copy(ch,Succ(xpos),sppos-xpos-1),degree,i);
          if (i<>0) or (degree>MAXPOL) then exit
        end;
        if sppos <= length(ch) then ch:=Copy(ch,sppos+1,length(ch)-sppos)
                               else ch:=''
      end;
      if signe=-1 then
      begin
        coef.p:=-coef.p;
        coef.value:=-coef.value
      end;
      if Not AddNumber(coef,P.coeff[degree],P.coeff[degree]) then exit;
      if degree>P.degree then P.degree:=degree;
      ExtractMonom:=TRUE
    End;

    Begin  {of EnterPolynom}
      EnterPolynom:=TRUE;
      for i:=1 to 3 do
      begin
        write(tx); fillchar(P,sizeof(P),0);
        readln(ch);
        for i:=1 to length(ch) do ch[i]:=Upcase(ch[i]);
        While (length(ch)>0) and ExtractMonom(ch,P) do;
        if length(ch)=0 then exit;
        MessageBeep(0)
      end;
      writeln; EnterPolynom:=FALSE
    End;

    {This procedure displays to screen the symbolic
     representation of a polynom of type POLYNOM  }
    Procedure DisplayPolynom(VAR P:POLYNOM);
    Var  i,vx: INTEGER;
    Begin
      vx:=WhereX;
      writeln; writeln; GotoXY(vx,WhereY);
      if (P.degree=0) and (P.coeff[0].value=0) then write(' 0')
      else
        write(' ');
        for i:=P.degree downto 0 do
          if P.coeff[i].value<>0 then
          begin
            if WhereX>60 then
            begin
              writeln; writeln;
              GotoXY(Succ(vx),WhereY)
            end;
            if P.coeff[i].value < 0 then write('- ');
            {do not write initial '+ '}
            if (P.coeff[i].value > 0) and (i<P.degree) then write('+ ');
            if P.coeff[i].is_real then
              Write(StrReal(abs(P.coeff[i].value),6,2))
            else
            begin
              if (abs(P.coeff[i].p)<>1) or (abs(P.coeff[i].q)<>1) or (i=0) then
                Write(abs(P.coeff[i].p));
              if abs(P.coeff[i].q)<>1 then
                Write('/',abs(P.coeff[i].q))
            end;
            if i>1 then
            begin
              write(' X');
              GotoXY(WhereX, Pred(WhereY)); write(i);
              GotoXY(WhereX, Succ(WhereY)); write(' ')
            end
            else
              if i=1 then write(' X ')
          end;
      GotoXY(WhereX,Succ(WhereY))
    End;


    {return the value y of a polynomial of type POLYNOM for argument=x
     x and y are of type NUMBER  (real, integer or fractional)        }
    Function  EvaluatePolynom(VAR P:POLYNOM; x:NUMBER; VAR y:NUMBER):Boolean;
    Var  i: INTEGER;
    Begin
      EvaluatePolynom:=FALSE;
      if P.degree > MAXPOL then exit;
      if Not SetNumber(y,'0') then exit;
      for i:=0 to P.degree do
      begin
        if Not MultNumber(y,x,y) then exit;
        if Not AddNumber(y,P.coeff[P.degree-i],y) then exit
      end;
      EvaluatePolynom:=TRUE
    End;

END.

{end of file polynom.pas}