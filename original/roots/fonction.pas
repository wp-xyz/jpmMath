{--------------------------------------------------------------------*
*      Definition of six test functions for pegasus algorithm        *
*      [BIBLI 11]                                                    *
*--------------------------------------------------------------------}
UNIT Fonction;

INTERFACE

  Var NumFunc: Integer;

  Function fct(x:Double): Double;


IMPLEMENTATION

  Function fct(x:Double): Double;
  Begin
    Case NumFunc of
      1: fct := 5*x -exp(x);
      2: fct := (((((x-6)*x+15)*x-20)*x+15)*x-6)*x+1;
      3: fct := sin(x);
      4: fct := 1.0 + sin(x);
      5: fct := exp(x) -(1.0 + x + x*x/2);
      6: fct := sqr(x-1)*(SIN(PI*x)-ln(2*x/(x+1)))
    end
  End;
  
END.
