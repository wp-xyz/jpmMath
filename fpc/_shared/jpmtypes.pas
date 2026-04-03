unit jpmTypes;

{$mode ObjFPC}{$H+}

interface

type
  float = double;

  TFloatArray = array of Float;

  TFunction1 = function(x: Float): float;
  TFunction2 = function(x, y: Float): float;

const
  NaN = 1.0/0.0;

implementation

end.

