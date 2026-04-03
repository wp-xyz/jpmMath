unit jpmTypes;

{$mode ObjFPC}{$H+}

interface

type
  float = double;

  TFloatArray = array of Float;

  TFunction1 = function(x: Float): float;
  TFunction2 = function(x, y: Float): float;
  TFunctionN = function(var x: TFloatArray): Float;

const
  NaN = 1.0/0.0;

implementation

end.

