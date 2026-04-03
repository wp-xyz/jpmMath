unit jpmTypes;

{$mode ObjFPC}{$H+}

interface

type
  float = double;

  TFloatArray = array of float;
  TFunction1 = function(x: float): float;
  TFunction2 = function(x, y: float): float;
  TFunction3 = function(x, y, z: float): float;
  TFunctionN = function(var x: TFloatArray): Float;

const
  NaN = 1.0/0.0;

implementation

end.

