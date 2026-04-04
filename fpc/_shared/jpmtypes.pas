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

  TIntArray    = array of integer;
  TFloatMatrix = array of TFloatArray;
  TMatrix      = array of TFloatArray;
  TFuncND      = function(var x: TFloatArray; n: integer): Float;

  { LM residual function: given params x (length n), fill fvec (length m) with residuals }
  TLMFunc = procedure(var x: TFloatArray; n: integer; var fvec: TFloatArray; m: integer);

{ Raise an exception with the given message. Call from self_test on failure. }
procedure SelfTestFail(const msg: string);
{ Assert cond; if false, raise with msg. }
procedure SelfTestCheck(cond: boolean; const msg: string);

const
  NaN = 1.0/0.0;

implementation

uses SysUtils;

procedure SelfTestFail(const msg: string);
begin
  raise Exception.Create('Self-test failure: ' + msg);
end;

procedure SelfTestCheck(cond: boolean; const msg: string);
begin
  if not cond then SelfTestFail(msg);
end;

end.

