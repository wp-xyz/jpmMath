unit main;

{$mode objfpc}{$H+}

interface

uses
  SysUtils, Classes, Forms, Graphics, Controls, ExtCtrls, Dialogs,
  TAGraph, TASeries,
  jpmTypes, jpmIntegration;

type
  TMainForm = class(TForm)
    Chart: TChart;
    RadioGroup1: TRadioGroup;
    RombergSeries: TLineSeries;
    DataSeries: TLineSeries;
    Panel1: TPanel;
    procedure btnNewDataClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure RadioGroup1Click(Sender: TObject);
  private
    X, Y: Array of Double;
    procedure NewData;
    procedure PlotData;

  public

  end;

var
  MainForm: TMainForm;

implementation

{$R *.lfm}

const
  MIN =  0.0;
  MAX = 10.0;

function SinFunc(x: Float): Float;
begin
  Result := sin(x);
end;

function CosFunc(x: Float): Float;
begin
  Result := cos(x);
end;

function ExpFunc(x: Float): Float;
begin
  Result := exp(-x);
end;

function SincFunc(x: Float): Float;
begin
  if x <> 0.0 then
    Result := sin(x) / x
  else
    Result := 1.0;
end;

procedure TMainForm.btnNewDataClick(Sender: TObject);
begin
  NewData;
  PlotData;
end;

procedure TMainForm.FormCreate(Sender: TObject);
begin
  NewData;
  PlotData;
end;

procedure TMainForm.RadioGroup1Click(Sender: TObject);
begin
  NewData;
  PlotData;
end;

procedure TMainForm.NewData;
var
  i, n: Integer;
begin
  n := 100;
  SetLength(X, n);
  SetLength(Y, n);
  for i := 0 to n - 1 do
  begin
    X[i] := MIN + (MAX - MIN) * i / (n-1);
    case RadioGroup1.ItemIndex of
      0: Y[i] := SinFunc(X[i]);
      1: Y[i] := CosFunc(X[i]);
      2: Y[i] := ExpFunc(X[i]);
      3: Y[i] := SincFunc(X[i]);
      else raise Exception.Create('Function not implemented.');
    end;
  end;
end;

procedure TMainForm.PlotData;
var
  i, n, nIt: Integer;
  prec: Float;
  xx, yy: Double;
  f: TFunction1;
begin
  DataSeries.Clear;
  for i := 0 to High(X) do
    DataSeries.AddXY(X[i], Y[i], IntToStr(i));

  case RadioGroup1.ItemIndex of
    0: f := @SinFunc;
    1: f := @CosFunc;
    2: f := @ExpFunc;
    3: f := @SincFunc;
  end;

  n := 100;  // Number of points in integrated series
  RombergSeries.Clear;
  for i := 0 to n-1 do
  begin
    xx := MIN + (MAX - MIN) * i / (n-1);
    yy := RombergIntegral(f, MIN, xx, 1E-9, prec{%H-}, nIt{%H-}, 0, 50);
    RombergSeries.AddXY(xx, yy);
  end;
end;

end.

