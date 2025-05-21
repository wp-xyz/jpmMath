unit main;

{$mode objfpc}{$H+}

interface

uses
  Classes, ExtCtrls, StdCtrls, SysUtils, Forms, Controls, Graphics, Dialogs,
  TAGraph, TASeries, TASources;

type
  TMainForm = class(TForm)
    Button1: TButton;
    Chart1: TChart;
    AkimaSeries: TLineSeries;
    CheckBox1: TCheckBox;
    DataSeries: TLineSeries;
    Panel1: TPanel;
    procedure Button1Click(Sender: TObject);
    procedure CheckBox1Change(Sender: TObject);
    procedure FormCreate(Sender: TObject);
  private
    X, Y: Array of Double;
    DataExtended: Boolean;
    procedure NewData;
    procedure PlotData;

  public

  end;

var
  MainForm: TMainForm;

implementation

{$R *.lfm}

uses
  jpmFunc;

const
  MIN = -10;
  MAX = +10;

procedure TMainForm.Button1Click(Sender: TObject);
begin
  NewData;
  PlotData;
end;

procedure TMainForm.CheckBox1Change(Sender: TObject);
begin
  PlotData;
end;

procedure TMainForm.FormCreate(Sender: TObject);
begin
  RandSeed := 1;
  NewData;
  PlotData;
end;

procedure TMainForm.NewData;
var
  i, n: Integer;
begin
  n := Random(10) + 20;
  SetLength(X, n);
  SetLength(Y, n);
  for i := 0 to n-1 do
  begin
    X[i] := MIN + (MAX - MIN) * i / (n-1);
    Y[i] := Random * 100;
  end;
end;

procedure TMainForm.PlotData;
var
  i, i1, i2, n: Integer;
  xx, yy: Double;
begin
  n := Length(X);
  if Checkbox1.Checked then
  begin
    SetLength(X, n+4);
    SetLength(Y, n+4);
    X[n+3] := X[n-1] + 4;   Y[n+3] := Y[n-1];
    X[n+2] := X[n-1] + 3;   Y[n+2] := Y[n-1];
    X[n+1] := X[n-1] + 2;   Y[n+1] := Y[n-1];
    X[n]   := X[n-1] + 1;   Y[n] := Y[n-1];
    for i := n-1 downto 0 do
    begin
      X[i+1] := X[i];
      Y[i+1] := Y[i];
    end;
    X[0] := X[1] - 1;
    Y[0] := Y[1];
    i1 := 1;
    i2 := n;
    DataExtended := true;
  end else
  begin
    if DataExtended then
    begin
      for i := 1 to n-4 do
      begin
        X[i-1] := X[i];
        Y[i-1] := Y[i];
      end;
      SetLength(X, n-4);
      SetLength(Y, n-4);
      n := Length(X);
      DataExtended := false;
    end;
    i1 := 0;
    i2 := n-1;
  end;

  DataSeries.Clear;
  for i := i1 to i2 do
    DataSeries.AddXY(X[i], Y[i]);

  AkimaSeries.Clear;
  n := 100;
  for i := 0 to n-1 do
  begin
    xx := MIN + (MAX - MIN) * i/(n-1);
    if AkimaInterpolation(X,Y, xx, yy) then
      AkimaSeries.AddXY(xx, yy);
  end;
end;

end.

