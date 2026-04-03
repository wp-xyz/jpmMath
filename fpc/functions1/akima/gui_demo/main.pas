unit main;

{$mode objfpc}{$H+}

interface

uses
  Classes, ExtCtrls, Spin, StdCtrls, SysUtils, Forms, Controls, Graphics,
  Dialogs, TAGraph, TASeries,
  jpmInterpolation;

type
  TMainForm = class(TForm)
    btnNewData: TButton;
    Chart: TChart;
    AkimaSeries: TLineSeries;
    DataSeries: TLineSeries;
    lblStartIndex: TLabel;
    lblEndIndex: TLabel;
    lblInfo: TLabel;
    Panel1: TPanel;
    rbRandom: TRadioButton;
    rbSin: TRadioButton;
    seStartIndex: TSpinEdit;
    seEndIndex: TSpinEdit;
    procedure btnNewDataClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure seStartIndexChange(Sender: TObject);
    procedure seEndIndexChange(Sender: TObject);
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
  MIN = -10;
  MAX = +10;

procedure TMainForm.btnNewDataClick(Sender: TObject);
begin
  NewData;
  PlotData;
end;

procedure TMainForm.FormCreate(Sender: TObject);
begin
  RandSeed := 1;
  NewData;
  PlotData;
end;

procedure TMainForm.seStartIndexChange(Sender: TObject);
begin
  PlotData;
end;

procedure TMainForm.seEndIndexChange(Sender: TObject);
begin
  PlotData;
end;

procedure TMainForm.NewData;
var
  i, n: Integer;
begin
  n := Random(10) + 20;
  SetLength(X, n);
  SetLength(Y, n);
  for i := 0 to n - 1 do
  begin
    X[i] := MIN + (MAX - MIN) * i / (n-1);
    if rbRandom.Checked then
      Y[i] :=  Random * 100.0
    else
      Y[i] := (sin(X[i]) + 1) * 50;
  end;
end;

procedure TMainForm.PlotData;
var
  i, n: Integer;
  xx, yy: Double;
  err: String;
begin
  DataSeries.Clear;
  for i := 0 to High(X) do
    DataSeries.AddXY(X[i], Y[i], IntToStr(i));

  n := 100;  // Number of interpolated data points
  AkimaSeries.Clear;
  for i := 0 to n-1 do
  begin
    xx := MIN + (MAX - MIN) * i / (n-1);
    yy := AkimaInterpolation(xx, X, Y, seStartIndex.Value, seEndIndex.Value, err);
    if err = '' then
      AkimaSeries.AddXY(xx, yy);
  end;
end;

end.

