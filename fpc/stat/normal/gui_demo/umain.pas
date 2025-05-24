unit uMain;

{$mode objfpc}{$H+}

interface

uses
  Classes, ComCtrls, ExtCtrls, StdCtrls, SysUtils, Forms, Controls, Graphics,
  Dialogs, TAGraph, TASeries, TAChartUtils, TAGeometry, TADrawUtils;

type
  TMainForm = class(TForm)
    Chart: TChart;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    LineSeries: TLineSeries;
    Panel1: TPanel;
    procedure ChartAfterDraw(ASender: TChart; ADrawer: IChartDrawer);
    procedure ChartMouseLeave(Sender: TObject);
    procedure ChartMouseMove(Sender: TObject; Shift: TShiftState; X, Y: Integer);
    procedure FormCreate(Sender: TObject);
  private
    FMousePos: TPoint;
    procedure CalculateNormalDist;
  public

  end;

var
  MainForm: TMainForm;

implementation

{$R *.lfm}

uses
  jpmStats;

procedure TMainForm.FormCreate(Sender: TObject);
begin
  CalculateNormalDist;
  Label1.Caption := ' ';
  Label2.Caption := ' ';
  Label3.Caption := ' ';
end;

procedure TMainForm.ChartMouseMove(Sender: TObject; Shift: TShiftState;
  X, Y: Integer);
var
  xx: Double;
  yy: Double;
begin
  if not Chart.ScalingValid then
    exit;
  FMousePos := Point(X, Y);
  xx := Chart.XImageToGraph(X);
  yy := NormalDist(xx);
  Label1.Caption := Format('Mouse at x = %.6f', [xx]);
  Label2.Caption := Format('NormalDist(x) = y = %.6f', [yy]);
  Label3.Caption := Format('InvNormalDist(y) = %.6f', [InvNormalDist(yy)]);
  Chart.Invalidate;
end;

procedure TMainForm.ChartMouseLeave(Sender: TObject);
begin
  Label1.Caption := ' ';
  Label2.Caption := ' ';
  Label3.Caption := ' ';
  FMousePos := Point(-999, -999);
end;

procedure TMainForm.ChartAfterDraw(ASender: TChart; ADrawer: IChartDrawer);
var
  ext: TDoubleRect;
  R: TRect;
  x: Double;
  y: Double;
  P: TPoint;
begin
  if FMousePos.X = -999 then
    exit;
  ext := Chart.LogicalExtent;
  R.TopLeft := Chart.GraphToImage(DoublePoint(ext.a.x, ext.b.y));
  R.BottomRight := Chart.GraphToImage(DoublePoint(ext.b.x, ext.a.y));
  if not IsPointInRect(FMousePos, R) then
    exit;
  x := Chart.XImageToGraph(FMousePos.X);
  y := NormalDist(x);
  P := Chart.GraphToImage(DoublePoint(x, y));
  ADrawer.SetPenParams(psDash, clRed, 1);
  ADrawer.Line(R.Left, P.Y, P.X, P.Y);
  ADrawer.Line(P.X, P.Y, P.X, R.Bottom);
end;

procedure TMainForm.CalculateNormalDist;
const
  XMIN = -3;
  XMAX = +3;
  N = 100;
var
  i: Integer;
  x, y: Double;
begin
  LineSeries.Clear;
  for i := 0 to N-1 do
  begin
    x := XMIN + (XMAX - XMIN) * i / (N-1);
    y := NormalDist(x);
    LineSeries.AddXY(x, y);
  end;
end;

end.

