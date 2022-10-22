function ax = doubleXAxis(targetAxes)
ax = axes();
ax.Units = "normalized";

targetAxes.Units = "normalized";
ax.Position = targetAxes.Position - [0.,0.03,0,0];
ax.Color = "none";
ax.YAxis.Visible = "off";
ax.XMinorTick = "on"
ax.LineWidth = 2;

end
