function ax = doubleXAxis(targetAxes)
% Gregorio Marchesini 
% 22.10.2022
%
% Creates second X axis below a given target axes object
% Returns an axes object that can be edited as any other normal axes object

ax = axes();
ax.Units = "normalized";

targetAxes.Units = "normalized";
ax.Position = targetAxes.Position - [0.,0.03,0,0];
ax.Color = "none";
ax.YAxis.Visible = "off";
ax.XMinorTick = "on"
ax.LineWidth = 2;

end
