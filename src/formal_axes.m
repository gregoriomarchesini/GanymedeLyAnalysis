function [ax]=formal_axes(ax)
% predefined set of specifications for plots

if nargin == 0
    ax = gca();
end

ax.XLabel.Interpreter  = 'latex';
ax.XLabel.FontSize     = 14;
ax.XLabel.FontWeight   = "Bold";


ax.YLabel.FontSize     = 14;
ax.YLabel.FontWeight   = "Bold";

ax.Title.FontSize     = 16;
ax.Title.FontWeight   = "Bold";



ax.GridColor          = [0,0,0];
ax.GridAlpha          = 0.1;
ax.GridLineStyle      = "--";
ax.Box                = "On";
ax.LineWidth          = 1.3;
ax.XMinorTick         = "on";
ax.YMinorTick         = "on";
ax.Layer              = "top";

end
