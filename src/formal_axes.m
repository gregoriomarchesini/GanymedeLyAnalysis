function [ax]=formal_axes(ax)
% predefined set of specifications for plots

if nargin == 0
    ax = gca();
end

ax.XLabel.FontSize     = 12;
ax.YLabel.FontSize     = 12;

ax.Title.FontSize     = 13;
ax.GridColor          = [0,0,0];
ax.GridAlpha          = 0.1;
ax.GridLineStyle      = "--";
ax.Box                = "On";
ax.LineWidth          = 1.3;
ax.XMinorTick         = "on";
ax.YMinorTick         = "on";
ax.Layer              = "top";

end
