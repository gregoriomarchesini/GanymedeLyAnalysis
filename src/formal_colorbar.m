function [colorbar]=formal_colorbar(colorbar)
% predefined set of specifications for plots

colorbar.XLabel.Interpreter = "latex";
colorbar.XLabel.FontSize    = 14;
colorbar.XLabel.FontWeight  = "bold";

colorbar.YLabel.Interpreter = "latex";
colorbar.YLabel.FontSize    = 14;
colorbar.YLabel.FontWeight  = "bold";

colorbar.TickLabelInterpreter = 'latex'

colorbar.Box                = "On";
colorbar.LineWidth          = 1.3;
colorbar.Location           = "eastoutside";

end