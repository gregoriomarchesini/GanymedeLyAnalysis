function shaded_error_bar(x,y,error_bar_values,options)

arguments
    x   (1,:) double
    y   (1,:) double
    error_bar_values (1,:) double
    options.Color    = "r"
    options.Alpha    double  {mustBeGreaterThanOrEqual(options.Alpha,0),mustBeLessThanOrEqual(options.Alpha,1)} =1
    options.EdgeThickness double  =0
    options.EdgeColor = "r"
end

if length(x)~=length(y) || length(x)~=length(error_bar_values)
    error("all input arrays must have consisten size")
end

ax = gca();
hold on
fill(ax ,[x ,fliplr(x)],[y+error_bar_values,fliplr(y-error_bar_values)],options.Color,"FaceAlpha", options.Alpha,'linestyle','none',"HandleVisibility","off")
if   options.EdgeThickness ~=0
  plot(ax ,x,y+error_bar_values,'Color',options.EdgeColor,'linewidth',options.EdgeThickness)
  plot(ax ,x,y-error_bar_values,'Color',options.EdgeColor,'linewidth',options.EdgeThickness)
end



end