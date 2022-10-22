function draw_circle(center,radius,ax,options)

arguments
    center (1,2) double  
    radius double {mustBePositive} 
    ax     = gca()
    options.LineWidth double {mustBePositive} =1
    options.Color  = "k";
    options.HandleVisibility string {mustBeMember(options.HandleVisibility,["on","off"])}= "off"
end

angle_range = linspace(0,2*pi,100);
x = center(1) + radius*cos(angle_range);
y = center(2) + radius*sin(angle_range);
plot(ax,x,y,Color=options.Color,LineWidth=options.LineWidth,HandleVisibility=options.HandleVisibility)

end