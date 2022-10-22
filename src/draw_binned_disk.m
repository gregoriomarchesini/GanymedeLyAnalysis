function draw_binned_disk(center,inner_radius,outer_radius,n_bins,ax,options)

arguments
    
    center (1,2) double  ;
    inner_radius double {mustBePositive} ;
    outer_radius double {mustBePositive} ;
    n_bins double {mustBePositive} ;
    ax     = gca();
    options.Color  = "k";
    options.LineWidth double {mustBePositive} = 1;
    options.HandleVisibility string {mustBeMember(options.HandleVisibility,["on","off"])}= "off";
end

if nargin <5
    error("Not enough inputs")
end
angle_range=linspace(-pi,pi,n_bins);
angle_range=angle_range(1:end-1);

hold on;
for angle=angle_range
    x =  center(1) + [inner_radius outer_radius]*cos(angle);
    y =  center(2) + [inner_radius outer_radius]*sin(angle);
    
    plot(x,y,"LineWidth",options.LineWidth,"color",options.Color,"HandleVisibility",options.HandleVisibility);
end
draw_circle(center,inner_radius,ax,Color=options.Color,LineWidth=options.LineWidth,HandleVisibility=options.HandleVisibility)
draw_circle(center,outer_radius,ax,Color=options.Color,LineWidth=options.LineWidth,HandleVisibility=options.HandleVisibility)
end