close all;
clear all;
addpath("../src")

which_observation = 3010;

% OBSERVATION oe9z03010
if which_observation == 3010
    observation.file = "./oe9z03011_flt.fits";
    observation.x_pcenter = 166;
    observation.y_pcenter = 349;

elseif which_observation == 4010
    observation.file = "./oe9z04011_flt.fits";
    observation.x_pcenter = 158;
    observation.y_pcenter = 361;
else
    error("uknown obseravton")
end

observations_dir   = "../ganymede_observations";
filename           = observation.file;

% extract fits data
GanymedeImage = FitsImageObject;
GanymedeImage.read_image(fullfile(observations_dir,filename));

%%  PLOT IMAGE

[rows,cols]              = size(GanymedeImage.image);     % find image size
x_pixel_range_full_image = 1:cols;                        % full image axis range  x
y_pixel_range_full_image = 1:rows;                        % full images axis range y
max_intensity            = max(max(GanymedeImage.image)); % define max intensity
cscale                   = [0,max_intensity];

% dark counts should remain for the analysis, but you don't want
% to see them in the picture

name ="standard_image.eps";
fig  = figure;
ax   = gca();
ax   = formal_axes(ax);
hold(ax);
imagesc(ax,GanymedeImage.image,cscale);
ax.XLim = [1,cols];
ax.YLim = [1,rows];
color_bar_instance = formal_colorbar(colorbar(ax));
north_pole_angle_direction = 27.0; %deg defined clockwise from horizonatal axis

% define subimage outer edge

x_pmin = 70 ; % pixel
x_pmax = 260; % pixel
y_pmin = 260; % pixel
y_pmax = 450; % pixel

% plot borders
x_square = [x_pmin,x_pmin,x_pmax,x_pmax];
y_square = [y_pmin,y_pmax,y_pmin,y_pmax];
scatter(ax,x_square,y_square,'green','filled',LineWidth=10)

%% RESIZE IMAGE
[ganymede_centred_subimage,sigma_matrix_ganymede_centred_subimage] = GanymedeImage.resize_image(x_pmin,x_pmax,y_pmin,y_pmax);
x_pixel_range_ganymede_centred_subimage     = x_pixel_range_full_image(x_pixel_range_full_image<=x_pmax & x_pixel_range_full_image>= x_pmin);
y_pixel_range_ganymede_centred_subimage     = y_pixel_range_full_image(y_pixel_range_full_image<=y_pmax & y_pixel_range_full_image>= y_pmin);

%% CENTER DEFINITION

% floor(diameter/2) and center by eye (referenced to the full image)
diameter_ganymede    = 70.373;   % define floor(diameter_ganymede/2) so that diameter_ganymede is an odd number (makes easy to have a center than)
x_pixel_center       = observation.x_pcenter ;      % define center of the image by eye on the big image
y_pixel_center       = observation.y_pcenter ;    % define center of the image by eye on the big image
box_radial_extension = 1.5;      % box around the ganymede to be eliminated 
                                 % from the fit expressed is ganymedes radii
                         
x_index_center_ganymede_centred_subimage = find(x_pixel_range_ganymede_centred_subimage == x_pixel_center); % index of the center in the x_range_sub
y_index_center_ganymede_centred_subimage = find(y_pixel_range_ganymede_centred_subimage == y_pixel_center); % index of the center in the y_range_sub

%% CONVERSION FROM COUNTS TO REYLIGHTS

exposition_time =  GanymedeImage.find_key("TEXPTIME"); %s
mx              =  0.0246;     % field of view x dierction [arsec]
my              =  0.0246;     % field of view y dierction [arsec]
A_mirror        =  45238.9342; % cm2

filter_data   = dlmread("HST_STIS_FUV.25MAMA_G140L.dat");
wavelength    = filter_data(:,1);
throughput    = filter_data(:,2);
throughput_Ly = interp1(wavelength,throughput,1216);
A_eff         = A_mirror*throughput_Ly;
Omega         = mx*my*(2*pi/3600/360)^2;

count2KRayleight = 4*pi/10^6/(exposition_time*Omega*A_eff)*10^-3;
ganymede_centred_subimage_reyleights =  ganymede_centred_subimage*count2KRayleight;
sigma_matrix_ganymede_centred_subimage_reyleights = sigma_matrix_ganymede_centred_subimage*count2KRayleight;


% assumed ganymede brightness
ganymede_assumed_brightness = 1.3; % KReyleights

% create grid for bightness extraction
[rows,cols] = size(ganymede_centred_subimage_reyleights);

%create pixel mask for ganymede_centred image with (0,0) at ganymede center
[X_gird_ganymede_centred_subimage,Y_gird_ganymede_centred_subimage] = meshgrid(x_pixel_range_ganymede_centred_subimage- x_pixel_center,y_pixel_range_ganymede_centred_subimage-y_pixel_center);

brightness_mask                             = sqrt(X_gird_ganymede_centred_subimage.^2 + Y_gird_ganymede_centred_subimage.^2) < diameter_ganymede /4;  % only for finding a brightness value take the brightness in half the radius of the ganymede
ganymede_mask_for_ganymede_centred_subimage = sqrt(X_gird_ganymede_centred_subimage.^2 + Y_gird_ganymede_centred_subimage.^2) < diameter_ganymede/2 ;   % mask covering the ganymede disk in ganymede_centred image


mean_brightness_ganymede = mean(mean(ganymede_centred_subimage_reyleights(brightness_mask)));
IPMandGEO                = mean_brightness_ganymede - ganymede_assumed_brightness;
IPMandGEO_counts         = IPMandGEO/count2KRayleight;

% adjust sigma values
sigma_matrix_ganymede_centred_subimage_reyleight   = sqrt(sigma_matrix_ganymede_centred_subimage.^2+IPMandGEO_counts)*count2KRayleight;
mean_brightness_ganymede              = mean_brightness_ganymede - IPMandGEO;
ganymede_centred_subimage_reyleights  = ganymede_centred_subimage_reyleights - IPMandGEO;

n_bins                      = 20;
r_ganymede                  = diameter_ganymede/2;
shift_start                 = 1.0;
shift                       = shift_start+0.2;
average_ring_plot_per_shift = zeros(length(shift),n_bins);
counter                     = 1;
 
% define inner and outer limb limit

r_inner = r_ganymede;
r_outer = r_ganymede*1.2;

for x = -1:1
    for y = -1:1
        [~,~,error_bar] = radial_bins_average(ganymede_centred_subimage_reyleights,x_index_center_ganymede_centred_subimage+x,y_index_center_ganymede_centred_subimage+y,r_inner,r_outer,n_bins,sigma_matrix_ganymede_centred_subimage_reyleight);

        radial_masks=radial_bins_mask(ganymede_centred_subimage_reyleights,x_index_center_ganymede_centred_subimage+x,y_index_center_ganymede_centred_subimage+y,r_inner,r_outer,n_bins);
        average_along_the_bins = zeros(1,n_bins);

        for ii =1:n_bins
            bin_pixels = ganymede_centred_subimage_reyleights(radial_masks(ii).mask);
            average_along_the_bins(ii) =  mean(bin_pixels);
        end

        average_ring_plot_per_shift(counter,:) =  average_along_the_bins;
        
        if counter ==1 
            figure("position",[10,10,1200,400]);
            ax1=formal_axes(subplot(1,2,1));
            hold on;
    
    
            imagesc(x_pixel_range_ganymede_centred_subimage,y_pixel_range_ganymede_centred_subimage,ganymede_centred_subimage_reyleights)
            draw_binned_disk([x_pixel_center  ,y_pixel_center  ],r_inner,r_outer,n_bins,ax1,"Color",[1,1,0])
    
    
            ax1.XLim = [x_pixel_range_ganymede_centred_subimage(1),x_pixel_range_ganymede_centred_subimage(end)];
            ax1.YLim = [y_pixel_range_ganymede_centred_subimage(1),y_pixel_range_ganymede_centred_subimage(end)];
            ax1.XLabel.String = "x [pixel]";
            ax1.YLabel.String = "y [pixel]";
            formal_colorbar(colorbar())
            
            ax2 = formal_axes(subplot(1,2,2));
            hold on;
        end
        
        if x==0 & y==0
         plot(ax2,[radial_masks(1:end).angle]*180/pi,average_along_the_bins,'LineWidth',4)
         errorbar(ax2,[radial_masks(1:end).angle]*180/pi,average_along_the_bins,error_bar)
        end
         plot(ax2,[radial_masks(1:end).angle]*180/pi,average_along_the_bins,'LineWidth',1)
         errorbar(ax2,[radial_masks(1:end).angle]*180/pi,average_along_the_bins,error_bar)


        ax2.XLabel.String = "angle [deg]";
        ax2.YLabel.String = "mean brightness [KR]";
        ax2.YGrid = "on";
        counter = counter+1;

    end
end
