close all;
clear all;
addpath("../src")
addpath("../utility_files")

% THIS SCRIPT IS USED TO START THE WHOLE ANALYSIS FRO EACH OBSERVATION 
% INDEPENDENTLY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN USER DEFINED SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

which_observation = 3010; % DEFINE THE OBSERVATION THAT YOU WANT TO SEE

if which_observation == 3010
    observation.file = "./oe9z03011_flt.fits";
    observation.x_pcenter  = 166;       % pixels
    observation.y_pcenter  = 349;       % pixels
    observation.best_n0    = 820;       % 1/cm^2
    observation.zc         = 5.22e-14 ; % cm^2
    observation.g_factor   = 7e-14;
    observation.poly_order = 2;
    obseravtion.g_factor   = 7e-14;     % kR
    observation.north_pole_angle_direction = 27.0; %deg defined anti-clockwise from horizonatal axis

elseif which_observation == 4010
    observation.file = "./oe9z04011_flt.fits";
    observation.x_pcenter  = 158;       % pixels
    observation.y_pcenter  = 361;       % pixels
    observation.best_n0    = 1040;      % 1/cm^2
    observation.zc         = 5.22e-14 ; % cm^2
    observation.g_factor   = 7e-14;
    observation.poly_order = 3;
    observation.north_pole_angle_direction = 24.4; %deg defined anti-clockwise from horizonatal axis
else
    error("uknown obseravton")
end


observations_dir    = "../ganymede_observations";
filename             = observation.file;
output_dir           = "./images";
output_dir_models    = "./images/models";

show_model           = false; % TRUE IF YOU WANT TO SEE THE BACKGROUND MODEL IMAGES
show_giorno          = false; % TRUE IF YOU WANT TO SEE THAT THE GIORNO METJOD GIVES THE CORRECT CENTER
show_corona_analysis = true;  % TRUE IF YOU WANT TO SEE THE CORONA ANLYSIS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTOMATIC SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
ax   = formal_axes(gca());
hold(ax);
imagesc(ax,GanymedeImage.image,cscale);
ax.XLim = [1,cols];
ax.YLim = [1,rows];
color_bar_instance = formal_colorbar(colorbar(ax));


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
ganymede_centred_subimage_reyleights = ganymede_centred_subimage*count2KRayleight;
sigma_matrix_ganymede_centred_subimage_reyleights = sigma_matrix_ganymede_centred_subimage*count2KRayleight;


%% FIT THE MODEL
% assumed ganymede brightness
ganymede_assumed_brightness = 1.3; % KReyleights

% decide max order for background polynomial surface
poly_order = observation.poly_order;

syms n0
syms a [poly_order+1,poly_order+1] 
syms x y

polynomial_surface = 0;
for ii=0:poly_order
    for jj =0:poly_order-ii
        polynomial_surface = polynomial_surface + a(ii+1,jj+1)*x^ii*y^jj;
    end
end

zc        = observation.zc;
g_factor  = obseravtion.g_factor;
Rg        = 2643.1e5;
R_ganymede_pixel  = diameter_ganymede/2;

% H corona absorbtion model outside ganymede disk
Nh  = n0*Rg*(R_ganymede_pixel./sqrt(x.^2+y.^2)).*pi;
tau = zc*Nh;
T   = exp(-tau);

% H corona emission model outside ganymede disk
I_H_corona_emission = Nh*g_factor;


% Define point spread function
PSF = imread("fuvmama_1216_00.fits"); % point spread function
[rows_psf,cols_psf] = size(PSF);

% Expand axis definition to eliminate noise after convolution
right_expansion =  max(x_pixel_range_ganymede_centred_subimage)+1:max(x_pixel_range_ganymede_centred_subimage)+cols_psf*2;   % 2 is just a sefty margin
left_expansion  =  (min(x_pixel_range_ganymede_centred_subimage)-cols_psf*2):min(x_pixel_range_ganymede_centred_subimage)-1; % 2 is just a sefty margin
up_expansion    =  max(y_pixel_range_ganymede_centred_subimage)+1:max(y_pixel_range_ganymede_centred_subimage)+rows_psf*2;   % 2 is just a sefty margin
down_expansion  =  (min(y_pixel_range_ganymede_centred_subimage)-cols_psf*2):min(y_pixel_range_ganymede_centred_subimage)-1; % 2 is just a sefty margin

x_expanded_axis = [left_expansion,x_pixel_range_ganymede_centred_subimage,right_expansion];
y_expanded_axis = [down_expansion,y_pixel_range_ganymede_centred_subimage,up_expansion];
x_original_mask = boolean([zeros(size(left_expansion)),ones(size(x_pixel_range_ganymede_centred_subimage)),zeros(size(right_expansion))]);
y_original_mask = boolean([zeros(size(down_expansion)),ones(size(y_pixel_range_ganymede_centred_subimage)),zeros(size(up_expansion))]);


% create grid for bightness extraction
[rows,cols] = size(ganymede_centred_subimage_reyleights);

%create pixel mask for ganymede_centred image with (0,0) at ganymede center
[X_gird_ganymede_centred_subimage,Y_gird_ganymede_centred_subimage] = meshgrid(x_pixel_range_ganymede_centred_subimage- x_pixel_center,y_pixel_range_ganymede_centred_subimage-y_pixel_center);

%create pixel mask for expanded image
[X_grid_model,Y_grid_model] = meshgrid(x_expanded_axis- x_pixel_center ,y_expanded_axis- y_pixel_center );

brightness_mask                          = sqrt(X_gird_ganymede_centred_subimage.^2 + Y_gird_ganymede_centred_subimage.^2) < R_ganymede_pixel/2; % only for finding a brightness value take the brightness in half the radius of the ganymede
ganymede_mask_for_ganymede_centred_subimage = sqrt(X_gird_ganymede_centred_subimage.^2 + Y_gird_ganymede_centred_subimage.^2) < R_ganymede_pixel;   % mask covering the ganymede disk in ganymede_centred image
ganymede_mask_for_model_image            = sqrt(X_grid_model.^2 + Y_grid_model.^2) < R_ganymede_pixel;   % this the mask covering the ganymede in the expanded image

mean_brightness_ganymede = mean(mean(ganymede_centred_subimage_reyleights(brightness_mask)));
IPMandGEO                = mean_brightness_ganymede - ganymede_assumed_brightness;
IPMandGEO_counts         = IPMandGEO/count2KRayleight;

% adjust sigma values
sigma_matrix_ganymede_centred_subimage_reyleight   = sqrt(sigma_matrix_ganymede_centred_subimage.^2+IPMandGEO_counts)*count2KRayleight;

mean_brightness_ganymede  = mean_brightness_ganymede - IPMandGEO;

ganymede_centred_subimage_reyleights = ganymede_centred_subimage_reyleights - IPMandGEO;

% obtain polynomial fit of the image
% the x and y axis must be centered with zero at the ganymede center.
% Otherwise the transmission model won't work
[x_fit,y_fit,z_fit,weights_fit] = prepareSurfaceData(X_gird_ganymede_centred_subimage,Y_gird_ganymede_centred_subimage,ganymede_centred_subimage_reyleights,1./sigma_matrix_ganymede_centred_subimage_reyleight.^2);

exclude_set = ~excludedata(x_fit,y_fit,'box',[-box_radial_extension*R_ganymede_pixel ,...
                                              +box_radial_extension*R_ganymede_pixel ,...
                                              -box_radial_extension*R_ganymede_pixel ,...
                                              + box_radial_extension*R_ganymede_pixel]);




% backgroumd surface model
total_model_sym     = polynomial_surface.*T  + I_H_corona_emission; % background surface plus the Transmission
total_model         = matlabFunction(total_model_sym);
fit_parameters      = string(symvar(total_model_sym));

% eliminate variables that are not to be fit in the model
eliminate_xyn0_mask = string(symvar(total_model_sym)) ~= 'x' & ...
                      string(symvar(total_model_sym)) ~= 'y' &  ...
                      string(symvar(total_model_sym)) ~= 'n0';

fit_parameters      = cellstr(fit_parameters(eliminate_xyn0_mask));

myfittype = fittype(total_model,...
                     'dependent',{'z'},'independent',{'x','y'},...
                     'coefficients',fit_parameters ,'problem','n0');


options = fitoptions(myfittype);
coefficients_upperbound = +10000;
coefficients_lowerbound = -10000;

options.Upper      = ones(1,length(fit_parameters))*coefficients_upperbound;
options.Lower      = ones(1,length(fit_parameters))*coefficients_lowerbound;
mid_point_bound    = (coefficients_upperbound+coefficients_lowerbound)/2;
options.Exclude    = exclude_set;
options.StartPoint = zeros(1,length(fit_parameters))*mid_point_bound;
% options.Weights    = weights_fit;

counter = 1;
decimation_error_bar = 5;


for n0_value=observation.best_n0
   
    model_fit = fit([x_fit,y_fit],z_fit,myfittype,options,"problem",n0_value);
    % create model image from fit
    bkg_only_model_image = double(model_fit(X_grid_model,Y_grid_model));
   
    
    % Check quality of the fit and excluded points
    if counter ==1
        figure
        ax=formal_axes(gca());
        plot(model_fit,[x_fit,y_fit],z_fit,'Exclude',exclude_set)
        ax.Title.String="Check FIT points";
        ax.XLabel.String="x [pixels]";
        ax.YLabel.String="y [pixels]";
    end
    
    % assume constant brightness at disk as constant
    row_model_with_ganymede = bkg_only_model_image;
    row_model_with_ganymede(ganymede_mask_for_model_image)  = mean_brightness_ganymede;  
    
    final_model_image_PSF                         = conv2(row_model_with_ganymede,PSF,'same');% apply PSF
    final_model_image_PSF                         = final_model_image_PSF(y_original_mask,x_original_mask);
    rsw_model_image                               = row_model_with_ganymede(y_original_mask,x_original_mask); % image size before the PSF convolution is applied
    
    
    % ASEES MODEL ON THE BACKGROUND
    if counter==1 
        figure("Position",[1,1,1200,700]);
        big_ax = subplot(4,2,1:4);
        hold on
        imagesc(big_ax,final_model_image_PSF);
        [rows,cols] =size(final_model_image_PSF); 
        line_box_width = 3;
        
        % lower horizontal definition
        ymin_box1 = 10;
        ymax_box1 = 40;
        ymid_point_box1 = (ymin_box1 + ymax_box1)/2;
        yline(ymin_box1,"-r","LineWidth",line_box_width);
        yline(ymax_box1,"-r","LineWidth",line_box_width);
        
        xlim([1,cols])
        ylim([1,rows])
        t = annotation("textbox");
        t.Parent = gca();
        t.String= "BOX 1";
        t.LineStyle = "none";
        t.Position=[x_index_center_ganymede_centred_subimage,ymid_point_box1,10,10];
    
    
        % upper horizontal box
        ymin_box2 = 150;
        ymax_box2 = 180;
        ymid_point_box2 = (ymin_box2+ ymax_box2)/2;
        yline(ymin_box2,"-r","LineWidth",line_box_width);
        yline(ymax_box2,"-r","LineWidth",line_box_width);
    
        t = annotation("textbox");
        t.Parent = gca();
        t.String= "BOX 2";
        t.LineStyle = "none";
        
        t.Position=[x_index_center_ganymede_centred_subimage,ymid_point_box2,10,10];
        
        % right vertical definition
        xmin_box3 = 10;
        xmax_box3 = 40;
        xmid_point_box1 = (xmin_box3+ xmax_box3)/2;
        xline(xmin_box3,"-r","LineWidth",line_box_width);
        xline(xmax_box3,"-r","LineWidth",line_box_width);
    
        t = annotation("textbox");
        t.Parent = gca();
        t.String= "BOX 3";
        t.LineStyle = "none";
        t.Position=[xmid_point_box1,y_index_center_ganymede_centred_subimage,10,10];
        
        % left vertical box
        xmin_box4 = 150;
        xmax_box4 = 180;
        xmid_point_box2 = (xmin_box4+ xmax_box4)/2;
        xline(xmin_box4,"-r","LineWidth",line_box_width);
        xline(xmax_box4,"-r","LineWidth",line_box_width);
    
        t = annotation("textbox");
        t.Parent = gca();
        t.String= "BOX 4";
        t.LineStyle = "none";
        t.Position=[xmid_point_box2,y_index_center_ganymede_centred_subimage,10,10];
        xlabel("pixel")
        ylabel("pixel")

    end 

    if counter==1
      box1_ax = subplot(4,2,5);
      hold on
      plot(box1_ax,mean(ganymede_centred_subimage_reyleights(ymin_box1:ymax_box1,:)),"DisplayName","STIS observation")  
      title("BOX 1")
      xlabel("pixel")
      ylabel("mean pixel value")
      
    end
    plot(box1_ax,mean(final_model_image_PSF(ymin_box1:ymax_box1,:)),"DisplayName",sprintf("n_0 = %.f",n0_value));
    

    if counter ==1
      box2_ax=subplot(4,2,6);
      hold on
      plot(box2_ax,mean(ganymede_centred_subimage_reyleights(ymin_box2:ymax_box2,:)),"DisplayName","STIS observation")
      title("BOX 2")
      xlabel("pixel")
      ylabel("mean pixel value")
    end
    plot(box2_ax,mean(final_model_image_PSF(ymin_box2:ymax_box2,:)),"DisplayName",sprintf("n_0 = %.f",n0_value));

    
    
    if counter ==1 
      box3_ax=subplot(4,2,7);
      hold on
      plot(box3_ax,mean(ganymede_centred_subimage_reyleights(:,xmin_box3:xmax_box3),2),"DisplayName","STIS observation")
      title("BOX 3")
      xlabel("pixel")
      ylabel("mean pixel value")
    end
    plot(box3_ax,mean(final_model_image_PSF(:,xmin_box3:xmax_box3),2),"DisplayName",sprintf("n_0 = %.f",n0_value));
    



    if counter ==1
       box4_ax= subplot(4,2,8);
       hold on
       plot(box4_ax,mean(ganymede_centred_subimage_reyleights(:,xmin_box4:xmax_box4),2),"DisplayName","STIS observation")
       title("BOX 4")
       xlabel("pixel")
       ylabel("mean pixel value")
    end
    plot(box4_ax,mean(final_model_image_PSF(:,xmin_box4:xmax_box4),2),"DisplayName",sprintf("n_0 = %.f",n0_value));
    sgtitle(sprintf("polynomail order : %i",poly_order))
    
    % compute chi2 at the rim only
    R_gird = sqrt(X_gird_ganymede_centred_subimage.^2 + Y_gird_ganymede_centred_subimage.^2);
    rim_mask = R_gird >=R_ganymede_pixel & R_gird <= R_ganymede_pixel*1.5;
    chi2_rim = sum(sum((final_model_image_PSF(rim_mask)-ganymede_centred_subimage_reyleights(rim_mask)).^2./sigma_matrix_ganymede_centred_subimage_reyleight(rim_mask).^2))/numel(ganymede_centred_subimage_reyleights(rim_mask));

    if show_model
        figure("Position",[100,100,1300,400])
        ax  = formal_axes(subplot(1,2,1));
        hold on 
        imagesc(x_pixel_range_ganymede_centred_subimage,y_pixel_range_ganymede_centred_subimage,final_model_image_PSF )
        ax.XLim               = [x_pmin , x_pmax];
        ax.YLim               = [y_pmin , y_pmax];
        ax.XLabel.String      = "pixel";
        ax.YLabel.String      = "pixel";
        ax.Title.String =sprintf("bkg+ganymede RAW  n_0 value : %.2f  Chi2 : %.4f",n0_value,chi2);
        colobar_instance                = colorbar(ax);
        colobar_instance                = formal_colorbar(colobar_instance);
        clim(ax,[0,max(max(ganymede_centred_subimage_reyleights))])
        colobar_instance.YLabel.String  = "[kR]"';
    
        name = "standard_reyleights.eps";
        
        ax  = formal_axes(subplot(1,2,2));
        ax  = formal_axes(ax);
        hold on
    
    
        imagesc(x_pixel_range_ganymede_centred_subimage,y_pixel_range_ganymede_centred_subimage,ganymede_centred_subimage_reyleights)
        ax.XLabel.String      = "pixel";
        ax.YLabel.String      = "pixel";
        ax.Title.String       = "Ganymede brightness";
        ax.XLim               = [x_pmin , x_pmax];
        ax.YLim               = [y_pmin , y_pmax];
    
    
        colobar_instance                = colorbar(ax);
        colobar_instance                = formal_colorbar(colobar_instance);
        clim(ax,[0,max(max(ganymede_centred_subimage_reyleights))])
        colobar_instance.YLabel.String  = "[kR]"';
    
    end
    
    
    
    %% GIORNO CENTERING METHOD
    R_ganymede_pixel; 
    
    
    inner_radius   = R_ganymede_pixel;
    outer_radius   = R_ganymede_pixel *1.2;
    number_of_bins = 20; 
    [best_col_index_center,best_row_index_center,epsilon_struct] = Giorno_method(ganymede_centred_subimage_reyleights,...
                                                                                 final_model_image_PSF,...
                                                                                 x_index_center_ganymede_centred_subimage,...
                                                                                 y_index_center_ganymede_centred_subimage,...
                                                                                 inner_radius,...
                                                                                 outer_radius, ...
                                                                                 number_of_bins);
    if show_giorno 
        figure
        hold on
        imagesc(epsilon_struct.x_shift,epsilon_struct.y_shift,epsilon_struct.giorno_matrix);
        ax = formal_axes(gca());
        ax.XLim = [epsilon_struct.x_shift(1),epsilon_struct.x_shift(end)];
        ax.YLim = [epsilon_struct.y_shift(1),epsilon_struct.y_shift(end)];
        ax.XLabel.String="offset from assumed center [pixel]";
        ax.YLabel.String="offset from assumed center [pixel]";
        h  = formal_colorbar(colorbar);
        h.YLabel.String = "$\epsilon$";

    end
        
    %% GENERATE RADIAL PLOTS
    
    thikness_anulus = 1;
    if counter ==1 
      ganymede_centred_subimage_vertical_box_plot   = generate_vertical_rectangle_plot(ganymede_centred_subimage_reyleights,x_index_center_ganymede_centred_subimage,diameter_ganymede);
      ganymede_centred_subimage_horizontal_box_plot = generate_horizontal_rectangle_plot(ganymede_centred_subimage_reyleights,y_index_center_ganymede_centred_subimage,diameter_ganymede);
      [ganymede_centred_subimage_radial_plot,~]     = generate_radial_plot(ganymede_centred_subimage_reyleights,x_index_center_ganymede_centred_subimage,y_index_center_ganymede_centred_subimage,R_ganymede_pixel*2.5,thikness_anulus);
    end

    [model_image_vertical_box_plot,error_bar_vertical_box_plot]     = generate_vertical_rectangle_plot(final_model_image_PSF,x_index_center_ganymede_centred_subimage,diameter_ganymede,sigma_matrix_ganymede_centred_subimage_reyleight );
    [model_image_horizontal_box_plot,error_bar_horizontal_box_plot] = generate_horizontal_rectangle_plot(final_model_image_PSF,y_index_center_ganymede_centred_subimage,diameter_ganymede,sigma_matrix_ganymede_centred_subimage_reyleight );
    [model_image_radial_plot,radial_range,error_bar_radial_plot]    = generate_radial_plot(final_model_image_PSF,x_index_center_ganymede_centred_subimage,y_index_center_ganymede_centred_subimage,R_ganymede_pixel*2.5,thikness_anulus,sigma_matrix_ganymede_centred_subimage_reyleight );
    
    chi2_rad     = sum((model_image_radial_plot         - ganymede_centred_subimage_radial_plot        ).^2./error_bar_radial_plot.^2  )/numel(model_image_radial_plot);
    chi2_vert    = sum((model_image_vertical_box_plot   - ganymede_centred_subimage_vertical_box_plot  ).^2./error_bar_vertical_box_plot.^2  )/numel(model_image_vertical_box_plot  );
    chi2_hor     = sum((model_image_horizontal_box_plot - ganymede_centred_subimage_horizontal_box_plot).^2./error_bar_horizontal_box_plot.^2 )/numel(model_image_horizontal_box_plot );
    
   
    if counter == 1
       boxin_fig = figure("Position",[1,1,600,800]);
       boxing_ax1 = formal_axes(subplot(3,1,1));
       hold on
       plot(boxing_ax1,(y_pixel_range_ganymede_centred_subimage-y_pixel_center)/R_ganymede_pixel,ganymede_centred_subimage_vertical_box_plot,"DisplayName","STIS observation")
       errorbar(boxing_ax1,(y_pixel_range_ganymede_centred_subimage((1:decimation_error_bar:end))-y_pixel_center)/R_ganymede_pixel,ganymede_centred_subimage_vertical_box_plot((1:decimation_error_bar:end)),error_bar_vertical_box_plot((1:decimation_error_bar:end)),'o',"HandleVisibility","off")
       boxing_ax1.XLabel.String = "$r_{\mathcal{G}}$";
       boxing_ax1.YLabel.String = "Brightness $[kR]$";
       boxing_ax1.Title.String = "Vertical box around the ganymede";
    end
    
    plot(boxing_ax1,(y_pixel_range_ganymede_centred_subimage-y_pixel_center)/R_ganymede_pixel,model_image_vertical_box_plot,"DisplayName",sprintf("n0 = %.2f, chi= %.4f",[n0_value,chi2_vert]))
   
    
    if counter == 1
       boxing_ax2 = formal_axes(subplot(3,1,2));
       hold on
       plot(boxing_ax2,(x_pixel_range_ganymede_centred_subimage-x_pixel_center)/R_ganymede_pixel,ganymede_centred_subimage_horizontal_box_plot,"DisplayName","STIS observation")
       errorbar(boxing_ax2,(x_pixel_range_ganymede_centred_subimage((1:decimation_error_bar:end))-x_pixel_center)/R_ganymede_pixel,ganymede_centred_subimage_horizontal_box_plot((1:decimation_error_bar:end)),error_bar_horizontal_box_plot((1:decimation_error_bar:end)),'o',"HandleVisibility","off")
       boxing_ax2.XLabel.String = "$r_{\mathcal{G}}$";
       boxing_ax2.YLabel.String = "Brightness $[kR]$";
       boxing_ax2.Title.String = "Horizontal box around the ganymede";
    end
    plot(boxing_ax2,(x_pixel_range_ganymede_centred_subimage-x_pixel_center)/R_ganymede_pixel,model_image_horizontal_box_plot,"DisplayName",sprintf("n0 = %.2f, chi= %.4f",[n0_value,chi2_hor]))
    
    
    if counter ==1
        boxing_ax3 = formal_axes(subplot(3,1,3));
        hold on
        plot(boxing_ax3,radial_range/R_ganymede_pixel,ganymede_centred_subimage_radial_plot,"DisplayName","STIS observation")
        errorbar(boxing_ax3,radial_range((1:decimation_error_bar:end))/R_ganymede_pixel,ganymede_centred_subimage_radial_plot((1:decimation_error_bar:end)),error_bar_radial_plot((1:decimation_error_bar:end)),'o',"HandleVisibility","off")
        boxing_ax3.XLabel.String = "$r_{\mathcal{G}}$";
        boxing_ax3.YLabel.String = "Brightness $[kR]$";
        boxing_ax3.Title.String = "Radial average value";
    end
    plot(boxing_ax3,radial_range/R_ganymede_pixel,model_image_radial_plot,"DisplayName",sprintf("n0 = %.2f, chi= %.4f",[n0_value,chi2_rad]))
    
    
    % Radial analysis of the images 
   
    inner_radius = R_ganymede_pixel*1.0;
    outer_radius = R_ganymede_pixel*1.2;
    number_of_bins   = 24; 
    xcenter          = x_index_center_ganymede_centred_subimage;    
    ycenter          = y_index_center_ganymede_centred_subimage;    
      
    difference_image = ganymede_centred_subimage_reyleights-final_model_image_PSF;
    [mean_along_bin_ganymede_centred_subimage,center_angles_range_dense,error_per_bin] = radial_bins_average(ganymede_centred_subimage_reyleights,xcenter,ycenter,inner_radius,outer_radius,number_of_bins,sigma_matrix_ganymede_centred_subimage_reyleight);
    [mean_along_bin_model_image,center_angles_range]   = radial_bins_average(final_model_image_PSF,xcenter,ycenter,inner_radius,outer_radius,number_of_bins);
    


    mean_normalization_factor =mean(mean_along_bin_ganymede_centred_subimage)/mean(mean_along_bin_model_image);
    mean_along_bin_model_image = mean_along_bin_model_image* mean_normalization_factor;
    center_angles_range = center_angles_range*180/pi-observation.north_pole_angle_direction; % center north pole at zero
    center_angles_range_dense = center_angles_range_dense*180/pi-observation.north_pole_angle_direction;
    
    if counter ==1 
        fig=figure("Position",[1,1,1200,400]);
        rim_plot_ax4 = formal_axes(subplot(2,1,1));
        hold on
        plot(rim_plot_ax4,center_angles_range_dense,mean_along_bin_ganymede_centred_subimage,"-ob","DisplayName","STIS observation")
        fill(rim_plot_ax4,[center_angles_range_dense ,fliplr(center_angles_range_dense )],[mean_along_bin_ganymede_centred_subimage+error_per_bin,fliplr(mean_along_bin_ganymede_centred_subimage-error_per_bin)],[1,0,0],"FaceAlpha",0.1,'linestyle','none',"HandleVisibility","off")
        xline(0,"-k","North",LabelVerticalAlignment="top",HandleVisibility="off")
        xline(-90,"-k","East",LabelVerticalAlignment="top",HandleVisibility="off")
        xline(90,"-k","West",LabelVerticalAlignment="top",HandleVisibility="off")
        xline(-180,"-k","South",LabelVerticalAlignment="top",HandleVisibility="off")
        rim_plot_ax4.XLabel.String = "Angle from North $[deg]$";
        rim_plot_ax4.YLabel.String = "Brightness $[kR]$";

        rim_plot_ax5 = formal_axes(subplot(2,1,2));
        hold on
        rim_plot_ax5.YLabel.String="$\frac{Obs-Mod}{\sigma}$";
        rim_plot_ax5.XLabel.String="Angle from North $[deg]$";
        xline(0,"-k","North",LabelVerticalAlignment="top",HandleVisibility="off")
        xline(-90,"-k","East",LabelVerticalAlignment="top",HandleVisibility="off")
        xline(90,"-k","West",LabelVerticalAlignment="top",HandleVisibility="off")
        xline(-180,"-k","South",LabelVerticalAlignment="top",HandleVisibility="off")
        rim_plot_ax5.YLim = [-3,3];

    end
    plot(rim_plot_ax4,center_angles_range,mean_along_bin_model_image,"-ok","DisplayName","normalised model")
    plot(rim_plot_ax5,center_angles_range,(mean_along_bin_ganymede_centred_subimage-mean_along_bin_model_image)./error_per_bin,"-ok");
        
    counter = counter +1;

 end


legend(rim_plot_ax4,"location","southeast");
legend(box4_ax,"location","bestoutside")
legend(box2_ax,"location","bestoutside")
legend(box1_ax,"location","bestoutside")
legend(box3_ax,"location","bestoutside")

    
    

