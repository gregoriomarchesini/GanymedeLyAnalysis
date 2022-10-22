close all;
clear all;
addpath("./src")
addpath("../utility_files")

% THIS SCRIPT IS USED TO FIND THE OPTIMAL CHI@ MEASURE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN USER DEFINED SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

which_observation = 3010;
music_intretainment = true;


% OBSERVATION oe9z03010
if which_observation   == 3010
    observation.file      = "./oe9z03011_flt.fits";
    observation.name      = 'oe9z03011';
    observation.x_pcenter = 166;
    observation.y_pcenter = 349;
    observation.zc        = 5.22e-14 ;
    obseravtion.g_factor  = 7e-14;

elseif which_observation == 4010
    observation.file = "./oe9z04011_flt.fits";
    observation.name = 'oe9z04011';
    observation.x_pcenter = 158;
    observation.y_pcenter = 361;
    observation.zc        = 5.22e-14 ;
    obseravtion.g_factor  = 7e-14;
else
    error("uknown obseravton")
end

% define here the polynomail order range that you want to check and 
% the no range that you want to optimise over

min_poly_order = 2;
max_poly_order = 2;
min_n0         = 4000;
max_n0         = 5000;
n0_step        = 10


% EASTER EGG TO BE REMOVED
% type :> clear sound         if you want to stop it 

if music_intretainment 
   [y,Fs] = audioread('german-national-anthem.mp3');
   sound(y,Fs)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTOMATIC SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poly_surface_order_range = 2:-1:2;
n0_range                 = min_n0 :n0_step:max_n0;
poly_orders_string       = string(poly_surface_order_range);

observations_dir   = "../ganymede_observations";
filename           = observation.file;
output_dir         = "../images";


% extract fits data
GanymedeImage = FitsImageObject;
GanymedeImage.read_image(fullfile(observations_dir,filename));

%%  PLOT IMAGE

[rows,cols]              = size(GanymedeImage.image);     % find image size
x_pixel_range_full_image = 1:cols;                        % full image axis range  x
y_pixel_range_full_image = 1:rows;                        % full images axis range y
max_intensity            = max(max(GanymedeImage.image)); % define max intensity
cscale                   = [0,max_intensity];



% define subimage outer edge

x_pmin = 70 ; % pixel
x_pmax = 260; % pixel
y_pmin = 260; % pixel
y_pmax = 450; % pixel

% plot borders
x_square = [x_pmin,x_pmin,x_pmax,x_pmax];
y_square = [y_pmin,y_pmax,y_pmin,y_pmax];


%% RESIZE IMAGE
[ganymede_centred_subimage,sigma_matrix_ganymede_centred_subimage] = GanymedeImage.resize_image(x_pmin,x_pmax,y_pmin,y_pmax);
x_pixel_range_ganymede_centred_subimage     = x_pixel_range_full_image(x_pixel_range_full_image<=x_pmax & x_pixel_range_full_image>= x_pmin);
y_pixel_range_ganymede_centred_subimage     = y_pixel_range_full_image(y_pixel_range_full_image<=y_pmax & y_pixel_range_full_image>= y_pmin);

%% CENTER DEFINITION


diameter_ganymede    = 70.373;   % from computation
x_pixel_center       = observation.x_pcenter;    % to be optimised for
y_pixel_center       = observation.y_pcenter;    % to be optimised for
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
ganymede_centred_subimage_reyleights= ganymede_centred_subimage*count2KRayleight;
sigma_matrix_ganymede_centred_subimage_reyleights = sigma_matrix_ganymede_centred_subimage*count2KRayleight;


%% FIT THE MODEL
% assumed ganymede brightness
ganymede_assumed_brightness = 1.3; % KReyleights
zc                = observation.zc;
Rg                = 2643.1e5; % radius ganymede in cm
R_ganymede_pixel  = diameter_ganymede/2;

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

brightness_mask                             = sqrt(X_gird_ganymede_centred_subimage.^2 + Y_gird_ganymede_centred_subimage.^2) < R_ganymede_pixel/2; % only for finding a brightness value take the brightness in half the radius of the ganymede
ganymede_mask_for_ganymede_centred_subimage = sqrt(X_gird_ganymede_centred_subimage.^2 + Y_gird_ganymede_centred_subimage.^2) < R_ganymede_pixel;   % mask covering the ganymede disk in ganymede_centred image
ganymede_mask_for_model_image               = sqrt(X_grid_model.^2 + Y_grid_model.^2) < R_ganymede_pixel;   % this the mask covering the ganymede in the expanded image

mean_brightness_ganymede = mean(mean(ganymede_centred_subimage_reyleights(brightness_mask)));
IPMandGEO                = mean_brightness_ganymede - ganymede_assumed_brightness;
IPMandGEO_counts         = IPMandGEO/count2KRayleight;

% adjust sigma values
sigma_matrix_ganymede_centred_subimage_reyleight   = sqrt(sigma_matrix_ganymede_centred_subimage.^2+IPMandGEO_counts)*count2KRayleight;
mean_brightness_ganymede   = mean_brightness_ganymede - IPMandGEO;
ganymede_centred_subimage_reyleights = ganymede_centred_subimage_reyleights - IPMandGEO;


% obtain polynomial fit of the image
% the x and y axis must be centered with zero at the ganymede center.
% Otherwise the transmission model won't work
[x_fit,y_fit,z_fit,weights_fit] = prepareSurfaceData(X_gird_ganymede_centred_subimage,Y_gird_ganymede_centred_subimage,ganymede_centred_subimage_reyleights,1./sigma_matrix_ganymede_centred_subimage_reyleight.^2);

exclude_set = ~excludedata(x_fit,y_fit,'box',[-box_radial_extension*R_ganymede_pixel ,...
    +box_radial_extension*R_ganymede_pixel ,...
    -box_radial_extension*R_ganymede_pixel ,...
    + box_radial_extension*R_ganymede_pixel]);



% decide max order for background polynomial surface
total_combination = length(n0_range)*length(poly_surface_order_range); 
fprintf("totol number of combination tested : %i\n",total_combination)
proceed_check=input("proceed [y/n]  ","s");

count_process = 1;
if proceed_check =='y'
    h = waitbar(0,'Please wait...');
    for order = poly_surface_order_range 
        
 
        chi2_rad     = zeros(length(n0_range),1);
        chi2_vert    = zeros(length(n0_range),1);
        chi2_hor     = zeros(length(n0_range),1);
    
        poly_order = order;
        syms n0
        syms a [poly_order+1,poly_order+1]
        syms x y
    
        polynomial_surface = 0;
        
        for ii=0:poly_order
            for jj =0:poly_order-ii
                polynomial_surface = polynomial_surface + a(ii+1,jj+1)*x^ii*y^jj;
            end
        end
    
        
    
        % H corona absorbtion model outside ganymede disk
        Nh  = n0*Rg*(R_ganymede_pixel./sqrt(x.^2+y.^2)).*pi;
        tau = zc*Nh;
        T   = exp(-tau);

        % H corona emission model outside ganymede disk
        g_factor  = obseravtion.g_factor;; %[kR/cm2]
        I_H_corona_emission = Nh*g_factor;
    
        % backgroumd surface model
        total_model_sym     = polynomial_surface.*T + I_H_corona_emission  ; % background surface plus the Transmission plus corona emission
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
        %options.Weights    = weights_fit; % still to be decided if we add them
        %or not
    
        counter = 1;

        
        for n0_value=n0_range
    
            model_fit = fit([x_fit,y_fit],z_fit,myfittype,options,"problem",n0_value);
            % create model image from fit
            bkg_only_model_image = double(model_fit(X_grid_model,Y_grid_model));
            
            
            
            % assume constant brightness at disk as constant
            row_model_with_ganymede = bkg_only_model_image;
            row_model_with_ganymede(ganymede_mask_for_model_image)  = mean_brightness_ganymede;
    
            final_model_image_PSF                         = conv2(row_model_with_ganymede,PSF,'same');% apply PSF
            final_model_image_PSF                         = final_model_image_PSF(y_original_mask,x_original_mask);
            rsw_model_image                               = row_model_with_ganymede(y_original_mask,x_original_mask); % image size before the PSF convolution is applied
    
            chi2 = sum(sum((final_model_image_PSF-ganymede_centred_subimage_reyleights).^2./final_model_image_PSF))/numel(final_model_image_PSF);
    
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
    
            chi2_rad(counter)     = sum((model_image_radial_plot         - ganymede_centred_subimage_radial_plot        ).^2./error_bar_radial_plot.^2  )/numel(model_image_radial_plot);
            chi2_vert(counter)    = sum((model_image_vertical_box_plot   - ganymede_centred_subimage_vertical_box_plot  ).^2./error_bar_vertical_box_plot.^2  )/numel(model_image_vertical_box_plot  );
            chi2_hor(counter)     = sum((model_image_horizontal_box_plot - ganymede_centred_subimage_horizontal_box_plot).^2./error_bar_horizontal_box_plot.^2 )/numel(model_image_horizontal_box_plot );
            counter = counter+1;
            count_process = count_process +1;
            waitbar(count_process/total_combination,h)
        
        end
        chi2_radial_result.(['poly',num2str(order)]) = chi2_rad;
        chi2_vert_result.(['poly',num2str(order)])  = chi2_vert;
        chi2_hor_result.(['poly',num2str(order)])   = chi2_hor;
    end

close(h)
    
fileID = fopen(['optimization_',observation.name,'.txt'],'w');
fprintf(fileID,'----------------------------------------\n');
fprintf(fileID,'Column Density Optimization Results\n');
fprintf(fileID,'----------------------------------------\n');
fprintf(fileID,'Observation ID : %s\n',observation.file);
fprintf(fileID,'----------------------------------------\n');
    
    
for order = poly_surface_order_range

    figure
    ax = formal_axes(gca());
    hold on
    plot(n0_range,chi2_hor_result.(['poly',num2str(order)]),DisplayName="horizonal cut")
    plot(n0_range,chi2_vert_result.(['poly',num2str(order)]),DisplayName="vertical cut")
    plot(n0_range,chi2_radial_result.(['poly',num2str(order)]),DisplayName="radial cut")

    [optimal_chi_radial,best_index_radial]         = min(chi2_radial_result.(['poly',num2str(order)]));
    [optimal_chi_horizontal,best_index_horizontal] = min(chi2_hor_result.(['poly',num2str(order)]));
    [optimal_chi_vertical,best_index_vertical]     = min(chi2_vert_result.(['poly',num2str(order)]));


    scatter(n0_range(best_index_horizontal),optimal_chi_horizontal,"b","filled",HandleVisibility="off")
    scatter(n0_range(best_index_vertical),optimal_chi_vertical,"b","filled",HandleVisibility="off")
    scatter(n0_range(best_index_radial),optimal_chi_radial,"b","filled",HandleVisibility="off")

    title(sprintf("polynomial surface order %s",num2str(order)))
    xlabel("colum density $[1/cm^3]$")
    ylabel("$\chi^2$")
    grid on
    formal_legend(legend())


    fprintf(fileID,'Polynomial surface order : %i\n',order);
    formatSpec = '%-10s cut  Chi2 -> %-4.4f | optimal colum density [1/cm3] -> %5.2f\n';
    fprintf(fileID, formatSpec ,  'radial',optimal_chi_radial,n0_range(best_index_radial));
    fprintf(fileID, formatSpec ,  'vertical',optimal_chi_vertical,n0_range(best_index_vertical));
    fprintf(fileID, formatSpec ,  'horizontal',optimal_chi_horizontal,n0_range(best_index_horizontal));
    image_file_name = sprintf(['../images/optimised_chi_measure_poly_%i_',observation.name,'.eps'],order);
    exportgraphics(ax,image_file_name)
end
end

if music_intretainment
 clear sound;
end




    

