function [best_col_index_center, best_row_index_center,epsilon_struct] = Giorno_method(observation_image,model_image,assumed_col_index_center,assumed_row_index_center,radius_inner,radius_outer,number_of_bins)

% developer : Gregorio Marchesini
% date : 26/04/2022
%
% description : 
% ------------
%
% given a moon observation image and a model image, it computes the best
% column and row center location based on the Giorno method. Note that the 
% model image is assumed to be centered correctly at assumed_col_index_center and 
% assumed_row_index_center pixel location.
%
%
%
%
%     
% parameters : 
% ------------
% observation_image        (array) : obseravtion image
% model_image              (array) : model image (same size as observation image)
% assumed_col_index_center int     : assumed initial column index center location
% assumed_row_index_center int     : assumed initial row index centerlocation
% radius_inner             double  : inner binning radius
% radius_outer             double  : outer binning radius
% number_of_bins           int     : number of bins to be used in the average computation                                             
%
% retruns :
% --------
%
% 
% best_col_center_location int    : best col center location based on the
%                                   Giorno centering method
% best_row_center_location int    : best row center location based on the
%                                  Giorno centering method
% epsilon_struct           (struct) : single row structure with field
                                    
                                    % 1) giorno_matrix
                                    % 2) x_shift
                                    % 3) y_shift 
                                    % 4) average_difference_per_bin_at_best_center
                                    % -> average difference for each bin
                                    %    at the best found center location
                                    % 5) bins angles center 
                                    % -> bins angular location
                                    %
                                    %
                                    % The giorno_matrix can be used for
                                    % plotting like
                                    % hold on
                                    % imagesc(epsilon_struct.x_shift,epsilon_struct.y_shift,epsilon_struct.giorno_matrix);
                                    % ax = gca();
                                    % ax.XLim = [epsilon_struct.x_shift(1),epsilon_struct.x_shift(end)];
                                    % ax.YLim = [epsilon_struct.y_shift(1),epsilon_struct.y_shift(end)];
                                    % ax.XLabel.String=["offset from assumed center [pixel]"];
                                    % ax.YLabel.String=["offset from assumed center [pixel]"];
                                    % h  = colorbar;
                                    % h.YLabel.String = "$\epsilon$";

arguments
    observation_image        (:,:) double 
    model_image              (:,:) double
    assumed_col_index_center double {mustBeInteger}
    assumed_row_index_center double {mustBeInteger}
    radius_inner             double  
    radius_outer             double  
    number_of_bins           double
end

if any(size(observation_image)~=size(model_image))
    error("model image and real image must have the same size")
end

max_pixel_shift  = 5 ;
x_shift          = -(max_pixel_shift) : max_pixel_shift;
y_shift          = x_shift; 
epsilon_matrix   = zeros(length(x_shift));

mask_ber_bin_model = radial_bins_mask(model_image,assumed_col_index_center,assumed_row_index_center,radius_inner,radius_outer,number_of_bins);
for jj = y_shift 
       for kk = x_shift
           average_difference_per_bin = zeros(1,number_of_bins);
           mask_ber_bin_observation = radial_bins_mask(observation_image,assumed_col_index_center+kk,assumed_row_index_center+jj,radius_inner,radius_outer,number_of_bins);
           
           for bin_number = 1:number_of_bins

               pixel_per_bin_observation  = observation_image(mask_ber_bin_observation(bin_number).mask);
               pixel_per_bin_model        = model_image(mask_ber_bin_model(bin_number).mask);
               average_difference_per_bin(bin_number) =  mean(pixel_per_bin_observation - pixel_per_bin_model);
           
           end
           epsilon_matrix(jj==y_shift,kk==x_shift) = std(average_difference_per_bin);
       end
end

[best_row_index_shift,best_col_index_shift] = find(epsilon_matrix==min(min( epsilon_matrix)));

epsilon_struct.giorno_matrix = epsilon_matrix;
epsilon_struct.x_shift = x_shift;
epsilon_struct.y_shift = y_shift;

if length(best_row_index_shift) ~=1 
    warning("multiple optimal centers found appling Giorno Method. Inspect epsilon matrix");
    best_row_index_shift=best_row_index_shift(1);best_col_index_shift=best_col_index_shift(1);
end

best_row_index_center = assumed_row_index_center +y_shift(best_row_index_shift);
best_col_index_center = assumed_col_index_center +x_shift(best_col_index_shift);

% final best radial average value
average_difference_per_bin = zeros(1,number_of_bins);
mask_ber_bin_observation   = radial_bins_mask(observation_image,best_col_index_center,best_row_index_center,radius_inner,radius_outer,number_of_bins);     
for bin_number = 1:number_of_bins
   pixel_per_bin_observation  = observation_image(mask_ber_bin_observation(bin_number).mask);
   pixel_per_bin_model        = model_image(mask_ber_bin_model(bin_number).mask);
   average_difference_per_bin(bin_number) =  mean(pixel_per_bin_observation - pixel_per_bin_model);
end

epsilon_struct.average_difference_per_bin_at_best_center = average_difference_per_bin;
epsilon_struct.bins_angles_center = [mask_ber_bin_observation(:).angle];


