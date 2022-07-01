function varargout = radial_bins_average(image,col_index_center,row_index_center,inner_radius,outer_radius,nbins,sigma_matrix)

% developer : Gregorio Marchesini
% date : 26/04/2022
%
% description : 
% ------------
% given an image of a planet from a STIS observation and its center, the 
% functions outputs a vector named stat_per_bin of length nbins.
% the element i-th of the vector stat_per_bin is the desired statistics value
% of the pixels in the i-th radial bin produced by creating a circualr
% sector around the planet center with inner radius "inner_radius"
% and outer radius "outer_radius". This circular sector is divided into
% nbins 
%             --------
%            /  |  |  \
%           | \/----\/ |
%           |- |    | -|
%           |- |    | -|
%           |  \----/  |
%           \ /  |   \ /
%            ----------
% 
%     
% parameters : 
% ------------
% image            (array) : image to be analysed
% col_index_center  (int)  : column index of the planet center in the image
% row_index_center  (int)  : row index of the planet center in the image
% inner_radius      (int)  : inner radius of the circular section
% outer_radius      (int)  : outer radius of the circular section
% nbin              (int)  : number of bins in which the sections
%                            should be divided
% stat              (str)  : statistic to be computed ("std"  : standard deviation, "mean" : mean       )
% sigma_matrix(array)      : indicates the standard deviation of each
%                            point in the matrix "image". If the parameter
%                            is given, the error bar of the final line plot
%                            is returned                                                   
% retruns :
% --------
%
% stat_per_bin         (array(1,nbins)) : average pixel value per bin  (or standard deviation)
% center_angles_range  (array(1,nbins)) : bin center angle to use for plotting
% error_bar            (array(1,nbins)) : error value for each bin (only outputted if sigma matrix is given)


arguments
    image             (:,:) double
    col_index_center  double {mustBePositive(col_index_center),mustBeGreaterThan(col_index_center,0)}
    row_index_center  double {mustBePositive(row_index_center ),mustBeGreaterThan(row_index_center,0)}
    inner_radius      double {mustBePositive(inner_radius),mustBeGreaterThan(inner_radius,0)}
    outer_radius      double {mustBePositive(outer_radius),mustBeGreaterThan(outer_radius,0)}
    nbins             double {mustBePositive(nbins),mustBeGreaterThan(nbins,0)}
    sigma_matrix      (:,:) double = zeros(size(image))

end

[rows,cols]    = size(image);
if col_index_center+outer_radius  >cols || row_index_center+outer_radius  >rows
    error("the outer radius is outside the matix dimension")
end

[Xmap,Ymap]    = meshgrid(1:cols,1:rows);
% shift the center 
Xmap    = Xmap - col_index_center;
Ymap    = Ymap - row_index_center;
R_map   = sqrt(Xmap.^2+Ymap.^2);
Phi_map = atan2(Ymap,Xmap);

angles_range           = linspace(-pi,pi,nbins+1); %n equispaced bins have n+1 nodes
angle_step             = angles_range(2) - angles_range(1) ;
center_angles_range    = angles_range(1:end-1); %n equispaced bins have n+1 nodes
center_angles_range    = center_angles_range+angle_step/2; % this tells where is the center of each bin in terms of angle
stat_per_bin           = zeros(1,nbins);

if nargin ==7
   error_per_bin       =zeros(1,nbins);
end

for angle_index = 1:length(angles_range)-1
   
    current_angle = angles_range(angle_index);
    phi_mask      = Phi_map>=current_angle & Phi_map <=(current_angle+angle_step);
    r_mask        = R_map >=inner_radius & R_map <=outer_radius;
   
    stat_per_bin(angle_index) = mean(image(r_mask & phi_mask));  
    
      
    if nargin ==7
       error_per_bin(angle_index) = sqrt(sum(sigma_matrix(r_mask & phi_mask).^2))/numel(sigma_matrix(r_mask & phi_mask));
    end


  
end

if nargin <7
    varargout={stat_per_bin,center_angles_range};
elseif nargin ==7
    varargout={stat_per_bin,center_angles_range,error_per_bin};
end
