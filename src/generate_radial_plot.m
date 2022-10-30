function varargout = generate_radial_plot(image,col_index_center,row_index_center,max_pixel_radius,ring_thikness,sigma_matrix)
% developer : Gregorio Marchesini
    % date : 26/04/2022
    %
    % description : 
    % ------------
    % given an image of a planet from a stis observation, the radial analus
    % average plot is produced. Starting from the planet center, a analus
    % of thikness ring_thikness is created and the average brightness is
    % calculated. This is repeated (max_pixel_radius/ring_thikness) times
    %
    % parameters : 
    % ------------
    % image (array)          : image to be analysed
    % col_index_center(int)  : column index of the planet center
    % row_index_center(int)  : row index of the planet center
    % max_pixel_radius(float): radial extension of the plot
    % ring_thikness(float)   : thikness of each anulus
    % sigma_matrix(array)    : indicates the standard deviation of each
    %                          point in the matrix "image". If the parameter
    %                          is give, the error bar of the final line plot
    %                          is returned
    %
    % retruns :
    % --------
    % (*) optional 
    %
    % radial_plot          (array(1,:)) : average brightness computed along
    %                                     the sequence of anulus
    % radial_range         (array(1,:)) : center radius of each annulus
    %                                     (r_min + r_max)/2
    % error_bar_values (*) (array(1,:)) : error bar derived by the given
    %                                     standard deviation matrix (only if weight matrix is given)
    
     arguments
        image            (:,:) double
        col_index_center double
        row_index_center double
        max_pixel_radius double
        ring_thikness    double
        sigma_matrix      (:,:) double = false
     end
    
     [rows,cols] = size(image);
    
    if nargin==6
        [rows_w,cols_w] = size(sigma_matrix);
        if rows ~= rows_w || cols~=cols_w
            error("sigma matrix is not same size as image matrix")
        end
    end
    

[rows,cols] = size(image);
[Xmap,Ymap] = meshgrid(1:cols,1:rows);
Xmap        = Xmap - col_index_center;
Ymap        = Ymap - row_index_center;
R_map       = sqrt(Xmap.^2 + Ymap.^2);

radius_range      = 0:ring_thikness:max_pixel_radius;
mid_points_range  = mid_points(radius_range);
image_radial_mean = zeros(1,length(radius_range)-1);
error_bar_width   = zeros(1,length(radius_range)-1);


for jj = 1:(length(radius_range)-1)
    
    r_mask = R_map>=radius_range(jj) & R_map<=(radius_range(jj+1));
    image_radial_mean(jj) = mean(image(r_mask));

    if nargin==6
      error_bar_width(jj)    = sqrt(sum(sigma_matrix(r_mask).^2))/numel(sigma_matrix(r_mask));
    end
end

 varargout{1} = image_radial_mean;
 varargout{2} = mid_points_range;
if nargin==6
  varargout{3}  = error_bar_width;
end

end