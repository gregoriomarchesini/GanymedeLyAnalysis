function varargout = generate_vertical_rectangle_plot(image,col_index_center,box_width,sigma_matrix)
    % developer : Gregorio Marchesini
    % date : 26/04/2022
    %
    % description : 
    % ------------
    % given an image of a planet from a stis observation, the image is 
    % cut vertically at the column "x_index_center" with a width eqaul t
    % "box_width".
    % After that, the average brightness is taken along the rows of the
    % resulting vertical matrix
    %
    % parameters : 
    % ------------
    % image (array)         : image to be analysed
    % col_index_center(int) : column index where the vertical matrix is
    %                         centered
    % box_width(float)       : width of the vertical matrix to be cut
    % sigma_matrix(array)   : indicates the standard deviation of each
    %                         point in the matrix "image". If the parameter
    %                         is give, the error bar of the final line plot
    %                         is returned
    %
    % returns :
    % --------
    %
    % line_plot        (array(1,n)) : average brightness along the verticla
    %                                  matrix width
    % error_bar_values (array(1,n)) : error bar derived by the given
    %                                 standard deviation matrix (only if weight matrix is given)
    
    
     arguments
        image            (:,:) double
        col_index_center double
        box_width        double
        sigma_matrix     (:,:) double = false
     end
    
    [rows,cols] = size(image);
    
    if nargin==4 
        [rows_w,cols_w] = size(sigma_matrix);
        if rows ~= rows_w || cols~=cols_w
            error("sigma matrix is not same size as image matrix")
        end
    end
    
    upper_limit = col_index_center + floor(box_width/2);
    lower_limit = col_index_center - floor(box_width/2);
    box = image(1:end,lower_limit:upper_limit);
    [~,N] = size(box);

    line_plot        = mean(box,2);
    varargout{1} = line_plot;

    % define error values
    if nargin==4 
       error_bar_values = sqrt(sum(sigma_matrix(1:end,lower_limit:upper_limit).^2,2))/N;
       varargout{2}      = error_bar_values;
       
    end
end