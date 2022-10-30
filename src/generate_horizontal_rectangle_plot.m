function varargout = generate_horizontal_rectangle_plot(image,col_index_center,row_index_center,box_width,sigma_matrix)
    % developer : Gregorio Marchesini
    % date : 26/04/2022
    %
    % description : 
    % ------------
    % given an image of a planet from a STIS observation, the image is 
    % cut horizontally at the row "row_index_center" with a width eqaul to
    % "box_width".
    % After that, the average brightness is taken along the columns of the
    % resulting horizontal matrix
    %
    % parameters : 
    % ------------
    % image (array)         : image to be analysed
    % col_index_center(int) : column index where the vertical matrix is
    %                         centered
    % row_index_center(int) : column index where the vertical matrix is
    %                         centered
    % box_widt(float)       : width of the vertical matrix to be cut
    % sigma_matrix(array)   : indicates the standard deviation of each
    %                         point in the matrix "image". If the parameter
    %                         is give, the error bar of the final line plot
    %                         is returned
    %
    % retruns :
    % --------
    %
    % line_plot        (array(1,n)) : average brightness along the
    %                                 horizontal  matrix width   
    % yrange           (array(1,n)) : x axis of the line_plot. Each
    %                                 yrange(i) corresponds to row index of the image 
    %                                 along which line_plot(i)
    %                                 is computed. the indices are centered
    %                                 at row_index_center so that xrange
    %                                 will be in the range [row_index_cente-box_width/2,row_index_cente+box_width/2]
    % error_bar_values (array(1,n)) : error bar derived by the given
    %                                 standard deviation matrix (only if weight matrix is given)
    
    
     arguments
        image             (:,:) double
        col_index_center  double
        row_index_center  double
        box_width         double
        sigma_matrix      (:,:) double = false
     end
    
    [rows,cols] = size(image);
    
    if nargin==5
        [rows_w,cols_w] = size(sigma_matrix);
        if rows ~= rows_w || cols~=cols_w
            error("sigma matrix is not same size as image matrix")
        end
    end
    
    upper_limit = row_index_center + floor(box_width/2);
    lower_limit = row_index_center - floor(box_width/2);
    box         = image(lower_limit:upper_limit,1:end);
    yrange      = (1:cols)-col_index_center;

    [N,~] = size(box);

    line_plot        = mean(box).';
    varargout{1}     = line_plot ;
    varargout{2}     = yrange;
    % define error values
    if nargin==5
        error_bar_values = sqrt(sum(sigma_matrix(lower_limit:upper_limit,1:end).^2)).'/N;
        varargout{3} = error_bar_values;
    end
end