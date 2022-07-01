function C = circular_kernel(d1,d2,magnitude,outer_magnitude)

% Developer: Gregorio Marchesini 
% Date: 1 March 2021
% Description
% –––––––––––
% This function creates a matrix with an anular section of pixels.
%
%
%  INPUTS 
%  –––––––
% 
%  Variable           Description                         Data Type
%  ––––––––           –––––––––––                         –––––––––
%
%     d1              diameter of the outer circle         scalar
%                     in pixels
%
%     d2              diameter of the inner circle         scalar
%                     in pixels
%
%
% ------------------------------------------------------------------------
%
%  OUTPUTS
%  –––––––
% 
%  Variable           Description                   Data Type
%  ––––––––           –––––––––––                   –––––––––
%
%   M                 matrix with radious           array  (rxr)
%                     containing a circular 
%                     kernel
%
% ------------------------------------------------------------------------
%
%  Options
%  –––––––
%
%   Name              Description                   Data Type
%  ––––––––           –––––––––––                   –––––––––
%  
%  'magnitude'        define the value to           scalar
%                     be given at the rim 
%                     of the anulus
%
%  'outer_magnitude'  define the value to           scalar
%                     be given at the center 
%                     of the anulus
% ------------------------------------------------------------------------
%  NOTE: d must be integer 
%

%% Parse the function

arguments
    d1              (1,1) double 
    d2              (1,1) double   {mustBeLessThan(d2,d1)}
    magnitude       (1,1) double  = 1
    outer_magnitude (1,1) double  = 0

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function

d1 = round(d1);
d2 = round(d2);
C    = zeros(d1);
[X,Y]= meshgrid(1:d1);                       
    % center of the matrix;
    
% case of even number radious
if  ~rem(d1,2)
     
    % the central pixel is half way
    x_center=d1/2 + 0.5;
    y_center=d1/2 + 0.5;

% case of odd number radious
else
   
    x_center=floor(d1/2)+1;
    y_center=floor(d1/2)+1;
end


X    = X-x_center;
Y    = Y-y_center;
rpow2  = (X.^2+Y.^2);


C(rpow2<=(d1/2)^2 & rpow2>=(d2/2)^2)   = magnitude;
C(rpow2<(d2/2)^2)                      = outer_magnitude;

end

