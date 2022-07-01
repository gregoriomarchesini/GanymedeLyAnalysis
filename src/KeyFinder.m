function [result,description]= KeyFinder(fits_info_obj,desired_key,varargin)

% Developer: Gregorio Marchesini 
% Date: 1 March 2021


% This function was create in order to easily obtain Keyword values from
% the Primary Header of a fits file. Using arrayfun it is possible to  give multiple
% words input to this function if necessary. Otherwise you can simply
% insert this ina loop function and obtain all the desired words

%  INPUTS
%  ––––––
% 
% 
%  Variable           Description                   Data Type
% ––––––––––          ––––––––––––                  ––––––––––
%
% fits_info_obj     object obtained reading            fits-object or string
%                   a fits file with fitsinfo
%                   OR you can put directly the
%                   (complete) name of the observation
%
%
% desired_key       Keyword to reseacrh in the          string
%                   Primary eader
%
% extension         (optional) parameter if you         string
%                   want a keyword that is not
%                   in the 'PrimaryHeader'.
%                   examples are (Binary,Image)
%                   if the field dosen't exist
%                   you will receive an error
%   
% print             (optional) prints all the          (str)
%                   available keywords
%
%-------------------------------------------------------------------------
%
%  OUTPUTS
%  –––––––
%
%  Variable           Description                     Data Type
% ––––––––––          ––––––––––––                    ––––––––––
%
%  result             String related to Keyword       string
%                     selected
% 
% description        description related to Keyword   string
%                    selected
%
% NOTE: A complete list of the possible headers can be found in the HST
% Doucumentation. The word must be inserted correctly.
% also note that the optional values must be inserted in the correct order
%-------------------------------------------------------------------------%


if ~isstruct(fits_info_obj)
    fits_info_obj = fitsinfo(fits_info_obj);
end


% extract cell array of all the the Keywords


if isempty(varargin) 
        table_key    = fits_info_obj.PrimaryData.Keywords; 

elseif strcmpi(varargin{1},'print')
        table_key    = fits_info_obj.PrimaryData.Keywords; 
 
    
elseif ~strcmpi(varargin{1},'print')
    
    extension   = varargin{1};
    4;
    try
     table_key   = fits_info_obj.(extension).Keywords; 
    catch
        
        try
           extension(1) = upper(extension(1));
           table_key    = fits_info_obj.(extension).Keywords; 
        catch
           error('The extension %s is not available',extension)
        end
        
     end
end

[r,~]        = size(table_key);
flag         = 0;

for i=1:r
    
    current_key=table_key{i,1};
    
    if any(strcmp(varargin,'print'))
       fprintf([current_key,'\n'])
    end
    
    if strcmp(desired_key,current_key)
        flag=1;
        result      = table_key{i,2};
        description = table_key{i,3};
    
    end  
end

if flag==0
    result      = 'Not found';
    description = 'Not found';
end

end
