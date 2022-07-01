classdef FitsImageObject < handle
    
    properties
        image        (:,:) double {mustBeReal, mustBeFinite}
        fitsinfo     (:,:) struct
        center       (1,2) double {mustBeReal, mustBeFinite}
        subimage_box (1,4) double {mustBeReal, mustBeFinite}
        weights  (:,:) double {mustBeReal, mustBeFinite}
    end
    
    methods 
        function [] = read_image(obj,filepath)
            
            [fitsdata,infofits]=fits_flt_read(filepath);
            obj.fitsinfo = infofits;
            obj.image    = fitsdata.IMAGE ;
            obj.weights  = fitsdata.WEIGHTED_IMAGE;
        end
        
        
        function [resized_image,resized_weights] = resize_image(obj,xmin,xmax,ymin,ymax)
            
            resized_image   = obj.image(ymin:ymax,xmin:xmax);
            resized_weights = obj.weights(ymin:ymax,xmin:xmax);
            obj.subimage_box = [xmin,xmax,ymin,ymax];
            
        end
        
        function [result,description]= find_key(obj,desired_key,options)
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
            %  Variable           Description                      Data Type
            % ––––––––––          ––––––––––––                     ––––––––––
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
            %
            %  Name/Value
            %  ––––––----
            %
            %  extension       extenion to look into                string
            %                  [Image,BinaryTable]
            %                  dafault : PrimaryHeader
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
            
            
            arguments
                
                obj                     FitsImageObject
                
                desired_key             (1,:) string
                
                options.extension (1,:) string = "PrimaryData"
                
            end
            
            % extract cell array of all the the Keywords
            
            table_key    = obj.fitsinfo.(options.extension).Keywords;
            [r,~]        = size(table_key);
            flag         = 0;
            
            for i=1:r
                
                current_key=table_key{i,1};
                
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
        
        
        function [x_new,y_new] = obtain_axis_conversions(obj)
            % obtain the axis conversion for plots
            % example - from pixel to spatial position
            %         - from pixel to wavelength
            
            [crpx1 ,~]  =obj.find_key("CRPIX1","extension","Image");
            [crpx2,~]   =obj.find_key("CRPIX2","extension","Image");
            [crval1,~]  =obj.find_key("CRVAL1","extension","Image");
            [cd22  ,~]  =obj.find_key("CD2_2","extension","Image");
            [cd11  ,~]  =obj.find_key("CD1_1","extension","Image");
            
            [rows,cols]  = size(obj.image);
            
            x = 1:cols;
            y = 1:rows;
            
            x_new = crval1 + (x-crpx1)*cd11;
            y_new = (y-crpx2)*cd22;
            
        end
        
        
        
        
        
    end
end





function [data_structure,info_fits]=fits_flt_read(fits_filename)
% author : gregorio marchesini
% date   : 02/2022
% mail   : gregorio.marchesini@gmail.com/gremar@kth.se
%
%  DESCRIPTION
%  -----------
%  reads a fits file with extension _flt.fits. Flat file images are
%  read according to the specifications of the fits files format.
%  In case of flat field data the three data matrices available are
%  a scientific image, a weighted image and a Good bit matrix
%  All this matrices are given in the output data_structure,
%  so that the output of the function has the following fields
%
%  data_structure.IMAGE              % scientific image
%  data_structure.WEIGHTED_IMAGE     % weighted image
%  data_structure.DQ                 %good bits analysis
%
%
%  refernce :
%  https://hst-docs.stsci.edu/acsdhb/chapter-2-acs-data-structure/2-2-acs-file-structure
%
%
%  INPUT
%  —————
%
%                             Description
%                             ———————————
%
%  fits_filename(string)      full/relative path of the observation
%                             file.
%
%                             note : the function works only with
%                             _flt.fits files
%
%
%  OUTPUT                      Description
%  —————                       ____________
%
%
%
% data_structure(strct)       structure constaining data
%                             described in :
%                             https://hst-docs.stsci.edu/acsdhb/chapter-2-acs-data-structure/2-2-acs-file-structure
%
%
%
% info_fits(struct)           info fits structure
%
%
%  example
% ---------
%
%  data_structure(strct)  = fits_flt_read(xxxxxx_flt.fits,1)
%
%  data_structure.IMAGE              % scientific image
%  data_structure.WEIGHTED_IMAGE     % weighted image
%  data_structure.DQ                 % good bits analysis
% ------------------------------------------------------------------------



info_fits = fitsinfo(fits_filename);

data_structure.IMAGE          = fitsread(fits_filename,'image',1);
data_structure.WEIGHTED_IMAGE = fitsread(fits_filename,'image',2);
data_structure.DQ             = fitsread(fits_filename,'image',3);


end

