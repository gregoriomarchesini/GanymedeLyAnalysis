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
