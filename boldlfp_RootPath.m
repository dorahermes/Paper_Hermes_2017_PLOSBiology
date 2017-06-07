function rootPath = boldlfp_RootPath()
% Return the path to the root boldlfp  directory
%
% This function must reside in the directory at the base of the boldlfp
% directory structure.  It is used to determine the location of various
% sub-directories.

rootPath=which('boldlfp_RootPath');

rootPath=fileparts(rootPath);

return
