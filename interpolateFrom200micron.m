% Copyright (C) 2016 VRVis.
% All rights reserved.
% Contact: VRVis Forschungs-GmbH (office@vrvis.at)
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 3. All advertising materials mentioning features or use of this software
%    must display the following acknowledgement:
%    This product includes software developed by the VRVis Forschungs-GmbH.
% 4. Neither the name of the VRVis Forschungs-GmbH nor the
%    names of its contributors may be used to endorse or promote products
%    derived from this software without specific prior written permission.
%
%interpolates 200 micron expression data to 100 micron expression data
%The code was inspired by "ACCESSING ANNOTATION AND GRIDDED PROJECTION DATA MAPPED TO CCF V2", MAY 2015 v.1, Page 6:
%http://download.alleninstitute.org/informatics-archive/october-2014/mouse_projection/Accessing_October_2014_projection_data%20.pdf
function [ expression_100 ] = interpolateFrom200micron( expression )
    [xi,yi,zi] = meshgrid(1:0.5:41,1:0.5:67,1:0.5:58); %note: matlab transposes x-y
    d = expression; d(d<0) = 0; % fill in missing data as zeroes
    expression_100 = interp3(d ,xi,yi,zi,'linear');

    % Handle "no data" (-1) voxels.
    % Create a mask of "data" vs "no data" voxels and apply linear interpolation
    m = zeros(size(expression));
    m(expression  >= 0) = 1; mi = interp3(m,xi,yi,zi,'linear');

    % Normalize data by dividing by interpolated mask. Assign value of "-1" to "no data" voxels.
    expression_100 = expression_100 ./ mi;
    expression_100( mi <= 0 ) = -1;

    expression_100=expression_100(1:(end-1),1:(end-1),1:(end-1));
end

