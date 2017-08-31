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
%Before you can run this code, download the NRRD Format File Reader from 
%http://www.mathworks.com/matlabcentral/fileexchange/34653-nrrd-format-file-reader/content/nrrdread.m
%and copy nrrdread.m into this folder

%#######################################################################
%This part only needs to be done in the beginning
%#######################################################################

% Estimating the mean and standard deviation for every grid point in the 
% Allen Mouse Common Coordinate Framework (CCF v3) (100 micron resolution)
% This is necessary since coronal and sagital data has different
% distributions for every grid point, and therefore need to be normalized
% separately.
% THIS CAN TAKE UP TO SEVERAL HOURS SINCE THE MEAN/STD COMPUTATION IS NOT
% THE FASTEST IMPLEMENTATION!
estimate_Mean_STD_of_GeneExpression(300)

% Downloads a certain amount of random selected genes, and stores them into
% a certain amount of files. There shouldn't be more than 100 genes per
% file, since R has a problem with reading large matlab files. 10 files 
% with 100 genes (=1000 random genes total) are sufficient.  Every grid 
% point will be standardized with the mean and standard deviation computed 
% by the estimate_Mean_STD_of_GeneExpression function. Results will be
% saved into "storage/random_genes"
get_gene_expression_for_random_genes(100,10)

%
% Downloads the injection based connectivity data from the
% Allen Brain Atlas API. The injections are mirrored to both hemispheres. The connectivity of
% every injection image is normalized by its injection volume.
%
get_connectivity_data()


%#######################################################################
%This part only needs to be done when one wants to compute a new dataset
%#######################################################################


% Downloads gene expression for the given testsets and stores them into
% separate files. A .csv file that contains gene expression sest needs to
% be in the following format (without spaces or other special characters!)
%                 name_of_set_1;
%                 genename_1;entrez_id_1;ratio(optional)
%                 genename_2;entrez_id_2;ratio(optional) 
%                 ;;
%                 name_of_set_2;
%                 genename_3;entrez_id_3;ratio(optional)
%                 genename_2;entrez_id_2;ratio(optional)
%                 genename_4;entrez_id_4;ratio(optional)
%                 ;;
get_gene_expression_for_genesets('test_genesets.csv')