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
% This function estimates the mean and standard deviation for every grid
% point in the Allen Mouse Common Coordinate Framework (CCF v3) (100 micron
% resolution), separately for coronal and sagital data, and saves them to
% storage/expressionMatrixMeans_brain.mat. The computation can take several
% hours for randAmount>300. "no data" values are excluded for mean/std
% calculation.
%
% params:
%          randAmount  Amount of genes that will be used to to compute the
%                      mean/std
%
function [] = estimate_Mean_STD_of_GeneExpression(randAmount)

    mkdir('storage')
    mkdir('storage/expressions')

    %download a list of all entrez-ids of all genes available from the
    %Allen brain atlas
    disp('Download all available Entrez-IDs from brain-map.org')
    urlwrite('http://api.brain-map.org/api/v2/data/query.xml?num_rows=30000&only=genes.entrez_id&criteria=model::Gene,rma::criteria,products[abbreviation$eq%27Mouse%27]','storage/allEntrezInDB.xml');
    root = xmlread('storage/allEntrezInDB.xml');

    elems = root.getElementsByTagName('entrez-id');
    entrezAll = {};

    for i = 0:(elems.getLength-1)
        entrezAll{i+1} =char(elems.item(i).getTextContent);
    end

    %generate an index for the entrez ids, that can be randomized to
    %ensure, that random genes are picked from the database
    entrezIndizes=[];
    for i = 1:length(entrezAll)
        if(length(entrezAll{i})>0)
            entrezIndizes=[entrezIndizes,i];
        end
    end

    disp(sprintf('%d Entrez-IDs downloaded!',length(entrezIndizes)))
    entrezIndizes=entrezIndizes(randperm(length(entrezIndizes))); %randomization of the index
    actEntrezIndizes=0;

    disp('Download 100 micron annotations from brain-map.org')

    %download the annotation of the Allen Mouse Common Coordinate Framework
    %in this case, only because we need the dimension-size of the space
    urlwrite('http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2015/annotation_100.nrrd','storage/annotation_100.nrrd');
    [newANOGD, metaDMASK] = nrrdread('storage/annotation_100.nrrd');
    ANOGD=permute(newANOGD,[2 1 3]);

    disp('Complete!')

    expressionIndexSag=[];
    expressionMatrixAllSag =sparse(size(ANOGD,1)*size(ANOGD,2)*size(ANOGD,3),randAmount);

    expressionIndexCor=[];
    expressionMatrixAllCor =sparse(size(ANOGD,1)*size(ANOGD,2)*size(ANOGD,3),randAmount);

    actLoaded=0;
    disp('Download random gene expression data: ')
    while(size(expressionIndexSag,2)<randAmount || size(expressionIndexCor,2)<randAmount)
        actLoaded=actLoaded+1;
        actEntrezIndizes=actEntrezIndizes+1;
        entrez=entrezAll{entrezIndizes(actEntrezIndizes)};

        %check if genes with a certain entrez id can be found in the
        %Allen Brain Atlas API
        urlwrite(strcat('http://api.brain-map.org/api/v2/data/query.xml?criteria=model::SectionDataSet,rma::criteria,[failed$eq%27false%27],products[abbreviation$eq%27Mouse%27],genes[entrez_id$eq%27',entrez,'%27]'), 'storage/expression.xml');
        root = xmlread('storage/expression.xml');

        if str2double(root.item(0).getAttribute('num_rows'))<=0
            disp(sprintf('Could not find %s',entrez))
        else

            ids=[];
            referenceSpace=[];
            try
                for j=1:2:100
                    selectionDataSet = root.item(0).item(0).getChildNodes.item(j);
                    [ attribute, childnodes ] = goToChildNode( selectionDataSet,'id',0 );
                    id = str2num(attribute);
                    ids=[ids,id];
                    [ attribute, childnodes ] = goToChildNode( selectionDataSet,'reference-space-id',0 );
                    rid = str2num(attribute);
                    referenceSpace=[referenceSpace,rid];
                end
            catch err
            end

            expressionMatrixCor =sparse(size(ANOGD,1)*size(ANOGD,2)*size(ANOGD,3),1);
            expressionMatrixSag =sparse(size(ANOGD,1)*size(ANOGD,2)*size(ANOGD,3),1);
            amountCor=sparse(size(ANOGD,1)*size(ANOGD,2)*size(ANOGD,3),1);
            amountSag=sparse(size(ANOGD,1)*size(ANOGD,2)*size(ANOGD,3),1);


            %Download the data from the Allen Brain Atlas API, or take
            %the data from the local hard-drive, if the gene-expression
            %data has already been stored to "storage/expressions"

            %if a gene has more than one sagital or more than one
            %coronal image, the average will be computed
            for act=1:size(ids,2)
                id=ids(act);
                rid=referenceSpace(act);

                if( (rid==9 && size(expressionIndexCor,2)<randAmount)||(rid==10 && size(expressionIndexSag,2)<randAmount))
                    try
                        if(not(exist(sprintf('storage/expressions/%d/density.raw',id), 'file')))
                            if(rid==9)
                                disp(sprintf('      Download %d. from brain-map.org. #coronal: %d/%d',actLoaded,size(expressionIndexCor,2)+1,randAmount))
                            else
                                disp(sprintf('      Download %d. from brain-map.org. #sagital: %d/%d',actLoaded,size(expressionIndexSag,2)+1,randAmount))
                            end
                            urlwrite(sprintf('http://api.brain-map.org/grid_data/download/%d?include=density',id), 'storage/temp.zip');
                            unzip('storage/temp.zip',sprintf('storage/expressions/%d',id));
                        else
                            if(rid==9)
                                disp(sprintf('      Load %d. from local. #coronal: %d/%d',actLoaded,size(expressionIndexCor,2)+1,randAmount))
                            else
                                disp(sprintf('      Load %d. from local. #sagital: %d/%d',actLoaded,size(expressionIndexSag,2)+1,randAmount))
                            end
                        end

                        geneGridSize = [67 41 58];
                        fid = fopen(sprintf('storage/expressions/%d/density.raw',id), 'r', 'l'  );
                        expression = fread( fid, prod(geneGridSize), 'float' );
                        fclose(fid);
                        expression = reshape( expression, geneGridSize );

                        %Interpolate 200 micron gene expression data to 100 micron resolution
                        expression_100 = interpolateFrom200micron(expression);

                        %mirror
                        expression_100_mir=expression_100(1:end,1:end,(end):-1:1);
                        expression_100(1:end,1:end,(floor(size(expression_100,3)/2)+1):end)=expression_100_mir(1:end,1:end,(floor(size(expression_100_mir,3)/2)+1):end);
                        exp = reshape(expression_100,size(ANOGD,1)*size(ANOGD,2)*size(ANOGD,3),1);

                        expBiggerZero=exp;

                        expBiggerZero(expBiggerZero<0)=0;


                        if(rid==9)
                            expressionMatrixCor=(expressionMatrixCor+expBiggerZero);
                            amountCor=amountCor+(exp>=0);
                        end
                        if((rid==10))
                            expressionMatrixSag=(expressionMatrixSag+expBiggerZero);
                            amountSag=amountSag+(exp>=0);
                        end
                    catch err
                        disp(sprintf('            Could not open %s',id))
                    end
                end
            end

            %this is the average computation if multiplie sagital/coronal
            %volumes are available. It excludes "no data" voxels
            if(sum(expressionMatrixCor(:)>0) && sum(expressionMatrixSag(:)>0))
                expressionMatrixCor(isnan(expressionMatrixCor))=0;
                expressionMatrixSag(isnan(expressionMatrixSag))=0;
                if(size(expressionIndexSag,2)<randAmount)
                    expressionMatrixAllSag(:,size(expressionIndexSag,2)+1)=(expressionMatrixSag./amountSag);
                    expressionIndexSag= [expressionIndexSag,size(expressionIndexSag,2)+1];
                end
                if(size(expressionIndexCor,2)<randAmount)
                    expressionMatrixAllCor(:,size(expressionIndexCor,2)+1)=(expressionMatrixCor./amountCor);
                    expressionIndexCor= [expressionIndexCor,size(expressionIndexCor,2)+1];
                end
            else if(sum(expressionMatrixCor(:)>0) && sum(expressionMatrixSag(:)==0))

                    if(size(expressionIndexCor,2)<randAmount)
                        expressionMatrixCor(isnan(expressionMatrixCor))=0;
                        expressionMatrixAllCor(:,size(expressionIndexCor,2)+1)=(expressionMatrixCor./amountCor);
                        expressionIndexCor= [expressionIndexCor,size(expressionIndexCor,2)+1];
                    end
                else if(sum(expressionMatrixCor(:)==0) && sum(expressionMatrixSag(:)>0))
                        if(size(expressionIndexSag,2)<randAmount)
                            expressionMatrixSag(isnan(expressionMatrixSag))=0;
                            expressionMatrixAllSag(:,size(expressionIndexSag,2)+1)=(expressionMatrixSag./amountSag);
                            expressionIndexSag= [expressionIndexSag,size(expressionIndexSag,2)+1];
                        end
                    end
                end
            end
        end
    end


    disp('Calculate mean and std for every gridpoint (for sagital experiments)... ')
    disp('WARNING: This can take up to several hourse without notable output')
    expressionMatrixCorMean=nanmean(expressionMatrixAllCor,2);
    %disp(sum(isnan(expressionMatrixCorMean)))
    expressionMatrixSagMean=nanmean(expressionMatrixAllSag,2);
    %disp(sum(isnan(expressionMatrixSagMean)))

    disp('Calculate mean and std for every gridpoint (for coronal experiments)... ')
    disp('WARNING: This can take up to several hourse without notable output')
    expressionMatrixCorSTD=nanstd(expressionMatrixAllCor,1,2);
    %disp(sum(isnan(expressionMatrixCorSTD)))
    expressionMatrixSagSTD=nanstd(expressionMatrixAllSag,1,2);
    %disp(sum(isnan(expressionMatrixSagSTD)))

    save(sprintf('storage/expressionMatrixMeans_brain.mat'),'expressionMatrixCorMean','');
    save(sprintf('storage/expressionMatrixMeans_brain.mat'),'expressionMatrixSagMean','-append');
    save(sprintf('storage/expressionMatrixMeans_brain.mat'),'expressionMatrixCorSTD','-append');
    save(sprintf('storage/expressionMatrixMeans_brain.mat'),'expressionMatrixSagSTD','-append');
    %clear('expressionMatrixAllSag');
    %clear('expressionMatrixAllCor');

    disp('Computation complete!')
end
