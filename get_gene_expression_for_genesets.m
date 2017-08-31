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
% This function downloads gene expression data for genesets, and stores 
% them into separate files. There shouldn't be more than 100 genes per
% set, since R has a problem with reading large matlab files. Every grid 
% point will be standardized with the mean and standard deviation computed 
% by the estimate_Mean_STD_of_GeneExpression function. Results will be
% saved into "storage/"
%
% params:
%          inputFile       A .csv file that contains gene expression set in
%                          the following format (without spaces or other
%                          special characters!)
%                          name_of_set_1;
%                          genename_1;entrez_id_1;ratio(optional)
%                          genename_2;entrez_id_2;ratio(optional) 
%                          ;;
%                          name_of_set_2;
%                          genename_3;entrez_id_3;ratio(optional)
%                          genename_2;entrez_id_2;ratio(optional)
%                          genename_4;entrez_id_4;ratio(optional)
%                          ;;
%          
%
function [] = get_gene_expression_for_genesets(inputFile)

    mkdir('storage')
    mkdir('storage/expressions')
    
    load('storage/expressionMatrixMeans_brain.mat')

    %download the annotation of the Allen Mouse Common Coordinate Framework
    %in this case, only because we need the dimension-size of the space
    urlwrite('http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2015/annotation_100.nrrd','storage/annotation_100.nrrd');
    [newANOGD, metaDMASK] = nrrdread('storage/annotation_100.nrrd');
    ANOGD=permute(newANOGD,[2 1 3]);

    outputDirectory = strcat('storage//',inputFile(1:(length(inputFile)-4)));

    mkdir(outputDirectory);

    fid = fopen(inputFile);
    out = textscan(fid,'%s%f%f','delimiter',';');
    fclose(fid);
    genenamesAll=out{1};
    entrezAll=out{2};
    ratiosAll=out{3};

    actFilePos=0;
    setnames=[];

    %traverses the input file and retrieve gene expression data for every
    %geneset in the input file separately
    while(actFilePos<size(entrezAll,1))
        entrez=[];
        genenames=[];
        setname='';
        ratios=[];

        if actFilePos==0
            actFilePos=actFilePos+1;
            setname=genenamesAll{actFilePos};
        end
        if(actFilePos<size(entrezAll,1))
            if(strcmp(genenamesAll{actFilePos},'')==1)
                actFilePos=actFilePos+1;
                setname=genenamesAll{actFilePos};
                actFilePos=actFilePos+1;
            end


            if(actFilePos<size(entrezAll,1))
                while strcmp(genenamesAll{actFilePos},'')==0 && actFilePos<=size(entrezAll,1)
                    if(entrezAll(actFilePos)>0)
                        entrez=[entrez,entrezAll(actFilePos)];
                        genenames=[genenames,genenamesAll(actFilePos)];
                        if(not(isnan(ratiosAll(actFilePos))))
                            ratios=[ratios,ratiosAll(actFilePos)];
                        else
                            ratios=[ratios,1];
                        end
                    end
                    actFilePos=actFilePos+1;

                end
            end
            setnames=[setnames;{setname}];
            disp(setname)
            
            mkdir('storage/expressions');
            expressionIndex=[];

            expressionMatrix =sparse(size(ANOGD,1)*size(ANOGD,2)*size(ANOGD,3),size(genenames,2));

            for i = 1:size(genenames,2)
                
                    urlwrite(strcat('http://api.brain-map.org/api/v2/data/query.xml?criteria=model::SectionDataSet,rma::criteria,[failed$eq%27false%27],products[abbreviation$eq%27Mouse%27],genes[entrez_id$eq%27',sprintf('%d',entrez(i)),'%27]'), 'storage/expression.xml');
                    %urlwrite(strcat('http://api.brain-map.org/api/v2/data/query.xml?criteria=model::SectionDataSet,rma::criteria,[failed$eq%27false%27],products[abbreviation$eq%27Mouse%27],plane_of_section[name$eq%27coronal%27],genes[entrez_id$eq%27',entrez{i},'%27]'), 'expression.xml');
                    root = xmlread('storage/expression.xml');

                    if str2double(root.item(0).getAttribute('num_rows'))<=0
                        disp(sprintf('Could not find %s',genenames{i}))
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
    
                        disp(sprintf('  Get data for %s (%d/%d)',genenames{i},size(expressionIndex,2)+1,size(genenames,2)))
                        for act=1:size(ids,2)
                            id=ids(act);
                            rid=referenceSpace(act);

                            try

                                if(not(exist(sprintf('storage/expressions/%d/density.raw',id), 'file')))
                                    disp(sprintf('      Download %d/%d images for gene from brain-map.org',act,size(ids,2)))
                                    urlwrite(sprintf('http://api.brain-map.org/grid_data/download/%d?include=density',id), 'storage/temp.zip');
                                    unzip('storage/temp.zip',sprintf('storage/expressions/%d',id));
                                else
                                    disp(sprintf('      Load %d/%d images for gene from local',act,size(ids,2)))
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
                                if(rid==10)
                                    expressionMatrixSag=(expressionMatrixSag+expBiggerZero);
                                    amountSag=amountSag+(exp>=0);
                                end
                            catch err
                                disp(sprintf('            Could not open %s: %s',id,getReport(err)))
                            end
                        end

                        %this is the average computation if multiplie sagital/coronal
                        %volumes are available. It excludes "no data" voxels
                        %before the average is computed, the gene expression for
                        %every grid point is standardized with the mean and std
                        %computed by estimate_Mean_STD_of_GeneExpression
                        if(sum(expressionMatrixCor(:)>0) && sum(expressionMatrixSag(:)>0))
                            expressionMatrixCor=(expressionMatrixCor./amountCor);
                            expressionMatrixSag=(expressionMatrixSag./amountSag);
                            expressionMatrixCor=(expressionMatrixCor-expressionMatrixCorMean)./expressionMatrixCorSTD;
                            expressionMatrixSag=(expressionMatrixSag-expressionMatrixSagMean)./expressionMatrixSagSTD;
                            expressionMatrixCor(expressionMatrixCorMean==0)=NaN;
                            expressionMatrixSag(expressionMatrixSagMean==0)=NaN;

                            isnanCor=isnan(expressionMatrixCor);
                            isnanSag=isnan(expressionMatrixSag);

                            divideBy=(2-(isnanCor+isnanSag));
                            expressionMatrixCor(isnanCor & not(isnanSag))=0;
                            expressionMatrixSag(isnanSag & not(isnanCor))=0;

                            expressionMatrixCor(isnanCor & (isnanSag))=NaN;
                            expressionMatrixSag(isnanSag & (isnanCor))=NaN;

                            expressionMatrix(:,i)=(expressionMatrixCor+expressionMatrixSag)./divideBy;
                            expressionIndex= [expressionIndex,i];
                        else if(sum(expressionMatrixCor(:)>0) && sum(expressionMatrixSag(:)==0))
                                expressionMatrixCor=(expressionMatrixCor./amountCor);
                                expressionMatrixCor=(expressionMatrixCor-expressionMatrixCorMean)./expressionMatrixCorSTD;
                                expressionMatrixCor(expressionMatrixCorMean==0)=NaN;
                                expressionMatrix(:,i)=expressionMatrixCor;
                                expressionIndex= [expressionIndex,i];
                            else  if(sum(expressionMatrixCor(:)==0) && sum(expressionMatrixSag(:)>0))
                                    expressionMatrixSag=(expressionMatrixSag./amountSag);
                                    expressionMatrixSag=(expressionMatrixSag-expressionMatrixSagMean)./expressionMatrixSagSTD;
                                    expressionMatrixSag(expressionMatrixSagMean==0)=NaN;
                                    expressionMatrix(:,i)=expressionMatrixSag;
                                    expressionIndex= [expressionIndex,i];
                                end
                            end
                        end
                    end 
            end

          
            save(sprintf('%s/%s.mat',outputDirectory,setname),'expressionMatrix','')
            save(sprintf('%s/%s.mat',outputDirectory,setname),'entrez','-append')
            save(sprintf('%s/%s.mat',outputDirectory,setname),'ratios','-append')
            save(sprintf('%s/%s.mat',outputDirectory,setname),'expressionIndex','-append')
            save(sprintf('%s/%s.mat',outputDirectory,setname),'genenames','-append')

        end
    end

    save(sprintf('%s.mat',outputDirectory),'setnames')
    disp('Computation complete!')
end
