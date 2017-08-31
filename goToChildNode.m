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

function [ attribute, childnodes ] = goToChildNode( node,childNodeName, start )
    act=0;
    for actnode=0:node.getChildNodes.getLength
   
        if strcmp(node.getChildNodes.item(actnode),sprintf('[%s: null]',childNodeName))
            if(act>=start)
               attributeText = char(node.getChildNodes.item(actnode).item(0));
               attribute = attributeText(8:end-1);
               childnodes = node.getChildNodes.item(actnode).getChildNodes;
               
               return;
            end
            act=act+1;
        end
    end
end

