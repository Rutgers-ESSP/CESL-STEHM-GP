function PX=SubsetDataStructure(PX,sub,subS)

% PX=SubsetDataStructure(PX,sub,subS)
%
% Construct a subset data structure from the sea-level data set
% where sub identifies the individual data points to include and
% subS the sites to include.
%
% EXAMPLE:
%   sub=find(PX.time1>=firsttime);
%   sub=intersect(sub,find(abs(PX.lat)<=90));
%   subS=find(abs(PX.sitecoords(:,1))<=90);
%   PX=SubsetDataStructure(PX,sub,subS);
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Nov 16 17:08:35 EST 2014

shortenfields={'datid','time1','time2','meantime','limiting','Y','dY','compactcorr','istg','lat','long','Y0','dt','obsGISfp','GIAproj'};
shortenfields=intersect(fieldnames(PX),shortenfields);
for jj=1:length(shortenfields)
    if isfield(PX,shortenfields{jj})
        PX.(shortenfields{jj}) =  PX.(shortenfields{jj})(sub);
    end
end
if ~isfield(PX,'Ycv')
    Ycv=diag(PX.dY.^2);
else
    PX.Ycv=sparse(PX.Ycv(sub,sub));
end

shortenfields={'siteid','sitenames','siteGISfp','siteGIA'};
for jj=1:length(shortenfields)
    if isfield(PX,shortenfields{jj})
        PX.(shortenfields{jj}) =  PX.(shortenfields{jj})(subS);
    end
    
end
PX.sitecoords=PX.sitecoords(subS,:);
PX.sitelen=zeros(length(PX.siteid),1);
for ii=1:length(PX.siteid)
    PX.sitelen(ii)=sum(PX.datid==PX.siteid(ii));
end