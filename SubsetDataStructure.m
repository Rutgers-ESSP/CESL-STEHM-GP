function PX=SubsetDataStructure(PX,sub,subS)

% PX=SubsetDataStructure(PX,sub,subS)
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Apr 20 18:54:45 EDT 2014

shortenfields={'datid','time1','time2','meantime','limiting','Y','dY','compactcorr','istg','lat','long','Y0'};
shortenfields=intersect(fieldnames(PX),shortenfields);
for jj=1:length(shortenfields)
    PX.(shortenfields{jj}) =  PX.(shortenfields{jj})(sub);
end
if ~isfield(PX,'Ycv')
    Ycv=diag(PX.dY.^2);
else
    PX.Ycv=sparse(PX.Ycv(sub,sub));
end

shortenfields={'siteid','sitenames'};
for jj=1:length(shortenfields)
    PX.(shortenfields{jj}) =  PX.(shortenfields{jj})(subS);
end
PX.sitecoords=PX.sitecoords(subS,:);
PX.sitelen=zeros(length(PX.siteid),1);
for ii=1:length(PX.siteid)
    PX.sitelen(ii)=sum(PX.datid==PX.siteid(ii));
end