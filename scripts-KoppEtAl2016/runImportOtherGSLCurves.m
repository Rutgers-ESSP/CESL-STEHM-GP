% Create pseudo-GSL curve from North Carolina proxy data
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 01 14:28:34 EST 2016

% extract the detrended North Carolina sea-level curve

ms=modelspecNC;
sitesub=find(strcmpi('North Carolina-Sand Point',wdataset.sitenames));
sitesub=union(sitesub,find(strcmpi('North Carolina-Tump Point',wdataset.sitenames)));
datsub=find(ismember(wdataset.datid,wdataset.siteid(sitesub)));

nmms=ones(size(thetNC)); nmms(3)=0;

testtNC = [-1000:20:2000 2010];

wtestsitedef.names={'NC_GSL'};
wtestsitedef.names2={'NC_GSL'};
wtestsitedef.firstage=-1000;
wtestsitedef.oldest=-1000;
wtestsitedef.youngest=2014;


wtestsitedef.sites(:,2:3)=mean(wdataset.sitecoords(sitesub,:),1);
[NC_pseudoGSL,NC_pseudoGSLsd,NC_pseudoGSLV,NClocs]=RegressHoloceneDataSets(wdataset,wtestsitedef,ms,thetNC,datsub,nmms,testtNC,refyear,3);
NC_yrs=NClocs.X(:,3);

[NC_pseudoGSL,NC_pseudoGSLV,NC_pseudoGSLsd,trendspec]=DetrendSLReconstruction(NC_pseudoGSL,NC_pseudoGSLV+ eye(size(NC_pseudoGSLV))*50^2,NClocs.sites,NClocs.reg,NC_yrs,0,1700,2000,100);
Mrefms=eye(length(NC_pseudoGSL));
subms=find((NC_yrs<=1800).*(NC_yrs>=0));
Mrefms(:,subms)=Mrefms(:,subms)-1/length(subms);
NC_pseudoGSL=Mrefms*NC_pseudoGSL;
NC_pesudoGSLsd=sqrt(diag(Mrefms*NC_pseudoGSLV*Mrefms'));

% output NC_pseudoGSL
[hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(NClocs.X(:,3),NClocs.reg,0,NC_pseudoGSL,NC_pseudoGSLV,colrs,NClocs.X(:,3),testt(end),0,100,{'NC_GSL'});

fid=fopen(['NCpseudoGSL.tsv'],'w');
fprintf(fid,outtable);


fclose(fid);


% output GSL covariance

fid=fopen(['NCpseudoGSL_cov.tsv'],'w');
fprintf(fid,'mm^2');
fprintf(fid,'\t%0.0f',NClocs.X(:,3));
fprintf(fid,'\n');
for ppp=1:length(NClocs.X)
    fprintf(fid,'%0.0f',NClocs.X(ppp,3));
    fprintf(fid,'\t%0.8e',NC_pseudoGSLV(ppp,:));
    fprintf(fid,'\n');
end

fclose(fid);


