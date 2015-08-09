% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Aug 09 11:36:57 EDT 2015

labl=labls{iii}; disp(labl);

firstyears0 = [0];
lastyears0 = [1800];
firstyears1=[ 0:200:1800 1900 ];
lastyears1= [200:200:2000 2000];
winvcv=[];
minpriorsd = 0.5; % mm/y
dynamicsd = 0.5; % mm/y;

% load fingerprint patterns
if ~exist('fpSTORED','var')
    pd=pwd;
    [fp,fpname,lo,la] = readFingerprint(fullfile(IFILES,'FPRINT'));
    cd(pd);
    lo(end+1)=360;
    fp(:,end+1,:)=fp(:,1,:);
    fp=fp*1000;
    fpSTORED=fp;
else
    fp=fpSTORED;
end

fpwtdat=importdata(fullfile(IFILES,'FPRINT/fingerprint_region_map.csv'));
clear fpprior wfpbasis fplabel;
for qqq=1:length(fpname)
    fpsub = find(strcmpi(fpname{qqq},fpwtdat.textdata(2:end,3)));
    if length(fpsub)==1
        fplabel{qqq}=fpwtdat.textdata{fpsub+1,2};
        fpprior(qqq)=fpwtdat.data(fpsub,2);
    else
        fpprior(qqq)=0;
    end
    
    wfpbasis(:,qqq)=interp2(lo',la,fp(:,:,qqq),mod(testsites(:,3),360),testsites(:,2));
end
subGSL=find(testsites(:,1)==0);
wfpbasis(subGSL,:)=1;
fppriormean=zeros(size(fpprior));

wad=dDist(testsites(:,2:3),testsites(:,2:3));
wtestcv=(wad<360).*kMat1(wad,[1 thetTGG{iii}(7)])*dynamicsd^2;

subfp=find(fpprior>.01*mean(fpprior));
fpprior=fpprior(subfp);
fplabel=fplabel(subfp);
wfpbasis=wfpbasis(:,subfp);

fpprior(end+1)=1;
fplabel{end+1}='uniform';
wfpbasis(:,end+1)=ones(size(wfpbasis(:,1)));

fpclusters={
    {'GIS'},
    {'WAIS','EAIS'},
    {'WAIS','EAIS','uniform'},
    {'GIS','Iceland','Baffin','Ellesmere','Svalbard','Alaska','Western_Can_US','NOVAYA_ZEMLYA','HMA_1'},
    {'WAIS','EAIS','Low_Lat_Andes','Patagonia','New_Zealand'},
    {'WAIS','EAIS','Low_Lat_Andes','Patagonia','New_Zealand','uniform'}
    {'Iceland','Baffin','Ellesmere','Svalbard','Alaska','Western_Can_US','NOVAYA_ZEMLYA','HMA_1','Low_Lat_Andes','Patagonia','New_Zealand'},
           };
fpclusterlabels={'GIS','AIS','AIS+uni','NH','SH','SH+uni','GICexGISAIS'};
fpclusters{end+1}=fplabel;
fpclusterlabels{end+1} = 'all';

Mcluster=zeros(length(fpclusters),length(fplabel));
for rrr=1:length(fpclusters)
    Mcluster(rrr,find(ismember(fplabel,fpclusters{rrr})))=1;
end

clear basisest basisestsd basisclustest basisclustestcv basisclustestsd GSLest GSLestsd priorclustmean priorclustsd priorbasissd;
for qqq=1:length(firstyears1)
    disp(sprintf('%0.0f--%0.0f',[firstyears1(qqq) lastyears1(qqq)]));
    firstyears=[firstyears0 firstyears1(qqq)];
    lastyears=[lastyears0 lastyears1(qqq)];

    wf=f2s{iii}(:,1);
    wV=V2s{iii}(:,:,1);
    [wfslope,wsdslope,wVslope,wfy,wly,wdosite,wfslopediff,wsdslopediff,wVslopediff,wdiffsite,wdiffplus,wdiffless]=SLRateCompareCrosssite(wf,wV,testsites,testreg,testX,firstyears,lastyears);
    
    sitesub=find((wsdslopediff>0));
    
    basissub=1:size(wfpbasis,2);
    priorbasiscv=diag(fpprior(basissub))*max(minpriorsd,wfslopediff(subGSL)).^2;
    [basisest(:,qqq),basisestcv,~,~,~,~,~,~,~,~,~,~,winvcv] = GaussianProcessRegressionWithBasis([],wfslopediff(sitesub),[],wtestcv(sitesub,sitesub)+wVslopediff(sitesub,sitesub),wtestcv(sitesub,sitesub),wtestcv(sitesub,sitesub),wfpbasis(sitesub,basissub)',wfpbasis(sitesub,basissub)',fppriormean(basissub),priorbasiscv);
    priorbasissd(:,qqq)=sqrt(diag(priorbasiscv));
    
    basisestsd(:,qqq)=sqrt(diag(basisestcv));
    basisclustest(:,qqq)=Mcluster*basisest(:,qqq);
    basisclustestcv=Mcluster*basisestcv*Mcluster';
    basisclustestsd(:,qqq)=sqrt(diag(basisclustestcv));
    
    GSLest(qqq) = wfslopediff(subGSL);
    GSLestsd(qqq) = wsdslopediff(subGSL);

    priorclustmean(:,qqq)=Mcluster*fppriormean(basissub)';
    priorclustcv=Mcluster*priorbasiscv*Mcluster';
    priorclustsd(:,qqq)=sqrt(diag(priorclustcv));
end



fid=fopen(['fingerprints' labl '.tsv'],'w');
fprintf(fid,'All based on sea-level rates relative to %0.0f--%0.0f CE\n',[firstyears0 lastyears0]);
for qqq=1:length(firstyears1);
    fprintf(fid,'\t%0.0f',firstyears1(qqq));
    fprintf(fid,'--%0.0f (mm/y)',lastyears1(qqq));
    fprintf(fid,'\t2s\tP>0');
end
fprintf(fid,'\n');

fprintf(fid,'GSL');
fprintf(fid,'\tNaN\tNaN\tNaN');
for qqq=1:length(firstyears1)
    fprintf(fid,'\t%0.3f',[GSLest(qqq) 2*GSLestsd(qqq) normcdf(GSLest(qqq)/GSLestsd(qqq))]);
end
fprintf(fid,'\n');

for rrr=1:length(fplabel)
    fprintf(fid,fplabel{rrr});
    for qqq=1:length(firstyears1)
        fprintf(fid,'\t%0.3g',[basisest(rrr,qqq) 2*basisestsd(rrr,qqq) normcdf(basisest(rrr,qqq)/basisestsd(rrr,qqq))]);
    end
    fprintf(fid,'\n');
end

for rrr=1:length(fpclusters)
    fprintf(fid,fpclusterlabels{rrr});
    for qqq=1:length(firstyears1)
        fprintf(fid,'\t%0.3f',[basisclustest(rrr,qqq) 2*basisclustestsd(rrr,qqq) normcdf(basisclustest(rrr,qqq)/basisclustestsd(rrr,qqq))]);
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\nPrior\n');
for rrr=1:length(fplabel)
    fprintf(fid,fplabel{rrr});
    for qqq=1:length(firstyears1)
        fprintf(fid,'\t%0.3g',[fppriormean(rrr) 2*priorbasissd(rrr,qqq) normcdf(fppriormean(rrr)/priorbasissd(rrr,qqq))]);
    end
    fprintf(fid,'\n');
end

for rrr=1:length(fpclusters)
    fprintf(fid,fpclusterlabels{rrr});
    for qqq=1:length(firstyears1)
        fprintf(fid,'\t%0.3f',[priorclustmean(rrr,qqq) 2*priorclustsd(rrr,qqq) normcdf(priorclustmean(rrr,qqq)/priorclustsd(rrr,qqq))]);
    end
    fprintf(fid,'\n');
end


fprintf(fid,'Local variability');
for qqq=1:length(firstyears1)
    fprintf(fid,'\t%0.3g',[0 2*sqrt(max(diag(wtestcv))) NaN]);
end
fprintf(fid,'\n');

fclose(fid);