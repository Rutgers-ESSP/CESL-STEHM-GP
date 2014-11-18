% Master script for Common Era proxy + tide gauge sea-level analysis
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Nov 09 09:50:00 EST 2014

% to do items:
% 1) MCMC

dosldecomp = 1;

pd=pwd;
%addpath('~/Dropbox/Code/TGAnalysis/MFILES');
addpath([pd '/MFILES']);
addpath([pd '/MFILES/scripts']);
IFILES=[pd '/IFILES'];
addpath(pd)
savefile='~/tmp/CESL';

WORKDIR='141117';
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);

GIAfiles=([pd '/../GIA/RSL4/rsl*.out.gz']);

%
firsttime=-1000;

runSetupHoloceneCovariance;
runImportHoloceneDataSets;

save(savefile,'datasets','modelspec');

% we will run: GLMW (5), LMW (2), GLW (4), GMW (9), LW (7), GW (8), MW (10)

trainspecs = [5 2 4 9 7 8 10];
trainsets = [1*ones(size(trainspecs))];
trainspecs = [trainspecs trainspecs];
trainfirsttime = -1000;

trainlabels={};
for ii=1:length(trainsets)
    trainlabels = {trainlabels{:}, [datasets{trainsets(ii)}.label '_' modelspec(trainspecs(ii)).label]};
end

runTrainModels;
% runCrossValidateModels;

save thetTGG thetTGG trainsubsubset
save(savefile);


%% now do a regression
% define prediction sites

clear oldest oldestnear youngest testsitedef;
wdataset = datasets{1};
distfrom = dDist(wdataset.sitecoords,wdataset.sitecoords);
for ii=1:length(wdataset.sitenames)
    sub = find(wdataset.datid==wdataset.siteid(ii));
    if length(sub)>0
        oldest(ii)=min(union(wdataset.time1(sub),wdataset.time2(sub)));
        youngest(ii)=max(union(wdataset.time1(sub),wdataset.time2(sub)));
    else
        oldest(ii)=2000;
        youngest(ii)=2014;
    end
end
for ii=1:length(wdataset.sitenames)
    sub=find(distfrom(ii,:)<5);
    oldestnear(ii) = min(oldest(sub));
end

clear testsitedef;
testsitedef.sites=[0 1e6 1e6];
testsitedef.names={'GSL'};
testsitedef.names2={'GSL'};
testsitedef.firstage=min(oldest);
testsitedef.oldest=min(oldest);
testsitedef.youngest=2014;

dosub=find(wdataset.siteid>0);
for ii=dosub(:)'
    sub=find(wdataset.datid==wdataset.siteid(ii));
    if length(sub)>0
        testsitedef.sites(end+1,:)=mean([wdataset.datid(sub) wdataset.lat(sub) wdataset.long(sub)],1);
        testsitedef.names2={testsitedef.names2{:}, wdataset.sitenames{ii}};
        
        sublett=setdiff(1:length(wdataset.sitenames{ii}),strfind(wdataset.sitenames{ii},' '));
        testsitedef.names={testsitedef.names{:}, wdataset.sitenames{ii}(sublett)};
        testsitedef.firstage = [testsitedef.firstage oldestnear(ii)];
        testsitedef.oldest = [testsitedef.oldest oldest(ii)];
        testsitedef.youngest = [testsitedef.youngest youngest(ii)];
    end
end

GISfpt.lat=GISfplat;
GISfpt.long=GISfplong;
GISfpt.fp=GISfp;

ICE5G.lat=ICE5Glat;
ICE5G.long=ICE5Glon;
ICE5G.gia=ICE5Ggia;

for ii=1:length(testsitedef.sites(:,1))
    testsitedef.GISfp = interp2(GISfpt.long,GISfpt.lat,GISfpt.fp,testsitedef.sites(:,3),testsitedef.sites(:,2),'linear');
    testsitedef.GISfp(find(testsitedef.sites(:,2)>100))=1;

    testsitedef.GIA = interp2(ICE5G.lat,ICE5G.long,ICE5G.gia,testsitedef.sites(:,2),testsitedef.sites(:,3),'linear');
    testsitedef.GIA(find(testsitedef.sites(:,2)>100))=0;

end
GISfpt.lat=GISfplat;
GISfpt.long=GISfplong;
GISfpt.fp=GISfp;

ICE5G.lat=ICE5Glat;
ICE5G.long=ICE5Glon;
ICE5G.gia=ICE5Ggia;

for ii=1:length(testsitedef.sites(:,1))
    testsitedef.GISfp = interp2(GISfpt.long,GISfpt.lat,GISfpt.fp,testsitedef.sites(:,3),testsitedef.sites(:,2),'linear');
    testsitedef.GISfp(find(testsitedef.sites(:,2)>100))=1;

    testsitedef.GIA = interp2(ICE5G.lat,ICE5G.long,ICE5G.gia,testsitedef.sites(:,2),testsitedef.sites(:,3),'linear');
    testsitedef.GIA(find(testsitedef.sites(:,2)>100))=0;

end


testt = [-1000:20:2000 2010];

% select regression parameters

regressparams=1:length(thetTGG);
regresssets=ones(size(regressparams));
clear regresslabels;
for i=1:length(regresssets)
    regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
end

% run predictions

for iii=1:length(regresssets)
    ii=regresssets(iii);
    jj=regressparams(iii);
    wmodelspec = modelspec(trainspecs(jj));
    

    noiseMasks = ones(4,length(thetTGG{trainspecs(jj)}));
    noiseMasks(1,[ wmodelspec.subampnoise]  )=0; %without linear
    noiseMasks(2,[wmodelspec.subamplinear wmodelspec.subampoffset wmodelspec.subampnoise]  )=0; %without linear
    noiseMasks(3,[setdiff(wmodelspec.subamp,wmodelspec.subamplinear)])=0; %only linear
    noiseMasks(4,[setdiff(wmodelspec.subamp,wmodelspec.subampregmat)])=0; %only regional Matern
    noiseMasklabels={'denoised','nonlin','linear','regmat'};

    wdataset=datasets{ii};

    labls{iii}=['_' regresslabels{iii}];
    disp(labls{iii});

    trainsub = find((wdataset.limiting==0)); % only index points
    wdataset.dY = sqrt(datasets{ii}.dY.^2 + (thetTGG{jj}(end)*wdataset.compactcorr).^2);
    wdataset.Ycv = datasets{ii}.Ycv + diag(thetTGG{jj}(end)*wdataset.compactcorr).^2;
    subtimes=find(testt>=min(union(wdataset.time1,wdataset.time2)));
    
    collinear=wmodelspec.subamplinear(1);
    [f2s{ii,jj},sd2s{ii,jj},V2s{ii,jj},testlocs{ii,jj},logp(ii,jj),passderivs,invcv]=RegressHoloceneDataSets(wdataset,testsitedef,wmodelspec,thetTGG{jj},trainsub,noiseMasks,testt(subtimes),refyear,collinear);

    labl=labls{iii}; disp(labl);
    
   runTableTGandProxyData;
    
    if dosldecomp; makeplots_sldecomp(wdataset,f2s{ii,jj},sd2s{ii,jj},V2s{ii,jj},testlocs{ii,jj},labl,[1 2],0); end
    
    testreg=testlocs{ii,jj}.reg;
    testsites=testlocs{ii,jj}.sites;
    testX=testlocs{ii,jj}.X;
    testnames2=testlocs{ii,jj}.names2;
    
    %%%%
    
    runMapField;
    runTableRates;
    runOutputGSL;
    runFingerprintAnalysis;
    % runGIARateComparison;
    
    if iii=1
        %runSensitivityTests;
        runSiteSensitivityTests;
        runPlotOtherGSLCurves;
    end
    

    %%%%
 
   save(savefile,'datasets','modelspec','f2s','sd2s','V2s', ...
         'testlocs','logp','testsitedef','trainspecs','thetTGG','GISfpt','ICE5G','noiseMasks','testt','refyear');

end

%runOutputForcingProxies;