% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Feb 19 18:41:08 EST 2015
%

dosldecomp = 0;

pd=pwd;
%addpath('~/Dropbox/Code/TGAnalysis/MFILES');
addpath([pd '/MFILES']);
addpath([pd '/MFILES/scripts']);
IFILES=[pd '/IFILES'];
addpath(pd)
savefile='~/tmp/CESL';

WORKDIR='150216';
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);

GIAfiles=([pd '/../GIA/RSL4/rsl*.out.gz']);

%
firsttime=-1000;

runImportHoloceneDataSets;
runSetupHoloceneCovariance;
runImportOtherGSLCurves;

trainspecs=[1 2 3 4 5];
trainsets = [1 1 1 1 1];
trainfirsttime = -1000;

trainlabels={};
for ii=1:length(trainsets)
    trainlabels = {trainlabels{:}, [datasets{trainsets(ii)}.label '_' modelspec(trainspecs(ii)).label]};
end

runTrainModels;
% runCrossValidateModels;

% add weakly informative prior
% $$$ trainspecs=[trainspecs 1];
% $$$ trainsets=[trainsets NaN];
% $$$ trainlabels={trainlabels{:},['wip_' modelspec(trainspecs(end)).label]};
% $$$ thetTGG{end+1}=[200 300 1.3 8 200 120 15 20 40 0.1];
% $$$ logp(end+1)=NaN;

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
    sub=find(distfrom(ii,:)<.5);
    if wdataset.sitecoords(ii,1)>360
        sub=ii;
    end
    oldestnear(ii) = min(oldest(sub));
end

clear testsitedef;
testsitedef.sites=[0 1e6 1e6];
testsitedef.names={'GSL'};
testsitedef.names2={'GSL'};
testsitedef.firstage=min(oldest);
testsitedef.oldest=min(oldest);
testsitedef.youngest=2014;

testsitedefGSL=testsitedef;

dosub=find(wdataset.siteid>0);
avail=ones(length(wdataset.siteid),1);
for ii=dosub(:)'
    if avail(ii)
        sub=find((wdataset.datid==wdataset.siteid(ii)));
        if length(sub)>0
            testsitedef.sites(end+1,:)=mean([wdataset.datid(sub) wdataset.lat(sub) wdataset.long(sub)],1);
            wdist = dDist(testsitedef.sites(end,2:3),wdataset.sitecoords);
            avail=avail.*(wdist>.2); % keep sites at least .2 degrees apart
            testsitedef.names2={testsitedef.names2{:}, wdataset.sitenames{ii}};
            
            sublett=setdiff(1:length(wdataset.sitenames{ii}),strfind(wdataset.sitenames{ii},' '));
            testsitedef.names={testsitedef.names{:}, wdataset.sitenames{ii}(sublett)};
            testsitedef.firstage = [testsitedef.firstage oldestnear(ii)];
            testsitedef.oldest = [testsitedef.oldest oldest(ii)];
            testsitedef.youngest = [testsitedef.youngest youngest(ii)];
        end
    end
    
end
% $$$ 
% $$$ GISfpt.lat=GISfplat;
% $$$ GISfpt.long=GISfplong;
% $$$ GISfpt.fp=GISfp;

ICE5G.lat=ICE5Glat;
ICE5G.long=ICE5Glon;
ICE5G.gia=ICE5Ggia;

for ii=1:length(testsitedef.sites(:,1))
% $$$     testsitedef.GISfp = interp2(GISfpt.long,GISfpt.lat,GISfpt.fp,testsitedef.sites(:,3),testsitedef.sites(:,2),'linear');
% $$$     testsitedef.GISfp(find(testsitedef.sites(:,2)>100))=1;

    testsitedef.GIA = interp2(ICE5G.lat,ICE5G.long,ICE5G.gia,testsitedef.sites(:,2),testsitedef.sites(:,3),'linear');
    testsitedef.GIA(find(testsitedef.sites(:,2)>100))=0;

end

testt = [-1000:20:2000 2010];

% select regression parameters

regressparams=[1 4 5 2 3];
regresssets=[1 1 1 1 1];
clear regresslabels;
for i=1:length(regresssets)
    regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
end

% run predictions

for iii=1:length(regresssets)
    ii=regresssets(iii);
    jj=regressparams(iii);
    wmodelspec = modelspec(trainspecs(jj));

    noiseMasks = ones(1,length(thetTGG{jj}));
    noiseMasks(1,[ wmodelspec.subampnoise]  )=0; %without noise
    noiseMasklabels={'denoised'};

    wdataset=datasets{ii};

    labls{iii}=['_' regresslabels{iii}];
    labl=labls{iii}; disp(labl);

    trainsub = find((wdataset.limiting==0)); % only index points
    wdataset.dY = sqrt(datasets{ii}.dY.^2 + (thetTGG{jj}(end)*wdataset.compactcorr).^2);
    wdataset.Ycv = datasets{ii}.Ycv + diag(thetTGG{jj}(end)*wdataset.compactcorr).^2;
    subtimes=find(testt>=min(union(wdataset.time1,wdataset.time2)));
    
    collinear=wmodelspec.subamplinear(1);
    [f2s{iii},sd2s{iii},V2s{iii},testlocs{iii},logp(iii),passderivs,invcv]=RegressHoloceneDataSets(wdataset,testsitedef,wmodelspec,thetTGG{jj},trainsub,noiseMasks,testt(subtimes),refyear,collinear);

    noiseMasksLin = ones(2,length(thetTGG{jj}));
    noiseMasksLin(1,[setdiff(wmodelspec.subamp,union(wmodelspec.subamplinear,wmodelspec.subampoffset))])=0; %only linear
    noiseMasksLin(2,[setdiff(wmodelspec.subamp,wmodelspec.subampoffset)])=0; %only offset
    noiseMaskLinlabels={'linear','offset'};
    [f2slin{iii},sd2slin{iii},V2slin{iii},testlocslin{iii}]=RegressHoloceneDataSets(wdataset,testsitedef,wmodelspec,thetTGG{jj},trainsub,noiseMasksLin,[0 1800],refyear,collinear,passderivs,invcv);

    
    runTableTGandProxyData;
    
    if dosldecomp; makeplots_sldecomp(wdataset,f2s{iii},sd2s{iii},V2s{iii},testlocs{iii},labl,[1 2],0); end
    
    testreg=testlocs{iii}.reg;
    testsites=testlocs{iii}.sites;
    testX=testlocs{iii}.X;
    testnames2=testlocs{iii}.names2;
    
    %%%%
    
    runTableRates;
    runOutputGSL;
    % runFingerprintAnalysis;
    % runGIARateComparison;
    runPlotOtherGSLCurves;

    %%%%
 
   save(savefile,'datasets','modelspec','f2s','sd2s','V2s', ...
         'testlocs','logp','testsitedef','trainspecs','thetTGG','ICE5G','noiseMasks','testt','refyear');

end


runLatexTables;
%runOutputForcingProxies;

for iii=1:length(regresssets)
    ii=regresssets(iii);
    jj=regressparams(iii);
    wmodelspec = modelspec(trainspecs(jj));

    noiseMasks = ones(1,length(thetTGG{jj}));
    noiseMasks(1,[ wmodelspec.subampnoise]  )=0; %without noise
    noiseMasklabels={'denoised'};

    wdataset=datasets{ii};

    labls{iii}=['_' regresslabels{iii}];
    labl=labls{iii}; disp(labl);
   
    runMapField;
    runSiteSensitivityTests;

end
runLatexTablesSiteSens;

ii=1;
wdataset=datasets{ii};
runMapSites;