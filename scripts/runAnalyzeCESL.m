% Master sea-level rise estimation script
% for Kopp et al., "Temperature-driven global sea-level variability in the Common Era"
%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Nov 28 20:46:33 EST 2015

dosldecomp = 0; % make plots for each site?

% set up  paths

pd=pwd;
addpath([pd '/MFILES']);
addpath([pd '/MFILES/scripts']);
IFILES=[pd '/IFILES'];
addpath(pd)
savefile='~/tmp/CESL';

WORKDIR='151128';
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);

GIAfiles=([pd '/../GIA/RSL4/rsl*.out.gz']);

% exclude all data before -2000 CE
firsttime=-2000;

% read in files and setup covariance model

runImportHoloceneDataSets;
runSetupHoloceneCovariance;
runImportOtherGSLCurves;

% set up and run hyperparameter optimization

trainspecs=[1 2 3 4 5];
trainsets = [2 2 2 2 2]; % use datasets{2}, which include TG, proxy, and flattener
trainfirsttime = -1000; % don't use data before -1000 CE for training

trainlabels={};
for ii=1:length(trainsets)
    trainlabels = {trainlabels{:}, [datasets{trainsets(ii)}.label '_' modelspec(trainspecs(ii)).label]};
end

runTrainModels;

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
avail=ones(length(wdataset.siteid),1); % make more restrictive if you want to exclude some data sets
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

ICE5G.lat=ICE5Glat;
ICE5G.long=ICE5Glon;
ICE5G.gia=ICE5Ggia;

for ii=1:length(testsitedef.sites(:,1))

    testsitedef.GIA = interp2(ICE5G.lat,ICE5G.long,ICE5G.gia,testsitedef.sites(:,2),testsitedef.sites(:,3),'linear');
    testsitedef.GIA(find(testsitedef.sites(:,2)>100))=0;

end

testt = [-1000:20:1800 1810:10:2010]; % ages for regression

% select regression parameters

regressparams=[1 4 5 2 3]; % which trained hyperparameters to use
regresssets=[2 2 2 2 2]; % which data set to use with each
clear regresslabels;
for i=1:length(regresssets)
    regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
end

% generate predictions

for iii=1:length(regresssets)
    ii=regresssets(iii);
    jj=regressparams(iii);
    
    wdataset=datasets{ii};

    wmodelspec = modelspec(trainspecs(jj));

    % create noiseMasks -- set to zero only for amplitudes of
    % terms you want to exclude; if you set to zero for a non-amplitude
    % term it is likely to create problems
    
    noiseMasks = ones(1,length(thetTGG{jj}));
    noiseMasks(1,[ wmodelspec.subampnoise ]  )=0; %without noise
    noiseMasklabels={'denoised'};

    labls{iii}=['_' regresslabels{iii}];
    labl=labls{iii}; disp(labl);

    trainsub = find((wdataset.limiting==0)); % only index points
    
    % add error for compactions, based on compactcorr and the last element of thetTGG{jjj}
    wdataset.dY = sqrt(datasets{ii}.dY.^2 + (thetTGG{jj}(end)*wdataset.compactcorr).^2);
    wdataset.Ycv = datasets{ii}.Ycv + diag(thetTGG{jj}(end)*wdataset.compactcorr).^2;
    
    subtimes=find(testt>=min(union(wdataset.time1,wdataset.time2)));    
    collinear=wmodelspec.subamplinear(1);
    [f2s{iii},sd2s{iii},V2s{iii},testlocs{iii},logp(iii),passderivs,invcv]=RegressHoloceneDataSets(wdataset,testsitedef,wmodelspec,thetTGG{jj},trainsub,noiseMasks,testt(subtimes),refyear,collinear);

    % extract linear rate terms
    noiseMasksLin = ones(2,length(thetTGG{jj}));
    noiseMasksLin(1,[setdiff(wmodelspec.subamp,union(wmodelspec.subamplinear,wmodelspec.subampoffset))])=0; %only linear
    noiseMasksLin(2,[setdiff(wmodelspec.subamp,wmodelspec.subampoffset)])=0; %only offset
    noiseMaskLinlabels={'linear','offset'};
    [f2slin{iii},sd2slin{iii},V2slin{iii},testlocslin{iii}]=RegressHoloceneDataSets(wdataset,testsitedef,wmodelspec,thetTGG{jj},trainsub,noiseMasksLin,[0 1800],refyear,collinear,passderivs,invcv);
    
    % output table of all data
    runTableTGandProxyData;
    
    % generate site plots if dosldecomp set
    if dosldecomp; makeplots_sldecomp(wdataset,f2s{iii},sd2s{iii},V2s{iii},testlocs{iii},labl,[1 2],0); end
    
    testreg=testlocs{iii}.reg;
    testsites=testlocs{iii}.sites;
    testX=testlocs{iii}.X;
    testnames2=testlocs{iii}.names2;
    
    %%%%
    
    runTableRates; % generate table of rates
    runOutputGSL; % output plots and table of GSL
    runPlotOtherGSLCurves; % plot tables of GSL curves from this and other sources

    %%%%
 
   save(savefile,'datasets','modelspec','f2s','sd2s','V2s', ...
         'testlocs','logp','testsitedef','trainspecs','thetTGG','ICE5G','noiseMasks','testt','refyear');

end


runLatexTables;
%runOutputForcingProxies;

% for preferred models, do some additional calculations

for iii=2
    ii=regresssets(iii);
    jj=regressparams(iii);
    wmodelspec = modelspec(trainspecs(jj));

    noiseMasks = ones(1,length(thetTGG{jj}));
    noiseMasks(1,[ wmodelspec.subampnoise ]  )=0; %without noise
    noiseMasklabels={'denoised'};

    wdataset=datasets{ii};
    collinear=wmodelspec.subamplinear(1);

    labls{iii}=['_' regresslabels{iii}];
    labl=labls{iii}; disp(labl);
 
    runMapField; % map spatial field of rates
    runSiteSensitivityTests; % perform calculations on sensitivity of prediction to site subsets
    runLatexTablesSiteSens; % latex table
    
    %calculate distribution minimum-to-maximum value between 0 and 1900 CE
    sub=find(testreg==0);
    samps=lhsnorm(f2s{iii}(sub,1),V2s{iii}(sub,sub,1),1000);
    sub2=find((testt<1900).*(testt>=0));
    mxtomn=max(samps(:,sub2),[],2)-min(samps(:,sub2),[],2);
    quantile(mxtomn,[.05 .95])
end

% map sites
ii=1;
wdataset=datasets{ii};
runMapSites;

% check for probability GSL is higher in 2000 than preceeding Common Era with 0.1 mm/yr linear trend


refyear=2000;
Mref = eye(size(testX,1));
for i=1:size(testsites,1)
    
    sub1=find(testreg==testsites(i,1));
    sub2=intersect(sub1,find(testX(:,3)==refyear));
    
    Mref(sub1,sub2)=Mref(sub1,sub2)-1;

end
Mref=sparse(Mref);

selsitenames={'GSL'};
sitesub=[];
for kk=1:length(selsitenames)
    q=find(strcmpi(selsitenames{kk},testsitedef.names));
    if length(q)>0
        sitesub=[sitesub q(1)];
    end
end
datsub=find(ismember(testreg,testsites(sitesub,1)));
selmask=1;

trend=0.1; % added trend in mm/yr
N=10000;
sub=find((testX(datsub,3)>=0).*(testX(datsub,3)<2000));

for iii=1:length(f2s)
    %wf=Mref(datsub,datsub)*f2s{iii}(datsub,selmask);
    %wV=Mref(datsub,datsub)*V2s{iii}(datsub,datsub,selmask)*Mref(datsub,datsub)';
    wf=f2s{iii}(datsub,selmask);
    wV=V2s{iii}(datsub,datsub,selmask);
    u=mvnrnd(wf(sub)+trend*(testX(datsub(sub),3)-2000),wV(sub,sub),N);
    sum(max(u,[],2)<0)/N
end
