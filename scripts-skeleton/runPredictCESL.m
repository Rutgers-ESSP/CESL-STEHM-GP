% Generate and output predictions
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-01-31 19:49:27 -0500

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
    collinear=[];
    if length(wmodelspec.subamplinear)>0
        collinear=wmodelspec.subamplinear(1);
    end
    [f2s{iii},sd2s{iii},V2s{iii},testlocs{iii},logp(iii),passderivs,invcv]=RegressHoloceneDataSets(wdataset,testsitedef,wmodelspec,thetTGG{jj},trainsub,noiseMasks,testt(subtimes),refyear,collinear);

    % extract linear rate terms
    noiseMasksLin = ones(2,length(thetTGG{jj}));
    noiseMasksLin(1,[setdiff(wmodelspec.subamp,union(wmodelspec.subamplinear,wmodelspec.subampoffset))])=0; %only linear
    noiseMasksLin(2,[setdiff(wmodelspec.subamp,wmodelspec.subampoffset)])=0; %only offset
    noiseMaskLinlabels={'linear','offset'};
    [f2slin{iii},sd2slin{iii},V2slin{iii},testlocslin{iii}]=RegressHoloceneDataSets(wdataset,testsitedef,wmodelspec,thetTGG{jj},trainsub,noiseMasksLin,[0 1800],refyear,collinear,passderivs,invcv);
    
    % output table of all data
    % runTableTGandProxyData;
    
    % generate site plots if dosldecomp set
    if dosldecomp; makeplots_sldecomp(wdataset,f2s{iii},sd2s{iii},V2s{iii},testlocs{iii},labl,1:size(noiseMasks,1),0); end
    
    testreg=testlocs{iii}.reg;
    testsites=testlocs{iii}.sites;
    testX=testlocs{iii}.X;
    testnames2=testlocs{iii}.names2;
    
    %%%%
    
    %runTableRates; % generate table of rates
    %runOutputGSL; % output plots and table of GSL
    runTableDecomp;

    %%%%
 
   save(savefile,'datasets','modelspec','f2s','sd2s','V2s', ...
         'testlocs','logp','testsitedef','trainspecs','thetTGG','ICE5G','noiseMasks','testt','refyear');

    if doMapField
        runMapField; % map spatial field of rates
    end
    runSitePlots;
    runSiteTables;
    
end


runTableMaxToMin;
runLatexTables;
runPlotGSLsWithDifferentPriors;
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
    
    runSitePlotArray; % demonstrate fits
    runSiteTables;
 
    if doMapField
        runMapField; % map spatial field of rates
    end
    
    runSiteSensitivityTests; % perform calculations on sensitivity of prediction to site subsets
    runLatexTablesSiteSens; % latex table
    
end

% check for probability GSL is higher in 2000 than preceeding Common Era with added linear trend


refyear2=2000;

% $$$ Mref = eye(size(testX,1));
% $$$ for i=1:size(testsites,1)
% $$$     
% $$$     sub1=find(testreg==testsites(i,1));
% $$$     sub2=intersect(sub1,find(testX(:,3)==refyear2));
% $$$     
% $$$     Mref(sub1,sub2)=Mref(sub1,sub2)-1;
% $$$ 
% $$$ end
% $$$ Mref=sparse(Mref);

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

trend=0.0; % added trend in mm/yr
N=10000;
sub=find((testX(datsub,3)>=0).*(testX(datsub,3)<=refyear2));

for iii=1:length(f2s)
    %wf=Mref(datsub,datsub)*f2s{iii}(datsub,selmask);
    %wV=Mref(datsub,datsub)*V2s{iii}(datsub,datsub,selmask)*Mref(datsub,datsub)';
    wf=f2s{iii}(datsub,selmask);
    wV=V2s{iii}(datsub,datsub,selmask);
    u=mvnrnd(wf(sub)+trend*(testX(datsub(sub),3)-2000),wV(sub,sub),N);
    u=bsxfun(@minus,u(:,1:end-1),u(:,end));
    sum(max(u,[],2)<0)/N
end
