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