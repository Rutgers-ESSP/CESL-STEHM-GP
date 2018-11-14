
% 2018-06-15 19:12:29 -0400

for iii=1:length(regressparams)

    ii=regresssets(iii);
    jj=regressparams(iii);

    wdataset=datasets{ii};

    wmodelspec = modelspec(trainspecs(jj));

    labls{iii}=['_' regresslabels{iii}];
    labl=labls{iii}; disp(labl);

    trainsub = find((wdataset.limiting==0)); % only index points

    % create noiseMasks -- set to zero only for amplitudes of
    % terms you want to exclude; if you set to zero for a non-amplitude
    % term it is likely to create problems

    noiseMasks = ones(6,length(thetTGG{jj}));
    noiseMasks(1,[ wmodelspec.subampnoise ]  )=0; %without noise
    noiseMasks(2,[ wmodelspec.subampoffset wmodelspec.subampnoise wmodelspec.subampregmat wmodelspec.subampglobal ]  )=0;
    noiseMasks(3,[ wmodelspec.subampoffset wmodelspec.subampnoise wmodelspec.subamplinear ]  )=0;
    noiseMasks(4,[ wmodelspec.subampoffset wmodelspec.subampnoise wmodelspec.subampregmat wmodelspec.subamplinear ]  )=0;
    noiseMasks(5,[ wmodelspec.subampoffset wmodelspec.subampnoise wmodelspec.subampglobal wmodelspec.subamplinear  wmodelspec.subampregmat2]  )=0; %without noise
    noiseMasks(6,[ wmodelspec.subampoffset wmodelspec.subampnoise wmodelspec.subampglobal wmodelspec.subamplinear  wmodelspec.subampregmat1]  )=0; %without noise
    noiseMasklabels={'denoised','lin','nonlin','global','reg-nonlin','local'};

    % add error for compactions, based on compactcorr and the last element of thetTGG{jjj}
    wdataset.dY = sqrt(datasets{ii}.dY.^2 + (thetTGG{jj}(end)*wdataset.compactcorr).^2);
    wdataset.Ycv = datasets{ii}.Ycv + diag(thetTGG{jj}(end)*wdataset.compactcorr).^2;

    subtimes=find(testt>=min(union(wdataset.time1,wdataset.time2)));    
    collinear=wmodelspec.subamplinear(1);
    [f2s{iii},sd2s{iii},V2s{iii},testlocs{iii},logp(iii),passderivs,invcv]=RegressHoloceneDataSets(wdataset,testsitedef,wmodelspec,thetTGG{jj},trainsub,noiseMasks,testt(subtimes),refyear,collinear);

    testreg=testlocs{iii}.reg;
    testsites=testlocs{iii}.sites;
    testX=testlocs{iii}.X;
    testnames2=testlocs{iii}.names2;

    save(savefile,'datasets','modelspec','f2s','sd2s','V2s', ...
         'testlocs','logp','testsitedef','trainspecs','thetTGG','ICE5G','noiseMasks','testt','refyear');

    % extract linear rate terms
    noiseMasksLin = ones(1,length(thetTGG{jj}));
    noiseMasksLin(1,[setdiff(wmodelspec.subamp,union(wmodelspec.subamplinear,wmodelspec.subampoffset))])=0; %only linear
    noiseMasksLin(2,[setdiff(wmodelspec.subamp,wmodelspec.subampoffset)])=0; %only offset
    noiseMaskLinlabels={'linear','offset'};
    [f2slin{iii},sd2slin{iii},V2slin{iii},testlocslin{iii}]=RegressHoloceneDataSets(wdataset,testsitedef,wmodelspec,thetTGG{jj},trainsub,noiseMasksLin,[0 1800],refyear,collinear,passderivs,invcv);
    
    % Output table of all data points, including model prediction at sites
    runTableTGandProxyData;

    
   % generate site plots if dosldecomp set
   
   if dosldecomp
        figure;
        set(gcf,'PaperOrientation','landscape','PaperPosition',[0.5 0.5 110 6]);
        makeplots_sldecomp(wdataset,f2s{iii},sd2s{iii},V2s{iii},testlocs{iii},labl,[1 2 3 4 5 6 7],0);
    end

    % decompsition table

    runTableDecomp;

    % For each site, make a plot
    runSitePlots;

    % For each site, make a table

    runSiteTables;

    % Generate maps of rate fields.
    runMapField;
    
    % Generate table of rates at each site, and differences in rate at each site
    runTableRates; 

    % Output table and plots of GSL.
    runOutputGSL;

end