% Master sea-level rise estimation script
% for Kopp et al., "Temperature-driven global sea-level variability in the Common Era"
%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Dec 25 09:52:09 EST 2015

dosldecomp = 0; % make plots for each site?

% set up  paths

pd=pwd;
addpath([pd '/MFILES']);
addpath([pd '/MFILES/scripts']);
IFILES=[pd '/IFILES'];
addpath(pd)
savefile='~/tmp/CESL';

WORKDIR='151225';
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

trainspecs=[1 2 3 4 5 9];
trainsets = [2 2 2 2 2 2]; % use datasets{2}, which include TG, proxy, and flattener
trainfirsttime = -1000; % don't use data before -1000 CE for training

trainlabels={};
for ii=1:length(trainsets)
    trainlabels = {trainlabels{:}, [datasets{trainsets(ii)}.label '_' modelspec(trainspecs(ii)).label]};
end

runTrainModels;

save thetTGG thetTGG trainsubsubset
save(savefile);

thetTGG0=thetTGG;

% map sites
runMapSites;



%% now do a regression
% define prediction sites

doMapField = 1;
runPredictCESL

doMapField = 0;
if ~exist('scaledamps','dir')
    mkdir('scaledamps');
end

cd('scaledamps');
% for sensitivity tests in which we scale prior for regional and global non-linear variability
multiplyamplitudes=3;
if multiplyamplitudes ~= 1
    for ii=1:length(trainsets)
        wms=modelspec(trainspecs(ii));
        toamp = [wms.subampglobal wms.subampregmat];
        thetTGG{ii}(toamp)= thetTGG{ii}(toamp)*multiplyamplitudes;
    end   
end
runPredictCESL;
cd('..');
