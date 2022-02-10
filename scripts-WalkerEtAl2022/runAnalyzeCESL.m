% Last updated by Bob Kopp, 2018-09-25 14:02:54 -0400

% set up  paths

rootdir='~/Dropbox/Code/CESL-STEHM-GP';

pd=pwd;
addpath([rootdir '/MFILES']); % add the path with the core MFILES
addpath([rootdir '/scripts']); % add the path with the script files
IFILES=[rootdir '/IFILES']; % point to the directory with the data files
IFILES2=[rootdir '/IFILES-working'];
addpath(pd);
savefile='~/tmp/CESL'; % point to the .mat file you want the analysis to be backed up into 

WORKDIR=[pd '/workdir-010621ToE']; % point to the working directory you want your tables and figures dumped from
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);

PXdatafile=fullfile(IFILES,'RSL_All_17Mar2020.tsv');
latlim=[-90 90]; longlim=[-180 180];

runImportCESLDataSets;
runMapSites;

runSetupCESLCovariance;

% set up and run hyperparameter optimization

%trainspecs=1:length(modelspec); % identify the different model specifications (as defined in runSetupCESLCovariance) 
                          % that will be trained
trainspecs=[9];
trainsets =ones(size(trainspecs))*1; % identify the different datasets (in the structure datasets created by runImportCESLDataSeta)
                           % that will be used for each specification
trainfirsttime = -2000; % don't use data before trainfirsttime (default: -1000 CE) for training

trainlabels={};
for ii=1:length(trainsets)
    trainlabels = {trainlabels{:}, [datasets{trainsets(ii)}.label '_' modelspec(trainspecs(ii)).label]};
end

runTrainModels;

save thetTGG thetTGG trainsubsubset
save(savefile);

thetTGG0=thetTGG;

% now do a prediction
testt = [0:20:1800 1810:10:2010]; % ages for regression

runSelectPredictionSites;

% select regression parameters

regressparams=[1]; % which trained hyperparameters to use
regresssets=ones(size(regressparams))*1; % which data set to use with each
clear regresslabels;
for i=1:length(regresssets)
    regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
end

dosldecomp = 0; % make sldecomp plots for each site
doMapField = 1; % make maps of the field

runPredictCESL;
