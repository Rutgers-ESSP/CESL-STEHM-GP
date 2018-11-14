% Last updated by Bob Kopp, 2018-09-25 14:02:54 -0400

% set up  paths

rootdir='~/Dropbox/Code/CESL-STEHM-GP';

pd=pwd; % script expects to start in a directory that contains MFILES and IFILES as a subdir
addpath([rootdir '/MFILES']); % add the path with the core MFILES
addpath([pd '/scripts']); % add the path with the script files
IFILES=[pd '/IFILES']; % point to the directory with the data files
IFILES2=[rootdir '/IFILES-working'];
addpath(pd)
savefile='~/tmp/CESL-NAtlantic'; % point to the .mat file you want the analysis to be backed up into 

WORKDIR=[pd '/workdir-180925-global']; % point to the working directory you want your tables and figures dumped from
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);

PXdatafile=fullfile(IFILES,'RSL_All_29Mar2018.tsv');
latlim=[-90 90]; longlim=[-180 180];

runImportCESLDataSets;
runMapSites;
runSetupCESLCovariance;

% set up and run hyperparameter optimization

%trainspecs=1:length(modelspec); % identify the different model specifications (as defined in runSetupCESLCovariance) 
                          % that will be trained
trainspecs=[9 2];
trainsets =ones(size(trainspecs))*2; % identify the different datasets (in the structure datasets created by runImportCESLDataSeta)
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

testt = [-2000:20:1800 1810:10:2010]; % ages for regression
runSelectPredictionSites;

% select regression parameters

regressparams=[1 2]; % which trained hyperparameters to use
regresssets=ones(size(regressparams))*2; % which data set to use with each
clear regresslabels;
for i=1:length(regresssets)
    regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
end

dosldecomp = 0; % make sldecomp plots for each site
doMapField = 1; % make maps of the field

runPredictCESL;
