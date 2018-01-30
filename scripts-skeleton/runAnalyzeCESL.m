% Skeleton sea-level rise estimation script for CESL-STEHM-GP 
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-01-29 19:00:41 -0500

% set up  paths

corepath = '~/Dropbox/Code/CESL-STEHM-GP';
coreIFILES=[corepath '/IFILES-working']; % point to the directory with the data files that aren't particular to this site, like the tide-gauge data set
addpath([corepath '/MFILES']); % add the path with the core MFILES

pd=pwd; % script expects to start in a directory that contains MFILES and IFILES as a subdir
addpath([pd '/scripts-KoppEtAl2016']); % add the path with the script files
IFILES=[pd '/IFILES-KoppEtAl2016']; % point to the directory with the data files
addpath(pd)
savefile='~/tmp/CESL-KoppEtAl2016'; % point to the .mat file you want the analysis to be backed up into 

WORKDIR='~/tmp/workdir-KoppEtAl2016-151230'; % point to the working directory you want your tables and figures dumped from
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);

PXdatafile=fullfile(IFILES,'RSL_All_19Feb2016.csv');

%%%%
% read in files and setup covariance model

runImportCESLDataSets; % this script handles importing data files into the expected structure
runSetupCESLCovariance; % this script sets up the prior covariance structure

%%%%
% map sites
runMapSites;

%%%%

% set up and run hyperparameter optimization

trainspecs=[1]; % identify the different model specifications (as defined in runSetupCESLCovariance) 
                          % that will be trained
trainsets = [2]; % identify the different datasets (in the structure datasets created by runImportCESLDataSeta)
                           % that will be used for each specification
trainfirsttime = -1000; % don't use data before trainfirsttime (default: -1000 CE) for training

trainlabels={};
for ii=1:length(trainsets)
    trainlabels = {trainlabels{:}, [datasets{trainsets(ii)}.label '_' modelspec(trainspecs(ii)).label]};
end

runTrainModels; % optimize the hyperparameters of the priors

save thetTGG thetTGG trainsubsubset
save(savefile);

thetTGG0=thetTGG;


%% now do a prediction
% define prediction sites

testt = [-1000:20:1800 1810:10:2010]; % ages for regression
runSelectPredictionSites;

% select regression parameters

regressparams=[1]; % which trained hyperparameters to use
regresssets=[1]; % which data set to use with each
clear regresslabels;
for i=1:length(regresssets)
    regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
end

dosldecomp = 0; % make plots for each site?
doMapField = 1; % make maps of the field?

runPredictCESL