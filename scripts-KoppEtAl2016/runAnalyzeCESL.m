% Master sea-level rise estimation script for 
%
%	Kopp, R. E., A. C. Kemp, K. Bittermann, J. P. Donnelly, W. R. Gehrels, 
%	C. C. Hay, J. X. Mitrovica, E. D. Morrow, S. Rahmstorf, and B. P. Horton 
%	(2016). Temperature-driven global sea level variability in the Common Era.
%	Proceedings of the National Academy of Sciences.
%	doi: 10.1073/pnas.1517056113.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Jan 21 16:40:45 EST 2016

dosldecomp = 0; % make plots for each site?

% set up  paths

pd=pwd;
addpath([pd '/MFILES']);
addpath([pd '/scripts-KoppEtAl2016']);
IFILES=[pd '/IFILES-KoppEtAl2016'];
addpath(pd)
savefile='~/tmp/CESL';

WORKDIR='~/tmp/workdir-151230';
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);

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

testt = [-1000:20:1800 1810:10:2010]; % ages for regression
runSelectPredictionSites;

% select regression parameters

regressparams=[1 4 5 2 3]; % which trained hyperparameters to use
regresssets=[2 2 2 2 2]; % which data set to use with each
clear regresslabels;
for i=1:length(regresssets)
    regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
end


runCalculatePriorRates;
runTablePriors;
runLatexTablePriors;

doMapField = 1;
runPredictCESL

doMapField=0;
runScaledAmpSensitivityTest

