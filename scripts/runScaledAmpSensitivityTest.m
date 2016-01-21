% run a sensitivity test in which the amplitude hyperparameters are scaled by a factor of three
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Jan 21 13:27:02 EST 2016

if ~exist('scaledamps','dir')
    mkdir('scaledamps');
end

thetTGG00=thetTGG;
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

runCalculatePriorRates;
runTablePriors;
runLatexTablePriors;
runPredictCESL;
cd('..'); thetTGG=thetTGG00;

