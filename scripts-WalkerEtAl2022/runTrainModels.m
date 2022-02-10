% Optimize hyperparameters for different model structures.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-06-15 19:09:51 -0400

trainrange=[100 100 100]; % optimize only using data with age errors < 100 yrs

clear thetTGG thethist trainsubsubset logp;

fid=fopen('thetTGG.tsv','a');
fprintf(fid,'set\ttraining data\tmodel\tlogp\tN\n');


for ii=1:length(trainspecs)

    disp(trainlabels{ii});
    ms = modelspec(trainspecs(ii));
    
    % first only fit ones without a compaction correction
    [thetTGG{ii},trainsubsubset{ii},logp(ii),thethist{ii}]= ...
        OptimizeHoloceneCovariance(datasets{trainsets(ii)}, ...
                                   ms,[3.4 3.0],trainfirsttime,trainrange,.01);

    % now add compaction correction factor
%    ms = modelspec(trainspecs(ii));
%    ms.thet0 = thetTGG{ii}(1:end-1);
%    ms.subfixed = 1:length(ms.thet0);
%    [thetTGG{ii},trainsubsubset{ii},logp(ii),thist]= ...
%        OptimizeHoloceneCovariance(datasets{trainsets(ii)}, ...
%                                   ms,[3.4 3.0],trainfirsttime(end),trainrange(end),1e6);   
%    thethist{ii}=[thethist{ii}; thist];

    
    % now final local optimization
    ms = modelspec(trainspecs(ii));
    ms.thet0 = thetTGG{ii}(1:end-1);
    startcompact = thetTGG{ii}(end);
    [thetTGG{ii},trainsubsubset{ii},logp(ii),thist]= ...
        OptimizeHoloceneCovariance(datasets{trainsets(ii)}, ...
                                   ms,[3.0],trainfirsttime(end),trainrange(end),1e6,startcompact);   
    thethist{ii}=[thethist{ii}; thist];

    
    for iii=1:length(thetTGG)
        fprintf(fid,[trainlabels{iii} '\t' datasets{trainsets(iii)}.label '\t' ...
                     modelspec(trainspecs(iii)).label]);
        fprintf(fid,['\t(%0.2f)'],logp(iii));
        fprintf(fid,'\t%0.0f',length(trainsubsubset{iii}));
        fprintf(fid,'\t%0.3f',thetTGG{iii});
        fprintf(fid,'\n');
    end
end

fclose(fid);
