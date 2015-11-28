% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Nov 28 15:36:57 EST 2015

trainrange=[100 100 100];
clear thetTGG thethist trainsubsubset logp;
for ii=1:length(trainspecs)

    ms = modelspec(trainspecs(ii));
    
    % first only fit ones without a compaction correction
    [thetTGG{ii},trainsubsubset{ii},logp(ii),thethist{ii}]= ...
        OptimizeHoloceneCovariance(datasets{trainsets(ii)}, ...
                                   modelspec(ms,[2.4 3.4 3.0],trainfirsttime,trainrange,.01);

    % now add compaction correction factor
    ms = modelspec(trainspecs(ii));
    ms.thet0 = thetTGG{ii}(1:end-1);
    ms.subfixed = 1:length(ms.thet0);
    [thetTGG{ii},trainsubsubset{ii},logp(ii),thist]= ...
        OptimizeHoloceneCovariance(datasets{trainsets(ii)}, ...
                                   ms,[3.4 3.0],trainfirsttime(end),trainrange(end),1e6);   
    thethist{ii}=[thethist{ii}; thist];

    
    % now final local optimization
    ms = modelspec(trainspecs(ii));
    ms.thet0 = thetTGG{ii}(1:end-1);
    startcompact = thetTGG{ii}(end);
    [thetTGG{ii},trainsubsubset{ii},logp(ii),thist]= ...
        OptimizeHoloceneCovariance(datasets{trainsets(ii)}, ...
                                   ms,[3.0],trainfirsttime(end),trainrange(end),1e6,startcompact);   
    thethist{ii}=[thethist{ii}; thist];

    
    fid=fopen('thetTGG.tsv','w');
    fprintf(fid,'set\ttraining data\tmodel\tlogp\tN\n');
    for iii=1:length(thetTGG)
        fprintf(fid,[trainlabels{iii} '\t' datasets{trainsets(iii)}.label '\t' ...
                     modelspec(trainspecs(iii)).label]);
        fprintf(fid,['\t(%0.2f)'],logp(iii));
        fprintf(fid,'\t%0.0f',length(trainsubsubset{iii}));
        fprintf(fid,'\t%0.3f',thetTGG{iii});
        fprintf(fid,'\n');
    end
    fclose(fid);
end