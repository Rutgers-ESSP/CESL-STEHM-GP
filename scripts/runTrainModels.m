% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Nov 06 15:05:42 EST 2014

trainrange=[100 100 2000 2000];
clear thetTGG thethist trainsubsubset logp;
for ii=1:length(trainspecs)
    % first only fit ones without a compaction correction
    [thetTGG{ii},trainsubsubset{ii},logp(ii),thethist{ii}]= ...
        OptimizeHoloceneCovariance(datasets{trainsets(ii)}, ...
                                   modelspec(trainspecs(ii)),[2.4 3.4 3.4 3.0],trainfirsttime,trainrange,.01);

    % now add compaction correction factor
    ms = modelspec(trainspecs(ii));
    ms.thet0 = thetTGG{ii}(1:end-1);
    ms.subfixed = 1:length(ms.thet0);
    [thetTGG{ii},trainsubsubset{ii},logp(ii),thist]= ...
        OptimizeHoloceneCovariance(datasets{trainsets(ii)}, ...
                                   ms,[3.4 3.0],trainfirsttime(end),trainrange(end),1e6);   
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