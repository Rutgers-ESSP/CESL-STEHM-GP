% Generate table of rates at each site, and difference sin rate at each site
%
% start and end year of each interval of interest

firstyears=[0 1000 0 1000 1400 1000 1200 1400 1600 1800 0 0 1000 1700 0 1000 1700 1800 1900 1940 0 500 1300 1600];
lastyears= [2000 2000 1000 1400 1800 1200 1400 1600 1800 2000 1800 2000 2000 2000 1700 1700 1800 1900 2000 2000 500 1300 1600 1800];

% rate calculation
[fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(f2s{iii}(:,1),V2s{iii}(:,:,1),testsites,testreg,testX(:,3),firstyears,lastyears);

% linear rate term
[fslopelin,sdslopelin]=SLRateCompare(f2slin{iii}(:,1),V2slin{iii}(:,:,1),testlocslin{iii}.sites,testlocslin{iii}.reg,testlocslin{iii}.X(:,3),0,1800);
for kk=1:size(testsites,1)
    sub=find(testlocslin{iii}.reg==testsites(kk));
    if length(sub)>0
        sub=sub(1);
        foffset(kk) = f2slin{iii}(sub,2);
        sdoffset(kk) = sd2slin{iii}(sub,2);
    else
        foffset(kk)=NaN;
        sdoffset(kk)=NaN;
    end
end

% output table

fid=fopen(['linrates' labl '.tsv'],'w');
fprintf(fid,['Rates (mm/y), ' labl '\n']);
fprintf(fid,'Site\tSiteID\tLat\tLong\tOldest\tYoungest\tICE5G VM2-90\tRate (linear)\t2s\tOffset (mm)\t2s');
for pp=1:length(firstyears)
    fprintf(fid,'\tRate (avg, %0.0f-%0.0f)\t2s\tP>0',[firstyears(pp) lastyears(pp)]);
end
for pp=1:length(diffplus)
    fprintf(fid,['\tRate Diff. (avg, %0.0f-%0.0f minus ' ...
                    '%0.0f-%0.0f)\t2s\tP>0'],[firstyears(diffplus(pp)) ...
                        lastyears(diffplus(pp)) firstyears(diffless(pp)) lastyears(diffless(pp))]);
end

fprintf(fid,'\n');
fprintf(fid,noiseMasklabels{1});
fprintf(fid,'\n');
for kk=1:size(testsites,1)
    fprintf(fid,testnames2{kk});
    fprintf(fid,'\t%0.2f',testsitedef.sites(kk,:));
    fprintf(fid,'\t%0.0f',testsitedef.oldest(kk));
    fprintf(fid,'\t%0.0f',testsitedef.youngest(kk));       
    fprintf(fid,'\t%0.2f',testsitedef.GIA(kk));
    fprintf(fid,'\t%0.2f',[fslopelin(kk) 2*sdslopelin(kk)]);
    fprintf(fid,'\t%0.2f',[foffset(kk) 2*sdoffset(kk)]);
    for pp=1:length(firstyears)
        fprintf(fid,'\t%0.2f',[fslopeavg(kk,pp) 2*sdslopeavg(kk,pp)]);
        fprintf(fid,'\t%0.3f',[normcdf(fslopeavg(kk,pp)/sdslopeavg(kk,pp))]);
    end
    for pp=1:length(diffplus)
        fprintf(fid,'\t%0.2f',[fslopeavgdiff(kk,pp) 2*sdslopeavgdiff(kk,pp)]);
        fprintf(fid,'\t%0.3f',[normcdf(fslopeavgdiff(kk,pp)/sdslopeavgdiff(kk,pp))]);
    end
    fprintf(fid,'\n');
end

fclose(fid);