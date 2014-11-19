% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Nov 18 18:44:53 EST 2014

firstyears=[0    0   400 800  1200 1600 1800 1900];
lastyears= [1800 400 800 1200 1600 1800 1900 2000];

[fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(f2s{ii,jj}(:,1),V2s{ii,jj}(:,:,1),testsites,testreg,testX(:,3),firstyears,lastyears);
[fslopereg,sdslopereg,fsloperegdiff,sdsloperegdiff,diffplus,diffless]=SLRateCompare(f2s{ii,jj}(:,4),V2s{ii,jj}(:,:,4),testsites,testreg,testX(:,3),firstyears,lastyears);
[fslopelin,sdslopelin]=SLRateCompare(f2s{ii,jj}(:,3),V2s{ii,jj}(:,:,3),testsites,testreg,testX(:,3),testt(end-1),testt(end));

fid=fopen(['linrates' labl '.tsv'],'w');
fprintf(fid,['Regional+Local Linear Rates (mm/y), ' labl '\n']);
fprintf(fid,'Site\tSiteID\tLat\tLong\tOldest\tYoungest\tICE5G VM2-90\tRate (linear)\t2s');
for pp=1:length(firstyears)
    fprintf(fid,'\tRate (avg, %0.0f-%0.0f)\t2s',[firstyears(pp) lastyears(pp)]);
end
for pp=1:length(diffplus)
    fprintf(fid,['\tRate Diff. (avg, %0.0f-%0.0f minus ' ...
                 '%0.0f-%0.0f)\t2s'],[firstyears(diffplus(pp)) ...
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
    for pp=1:length(firstyears)
        fprintf(fid,'\t%0.2f',[fslopeavg(kk,pp) 2*sdslopeavg(kk,pp)]);
    end
    for pp=1:length(diffplus)
        fprintf(fid,'\t%0.2f',[fslopeavgdiff(kk,pp) 2*sdslopeavgdiff(kk,pp)]);
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n');
fprintf(fid,noiseMasklabels{4});
fprintf(fid,'\n');
for kk=1:size(testsites,1)
    fprintf(fid,testnames2{kk});
    fprintf(fid,'\t%0.2f',testsitedef.sites(kk,:));
    fprintf(fid,'\t%0.0f',testsitedef.oldest(kk));
    fprintf(fid,'\t%0.0f',testsitedef.youngest(kk));       
    fprintf(fid,'\t%0.2f',testsitedef.GIA(kk));
    fprintf(fid,'\t%0.2f',[fslopelin(kk) 2*sdslopelin(kk)]);
    for pp=1:length(firstyears)
        fprintf(fid,'\t%0.2f',[fslopereg(kk,pp) 2*sdslopereg(kk,pp)]);
    end
    for pp=1:length(diffplus)
        fprintf(fid,'\t%0.2f',[fsloperegdiff(kk,pp) 2*sdsloperegdiff(kk,pp)]);
    end
    fprintf(fid,'\n');
end

fclose(fid);