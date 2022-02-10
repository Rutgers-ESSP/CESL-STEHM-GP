% last updated by Bob Kopp, 2017-06-30 22:52:06 -0400

fid=fopen(['sldecomp' labl '.tsv'],'w');

fprintf(fid,'Site\tLat\tLong\tYear');
for ssss=1:length(noiseMasklabels);
    fprintf(fid,['\t' noiseMasklabels{ssss} ' (mm)\t1s']);
end
fprintf(fid,'\n');

for rrrr=1:size(f2s{iii},1)
    subname=find(testsites(:,1)==testreg(rrrr));
    fprintf(fid,testnames2{subname});
    fprintf(fid,'\t%0.2f',testX(rrrr,1:2));
    fprintf(fid,'\t%0.0f',testX(rrrr,3));
    for ssss=1:size(f2s{iii},2)
        fprintf(fid,'\t%0.2f',f2s{iii}(rrrr,ssss));
        fprintf(fid,'\t%0.2f',sd2s{iii}(rrrr,ssss));        
    end
    fprintf(fid,'\n');
end

fclose(fid);