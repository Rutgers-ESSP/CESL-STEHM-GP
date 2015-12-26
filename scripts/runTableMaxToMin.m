% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Dec 26 12:19:07 EST 2015

clear mxtomnq;
qlevsmm=[.05 .167 .5 .833 .95];
for iii=1:length(regresssets)
   %calculate distribution minimum-to-maximum value between 0 and 1900 CE
    sub=find(testreg==0);
    samps=lhsnorm(f2s{iii}(sub,1),V2s{iii}(sub,sub,1),1000);
    sub2=find((testt<1900).*(testt>=0));
    mxtomn=max(samps(:,sub2),[],2)-min(samps(:,sub2),[],2);
    mxtomnq(iii,:) = quantile(mxtomn,qlevsmm);
end

fid=fopen('maxtomin.tsv','w');
fprintf(fid,'Maximum-To-Minimum, 0-1900 CE (mm)\n\n');
fprintf(fid,'\t%0.1f',qlevsmm*100);
fprintf(fid,'\n');
for iii=1:length(regresssets)
    fprintf(fid,labls{iii});
    fprintf(fid,'\t%0.0f',mxtomnq(iii,:));
    fprintf(fid,'\n');
end
fclose(fid);
    
