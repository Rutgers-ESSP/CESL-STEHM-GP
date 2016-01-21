% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Dec 26 13:04:14 EST 2015

% Table of Prior Rates

fid=fopen('priorrates.tsv','w');
for pp=1:length(firstyearspriors)
    fprintf(fid,'\tRate (avg, %0.0f-%0.0f)\t2s',[firstyearspriors(pp) lastyearspriors(pp)]);
end
for qqq=1:length(thetTGG)
    fprintf(fid,'\n');
    fprintf(fid,trainlabels{qqq});
    for pp=1:length(firstyearspriors)
        fprintf(fid,'\t%0.2f',[priorslopef(qqq,pp) 2*priorslopesd(qqq,pp)]);
    end
    fprintf(fid,'\n');
end
fclose(fid);


%%%%%

% Table of Prior Amplitudes


fid=fopen('priormaxtomin.tsv','w');
fprintf(fid,'Maximum-To-Minimum, 0-1900 CE (mm)\n\n');
fprintf(fid,'\t%0.1f',qlevsmmpriors*100);
fprintf(fid,'\n');
for iii=1:length(thetTGG)
    fprintf(fid,trainlabels{iii});
    fprintf(fid,'\t%0.0f',priormxtomnq(iii,:));
    fprintf(fid,'\n');
end
fclose(fid);
    
