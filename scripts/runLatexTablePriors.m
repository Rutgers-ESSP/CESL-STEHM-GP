% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Jan 04 16:44:47 EST 2016

% Table of Prior Rates

fid=fopen('priorrates.tex','w');

fprintf(fid,'Model');

wfly=[firstyearspriors(1:end) ; lastyearspriors(1:end)];
fprintf(fid,' & %0.0f--%0.0f',wfly);
fprintf(fid,' \\\\ \n');

for qqq=1:length(thetTGG)
    fprintf(fid,[trainlabels{qqq}]);
    for pp=1:length(firstyearspriors)
        fprintf(fid,' & $%0.2f \\pm %0.2f$',[priorslopef(qqq,pp) 2*priorslopesd(qqq,pp)]);
    end
    fprintf(fid,' \\\\ \n');
end

thetorder=[4 1 5 3 2];
fprintf(fid,'\n\n\n');
for qqq=thetorder
    fprintf(fid,['& ' trainlabels{qqq}]);
end
for pp=1:length(firstyearspriors)
    fprintf(fid,'\\\\ \n %0.0f--%0.0f',wfly(:,pp));
    for qqq=thetorder
       fprintf(fid,' & $%0.2f \\pm %0.2f$',[priorslopef(qqq,pp) 2*priorslopesd(qqq,pp)]);
    end
end

ia50=find(qlevsmmpriors==0.5);
[jk,ib,ia] = intersect([0.05 0.95],qlevsmmpriors);

fprintf(fid,'\\\\ \n 0--1900 amplitude');
for qqq=thetorder
    fprintf(fid,' & $\\pm %0.0f$ ($%0.0f$--$%0.0f$)',priormxtomnq(qqq,[ia50(:)' ia(:)'])/2/10);
end
fclose(fid);



%%%%%

% Table of Prior Amplitudes

ia50=find(qlevsmmpriors==0.5);
[jk,ib,ia] = intersect([0.05 0.95],qlevsmmpriors);
   
fid=fopen('priorvariabilityamplitude_0_1900.tex','w');
fprintf(fid,'Model & \\\\ \n');
for qqq=1:length(thetTGG)
    fprintf(fid,trainlabels{qqq});
    fprintf(fid,' & %0.0f (%0.0f--%0.0f)',priormxtomnq(qqq,[ia50(:)' ia(:)'])/2/10);
    fprintf(fid,' \\\\ \n');
end
fclose(fid);