% Output table of all data points, including model prediction at sites
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-06-15 19:13:27 -0400

% Get model prediction for sites

u=unique(wdataset.datid);
clear fp sdp;
clear testsitedefp;
subps=[];
for pp=1:length(u)
    subp=find(wdataset.datid==u(pp));
    subq=find(wdataset.siteid==u(pp));
    subq=subq(1);
    if length(subp)>0
        testtsp{pp}=wdataset.meantime(subp);
        testsitedefp.sites(pp,:)=[wdataset.siteid(subq) ...
                            wdataset.sitecoords(subq,:)];
        testsitedefp.names(pp)=wdataset.sitenames(subq);
        testsitedefp.names2(pp)=wdataset.sitenames(subq);
        testsitedefp.firstage(pp)=min(wdataset.meantime(subp));
        %        testsitedefp.GISfp(pp)=wdataset.siteGISfp(subq);
        testsitedefp.GIA(pp)=wdataset.siteGIA(subq);
    end
    subps=[subps ; subp];
end
[fp(subps),sdp(subps),~,testlocp,~,derivsp]=RegressHoloceneDataSets(wdataset,testsitedefp,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks(1,:),testtsp,refyear,collinear,passderivs,invcv);

% output table

fid=fopen(['TGandProxyData' labl '.tsv'],'w');
fprintf(fid,['Site\tID\ttime1\ttime2\tlimiting\tGIAproj\tY-GIAproj\tY\' ...
             'tdY\tcompactcorr\tistg\tlat\tlong\tsite GIA\tmean time\tage-induced error\tYposterior\tdYposterior\tneterrorratio\n']);
for i=1:size(wdataset.datid,1)
    subq=find(wdataset.siteid==wdataset.datid(i));
    compactcorr=full(wdataset.compactcorr);
    subq=find(wdataset.siteid==wdataset.datid(i));
    fprintf(fid,[wdataset.sitenames{subq} '\t']);
    fprintf(fid,'%d\t',wdataset.datid(i));
    fprintf(fid,'%0.1f\t',wdataset.time1(i));
    fprintf(fid,'%0.1f\t',wdataset.time2(i));
    fprintf(fid,'%d\t',wdataset.limiting(i));
    fprintf(fid,'%0.2f\t',wdataset.GIAproj(i));
    fprintf(fid,'%0.2f\t',wdataset.Y(i));
    fprintf(fid,'%0.2f\t',wdataset.Y0(i));
    fprintf(fid,'%0.2f\t',wdataset.dY(i));
    fprintf(fid,'%0.2f\t',compactcorr(i));
    fprintf(fid,'%d\t',wdataset.istg(i));
    fprintf(fid,'%0.2f\t',wdataset.lat(i));
    fprintf(fid,'%0.2f\t',wdataset.long(i));
    fprintf(fid,'%0.2f\t',wdataset.siteGIA(subq(1)));
    fprintf(fid,'%0.1f\t',wdataset.meantime(i));
    fprintf(fid,'%0.2f\t',sqrt(derivsp.dK(i)));
    fprintf(fid,'%0.2f\t',fp(i));
    fprintf(fid,'%0.2f\t',sdp(i)); fprintf(fid,'%0.2f\n',(wdataset.Y0(i)-fp(i))/sqrt(derivsp.dK(i)^2+wdataset.dY(i)^2+sdp(i)^2));
end
fclose(fid);
