function [fslopeavgF,sdslopeavgF,fsF,sdsF] = RegressRateField(PX,modelspec,thetL,noiseMasks,Flat,Flong,firstyears,lastyears)

%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Oct 28 13:17:53 EDT 2014

defval('firstyears',0);
defval('lastyears',1800);
defval('Flat',35:.5:52);
defval('Flong',232:.5:239);
defval('ICE5G',[]);
defval('noiseMasks',ones(length(thetL)));

clear testsitedefF;
[FLONG,FLAT]=meshgrid(Flong,Flat);
Fsite=[1:length(FLAT(:))]';
testsitedefF.sites=[Fsite FLAT(:) FLONG(:)];
for ii=1:size(testsitedefF.sites,1)
    testsitedefF.names2=[num2str(testsitedefF.sites(ii,1)) '-' num2str(testsitedefF.sites(ii,2)) '-' num2str(testsitedefF.sites(ii,3))];
    testsitedefF.names=testsitedefF.names2;
end

if length(ICE5G)>0
    testsitedefF.refGIA = interp2(ICE5G.lat,ICE5G.long,ICE5G.gia,testsitedefF.sites(:,2),testsitedefF.sites(:,3),'linear');
    testsitedefF.refGIA(find(testsitedefF.sites(:,2)>100))=0;
end


testtF = union(firstyears,lastyears);
[fsF,sdsF,VsF,testlocsF]=RegressHoloceneDataSets(PX,testsitedefF,modelspec,thetL,[],noiseMasks,testtF);
[fslopeavgF,sdslopeavgF]=SLRateCompare(fsF,VsF,testlocsF.sites,testlocsF.reg,testlocsF.X(:,3),firstyears,lastyears);

