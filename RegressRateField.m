function [fslopeavgF,sdslopeavgF,fsF,sdsF,fslopeavgdiffF,sdslopeavgdiffF,diffplusF,difflessF,passderivs,invcv] = RegressRateField(PX,modelspec,thetL,noiseMasks,Flat,Flong,firstyears,lastyears,trainsub,ICE5G,passderivs,invcv)

%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Nov 05 23:18:44 EST 2014

defval('firstyears',0);
defval('lastyears',1800);
defval('Flat',35:.5:52);
defval('Flong',232:.5:239);
defval('ICE5G',[]);
defval('noiseMasks',ones(length(thetL)));
defval('trainsub',[]);
defval('passderivs',[]);
defval('invcv',[]);

clear testsitedefF;
[FLONG,FLAT]=meshgrid(Flong,Flat);
Fsite=[1:length(FLAT(:))]';
allsites = [Fsite FLAT(:) FLONG(:)];

% partition into groups of 1000 points
lpt = length(unique(union(firstyears,lastyears)));
fsF=zeros(lpt*length(Fsite),1);
sdsF=zeros(lpt*length(Fsite),1);
VsF=sparse(lpt*length(Fsite),lpt*length(Fsite));
testlocsites=[];
testlocreg=[];
testloct=[];
counter=1;
for qqq=1:1000:length(Fsite)
    dosub=qqq:min(qqq+999,length(Fsite));
    dosub2=counter:(counter-1+length(dosub)*lpt);
    counter=dosub2(end)+1;
    clear testsitedefF;
    testsitedefF.sites=allsites(dosub,:);
    for ii=1:size(testsitedefF.sites,1)
        testsitedefF.names2=[num2str(testsitedefF.sites(ii,1)) '-' num2str(testsitedefF.sites(ii,2)) '-' num2str(testsitedefF.sites(ii,3))];
        testsitedefF.names=testsitedefF.names2;
    end

    if length(ICE5G)>0
        ux=mod(ICE5G.long,360);
        [sux,suxi]=sort(ux);
        testsitedefF.GIA = interp2(ICE5G.lat,sux,ICE5G.gia(suxi,:),testsitedefF.sites(:,2),mod(testsitedefF.sites(:,3),360),'linear');
        testsitedefF.GIA(find(testsitedefF.sites(:,2)>100))=0;
    end

    testtF = union(firstyears,lastyears);
    [fsF(dosub2),sdsF(dosub2),VsF(dosub2,dosub2),testlocsF,~,passderivs,invcv]=RegressHoloceneDataSets(PX,testsitedefF,modelspec,thetL,trainsub,noiseMasks,testtF,[],[],passderivs,invcv);
    testlocsites=[testlocsites ; testlocsF.sites];
    testlocreg=[testlocreg ; testlocsF.reg];
    testloct=[testloct ; testlocsF.X(:,3)];
end

[fslopeavgF,sdslopeavgF,fslopeavgdiffF,sdslopeavgdiffF,diffplusF,difflessF]=SLRateCompare(fsF,full(VsF),testlocsites,testlocreg,testloct,firstyears,lastyears);

