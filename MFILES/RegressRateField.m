function [fslopeavgF,sdslopeavgF,fsF,sdsF,fslopeavgdiffF,sdslopeavgdiffF,diffplusF,difflessF,passderivs,invcv,testlocsFcum] = RegressRateField(PX,modelspec,thetL,noiseMasks,Flat,Flong,firstyears,lastyears,trainsub,ICE5G,passderivs,invcv,chunksize)

%  [fslopeavgF,sdslopeavgF,fsF,sdsF,fslopeavgdiffF,sdslopeavgdiffF,diffplusF,difflessF,
%       passderivs,invcv] = RegressRateField(PX,modelspec,thetL,[noiseMask],[Flat],[Flong],
%       firstyear,lastyear,[trainsub],[GIA],[passderivs],[invcv])
%
% Run RegressHoloceneDataSets to calculate sea level on the grid specified
% by Flat and Flong, from dataset PX with model specification modelspec,
% hyperparameters thetL, noise mask noiseMask (default = no noise mask),
% training subset trainsub (default = all data points), GIA rate structure GIA,
% and optionally derivative structure passderiv and inverse covariance invcv,
% then convert to a rate field based on the difference between firstyear and
% lastyear.
%
% EXAMPLE:
%
%   Flat=-85:10:85;
%   Flong=0:20:360;
%
%   [fslopeF,sdslopeF,~,~,~,~,~,~,passderivs,invcv] = ...
%      RegressRateField(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:), ...
%                       Flat,Flong,0,1700,trainsub,[]);    
%
%   Flat1=min(Flat):max(Flat);
%   Flong1=min(Flong):max(Flong);
% 
%   [FLONG,FLAT]=meshgrid(Flong,Flat);
%   [FLONG1,FLAT1]=meshgrid(Flong1,Flat1);
%   mapped = griddata(FLONG(:),FLAT(:),fslopeF,Flong1,Flat1(:),'linear');
%   sdmapped = griddata(FLONG(:),FLAT(:),sdslopeF,Flong1,Flat1(:),'linear');
%
%   clf;
%   ax = worldmap('World');
%   setm(ax, 'Origin',[0 -90 0],'meridianlabel','off','parallellabel','off','flinewidth',3);
%   land = shaperead('landareas', 'UseGeoCoords', true);
%   geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
%   hold on;
%   hs1=scatterm(FLAT1(:),FLONG1(:),10,mapped(:),'filled','marker','s');
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2019-11-27 11:08:29 -0500

defval('firstyears',0);
defval('lastyears',1800);
defval('Flat',35:.5:52);
defval('Flong',232:.5:239);
defval('ICE5G',[]);
defval('noiseMasks',ones(length(thetL)));
defval('trainsub',[]);
defval('passderivs',[]);
defval('invcv',[]);
defval('chunksize',1000)

clear testsitedefF;
[FLONG,FLAT]=meshgrid(Flong,Flat);
Fsite=[1:length(FLAT(:))]';
allsites = [Fsite FLAT(:) FLONG(:)];

% partition into groups of chunksize points
lpt = length(unique(union(firstyears,lastyears)));
fsF=zeros(lpt*length(Fsite),1);
sdsF=zeros(lpt*length(Fsite),1);
VsF=sparse(lpt*length(Fsite),lpt*length(Fsite));
testlocsites=[];
testlocreg=[];
testloct=[];
testlocX=[];
counter=1;
for qqq=1:chunksize:length(Fsite)
    dosub=qqq:min(qqq+(chunksize-1),length(Fsite));
    dosub2=counter:(counter-1+length(dosub)*lpt);
    counter=dosub2(end)+1;
    clear testsitedefF;
    testsitedefF.sites=allsites(dosub,:);
    for ii=1:size(testsitedefF.sites,1)
        testsitedefF.names2=[num2str(testsitedefF.sites(ii,1)) '-' num2str(testsitedefF.sites(ii,2)) '-' num2str(testsitedefF.sites(ii,3))];
        testsitedefF.names=testsitedefF.names2;
    end

    if length(ICE5G)>1
        ux=mod(ICE5G.long,360);
        [sux,suxi]=sort(ux);
        testsitedefF.GIA = interp2(ICE5G.lat,sux,ICE5G.gia(suxi,:),testsitedefF.sites(:,2),mod(testsitedefF.sites(:,3),360),'linear');
        testsitedefF.GIA(find(testsitedefF.sites(:,2)>100))=0;
    end

    testtF = union(firstyears,lastyears);

    collinear=[];
    if isfield(modelspec,'subamplinear')
        if length(modelspec.subamplinear)>0
            collinear=modelspec.subamplinear(1);
        end
    end
    [fsF(dosub2),sdsF(dosub2),VsF(dosub2,dosub2),testlocsF,~,passderivs,invcv]=RegressHoloceneDataSets(PX,testsitedefF,modelspec,thetL,trainsub,noiseMasks,testtF,[],collinear,passderivs,invcv);
    testlocsites=[testlocsites ; testlocsF.sites];
    testlocreg=[testlocreg ; testlocsF.reg];
    testloct=[testloct ; testlocsF.X(:,3)];
    testlocX=[testlocX ; testlocsF.X];
end

testlocsFcum.sites=testlocsites;
testlocsFcum.reg=testlocreg;
testlocsFcum.t=testloct;
testlocsFcum.X=testlocX;

[fslopeavgF,sdslopeavgF,fslopeavgdiffF,sdslopeavgdiffF,diffplusF,difflessF]=SLRateCompare(fsF,VsF,testlocsites,testlocreg,testloct,firstyears,lastyears);

