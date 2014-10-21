function [thetTGG,trainsub,logp,thethist]=OptimizeHoloceneCovariance(dataset,modelspec,optimizesteps,mintime,maxageerror,maxcompactcorrallowed,startcompact,maxcompactcorrfactor)

% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Oct 21 14:16:19 EDT 2014

defval('optimizesteps',[1.1 2.11]);
defval('mintime',-Inf);

defval('maxcompactcorrfactor',1);
defval('maxageerror',ones(size(optimizesteps))*Inf);
defval('maxcompactcorrallowed',Inf);
defval('startcompact',.1);

istg = dataset.istg;
lat=dataset.lat;
long=dataset.long;
Y=dataset.Y;
Ycv=dataset.Ycv;
limiting=dataset.limiting;
compactcorr=dataset.compactcorr;
time1=dataset.time1;
time2=dataset.time2;
dY=dataset.dY;
if isfield(dataset,'meantime')
    meantime=dataset.meantime;
else
    meantime=mean([dataset.time1(:) dataset.time2(:)],2);
end
if isfield(dataset,'obsGISfp')
    obsGISfp=dataset.obsGISfp;
else
    obsGISfp=ones(size(lat));
end

if (length(maxageerror)==1).*(length(maxageerror)<length(optimizesteps))
    maxageerror=repmat(maxageerror,1,length(optimizesteps));
end

if (length(maxcompactcorrallowed)==1).*(length(maxcompactcorrallowed)<length(optimizesteps))
    maxcompactcorrallowed=repmat(maxcompactcorrallowed,1,length(optimizesteps));
end

cvfuncTGG=modelspec.cvfunc;
traincvTGG=modelspec.traincv;
thetTGG=modelspec.thet0;
ubTGG=modelspec.ub;
lbTGG=modelspec.lb;
subfixed=modelspec.subfixed;
sublength=modelspec.sublength;
doneglob=0;

angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));

dYears=@(years1,years2) abs(bsxfun(@minus,years1',years2));
dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

donetg=0;
addedcompactcorr=0;

for nnn=1:length(optimizesteps)

    if floor(optimizesteps(nnn))==1
        % first optimize only tide gauge data

        subnotfixed=setdiff(1:length(thetTGG),subfixed);
        Mfixed=sparse(length(subnotfixed),length(thetTGG));
        for i=1:length(subnotfixed)
            Mfixed(i,subnotfixed(i))=1;
        end
        fixedvect = 0*thetTGG; fixedvect(subfixed)=thetTGG(subfixed);

        trainsub=find(istg==1);
        if length(trainsub)>0
            dt1t1=dYears(meantime(trainsub),meantime(trainsub));
            dy1y1 = dDist([lat(trainsub) long(trainsub)],[lat(trainsub) long(trainsub)]);
            fp1fp1=bsxfun(@times,obsGISfp(trainsub)-1,obsGISfp(trainsub)'-1);
            opttype = floor((optimizesteps(nnn)-1+1e-6)*10)/10;
            if (opttype==0.1)
                doglobs=[0 1];
                doneglob=1;
            elseif opttype==0.3
                doglobs=2;
                doneglob=1;
            elseif opttype==0.4
                doglobs=3;
                doneglob=1;
            else
                doglobs=0;
            end
            for doglob=doglobs
                [thetTGG(subnotfixed)] = SLGPOptimize(Y(trainsub),@(x) traincvTGG(meantime(trainsub),meantime(trainsub),dt1t1,x*Mfixed+fixedvect,Ycv(trainsub,trainsub),dy1y1,fp1fp1),thetTGG(subnotfixed),lbTGG(subnotfixed),ubTGG(subnotfixed),doglob);
                disp(sprintf('%0.3f ',thetTGG));
            end

            donetg=1;
            
        end

     elseif floor(optimizesteps(nnn))>=2

        % optimize ignoring geochronological uncertainty
%        if donetg
%            lbTGG(1:5) = thetTGG(1:5); % set lower bound of amplitudes and temporal scale
%        end
        if ~addedcompactcorr
            thetTGG = [thetTGG startcompact];
            ubTGG = [ubTGG maxcompactcorrfactor];
            lbTGG = [lbTGG 1e-6];
            addedcompactcorr = 1;
        end

        if donetg
            subfixed=union(subfixed,sublength); % fix length scales at those determined from the tide gauge data
        end
        subnotfixed=setdiff(1:length(thetTGG)-1,subfixed);
        Mfixed=sparse(length(subnotfixed),length(thetTGG)-1);
        for i=1:length(subnotfixed)
            Mfixed(i,subnotfixed(i))=1;
        end
        fixedvect = 0*thetTGG(1:end-1); fixedvect(subfixed)=thetTGG(subfixed);
        subnotfixed = unique([subnotfixed length(thetTGG)]);
        trainsub = find((limiting==0));
        trainsub=intersect(trainsub,find((meantime>=mintime).*(abs(time2-time1)<maxageerror(nnn)).*(abs(compactcorr)<=maxcompactcorrallowed(nnn))));
        %trainsub=intersect(trainsub,find(datid<1e6));
        
%        % remove extreme points
%        trainsub=intersect(trainsub,find(abs((Y-mean(Y))/std(Y))<3));
        
        dt1t1=dYears(meantime(trainsub),meantime(trainsub));
        dy1y1 = dDist([lat(trainsub) long(trainsub)],[lat(trainsub) long(trainsub)]);
        fp1fp1=bsxfun(@times,obsGISfp(trainsub)-1,obsGISfp(trainsub)'-1);

        % move points away from edges        
        dublb=ubTGG-lbTGG;
        subadj=find(abs(thetTGG(subnotfixed)-ubTGG(subnotfixed))<min(.01,.01*dublb(subnotfixed)));
        if length(subadj)>0
            thetTGG(subnotfixed(subadj))=thetTGG(subnotfixed(subadj))-min(.01,.01*dublb(subnotfixed(subadj)));
        end
        subadj=find(abs(thetTGG(subnotfixed)-lbTGG(subnotfixed))<min(.01,.01*dublb(subnotfixed)));
        if length(subadj)>0
            thetTGG(subnotfixed(subadj))=thetTGG(subnotfixed(subadj))+min(.01,.01*dublb(subnotfixed(subadj)));
        end
        
        opttype = floor((optimizesteps(nnn)-2+1e-6)*10)/10;
        opttype = round((opttype-floor(opttype))*10)/10;
        
        if opttype==0.1
            if doneglob
                doglobs=0;
            else
                doglobs=[0 1];
            end
            doneglob=1;
        elseif opttype==0.2
            doglobs=[0 1];
            doneglob=1;
        elseif opttype==0.3
            doglobs=2;
            doneglob=1;
       elseif opttype==0.4
            doglobs=3;
            doneglob=1;
        elseif opttype==0.5
            doglobs=[];
        else
            doglobs=0;
        end

        if floor(optimizesteps(nnn))==3
            dt = abs(time2-time1)/4;
         for doglob=doglobs
             [thetTGG(subnotfixed),logp] = ...
                    SLNIGPOptimize(meantime(trainsub),Y(trainsub), ...
                                   dt(trainsub),dY(trainsub),@(thet) ...
                                   diag(thet(end)*compactcorr(trainsub)).^2, @(x1,x2,thet,xxx,yyy) cvfuncTGG(x1,x2,dYears(x1,x2),thet(1:end-1)*Mfixed+fixedvect,dy1y1(xxx,yyy)',fp1fp1(xxx,yyy)'),thetTGG(subnotfixed),lbTGG(subnotfixed),ubTGG(subnotfixed),doglob,[1:length(trainsub)]');
                disp(sprintf('%0.3f ',thetTGG));
            end
        else
            
            for doglob=doglobs
                [thetTGG(subnotfixed),logp] = SLGPOptimize(Y(trainsub),@(x) traincvTGG(meantime(trainsub),meantime(trainsub),dt1t1,x(1:end-1)*Mfixed+fixedvect,Ycv(trainsub,trainsub)+diag(x(end)*compactcorr(trainsub)).^2,dy1y1,fp1fp1),thetTGG(subnotfixed),lbTGG(subnotfixed),ubTGG(subnotfixed),doglob);
                disp(sprintf('%0.3f ',thetTGG));
            end

            if mod(optimizesteps(nnn),0.1)>0.001
                % now include geochronological uncertainty, one iteration

                wcvfunc = @(x1,x2,thet,xxx,yyy) cvfuncTGG(x1,x2,dYears(x1,x2),thet,dy1y1(xxx,yyy)',fp1fp1(xxx,yyy)');

            dt = abs(time2-time1)/4;
 

                [dK,df,d2f,yoffset] = GPRdx(meantime(trainsub),Y(trainsub),dt(trainsub),sqrt(dY(trainsub).^2+(thetTGG(end)*compactcorr(trainsub)).^2),@(x1,x2,r1,r2) wcvfunc(x1,x2,thetTGG,r1,r2),1,[1:length(trainsub)]');

                doglobs=0;
                if opttype==0.4
                    doglobs=3;
                end
                
                for doglob=doglobs
                    [thetTGG(subnotfixed),logp] = SLGPOptimize(Y(trainsub)-yoffset,@(x) traincvTGG(meantime(trainsub),meantime(trainsub),dt1t1,x(1:end-1)*Mfixed+fixedvect,Ycv(trainsub,trainsub)+diag(x(end)*compactcorr(trainsub)).^2+diag(dK),dy1y1,fp1fp1),thetTGG(subnotfixed),lbTGG(subnotfixed),ubTGG(subnotfixed),doglob);
                    disp(sprintf('%0.3f ',thetTGG));

                end
            end
        end
        
    end
    thethist(nnn,:)=thetTGG;
end

% thetTGG2 = [138.548 37.565 59.492 44.370 360.619 1.195 3.435 0.482 0.178 2.020 0.009 1.000 1.000 0.713 ]
