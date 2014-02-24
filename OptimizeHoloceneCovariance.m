function [thetTGG,trainsub]=OptimizeHoloceneCovariance(dataset,modelspec,optimizesteps,mintime)

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Feb 18 09:08:43 EST 2014

istg = dataset.istg;
meantime=dataset.meantime;
lat=dataset.lat;
long=dataset.long;
Y=dataset.Y;
Ycv0=dataset.Ycv0;
bedmsk=dataset.bedmsk;
limiting=dataset.limiting;
obsGISfp=dataset.obsGISfp;
compactcorr=dataset.compactcorr;
time1=dataset.time1;
time2=dataset.time2;
dY=dataset.dY;

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

defval('optimizesteps',[1.1 2.1]);
defval('mintime',-1000);

donetg=0;

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
            if optimizesteps(nnn)>=1.1
                doglobs=[0 1];
                doneglob=1;
            else
                doglobs=0;
            end
            for doglob=doglobs
                [thetTGG(subnotfixed),minima] = SLGPOptimize(Y(trainsub),@(x) traincvTGG(meantime(trainsub),meantime(trainsub),dt1t1,x*Mfixed+fixedvect,Ycv0(trainsub,trainsub),dy1y1,bedmsk(trainsub,trainsub),fp1fp1),thetTGG(subnotfixed),lbTGG(subnotfixed),ubTGG(subnotfixed),doglob);
                disp(sprintf('%0.3f ',thetTGG));
            end

            lbTGG(1:5) = thetTGG(1:5); % set lower bound of amplitudes and temporal scale
            donetg=1;
            
        end

     elseif floor(optimizesteps(nnn))==2

        % optimize ignoring geochronological uncertainty

        thetTGG = [thetTGG .1];
        ubTGG = [ubTGG 5];
        lbTGG = [lbTGG 0];

        if donetg
            subfixed=union(subfixed,sublength); % fix length scales at those determined from the tide gauge data
        end
        subnotfixed=setdiff(1:length(thetTGG),subfixed);
        Mfixed=sparse(length(subnotfixed),length(thetTGG));
        for i=1:length(subnotfixed)
            Mfixed(i,subnotfixed(i))=1;
        end
        fixedvect = 0*thetTGG; fixedvect(subfixed)=thetTGG(subfixed);

        subnotfixed = [subnotfixed length(thetTGG)];
        trainsub = find((limiting==0));
        trainsub=intersect(trainsub,find(meantime>=mintime));
        %trainsub=intersect(trainsub,find(datid<1e6));
        dt1t1=dYears(meantime(trainsub),meantime(trainsub));
        dy1y1 = dDist([lat(trainsub) long(trainsub)],[lat(trainsub) long(trainsub)]);
        fp1fp1=bsxfun(@times,obsGISfp(trainsub)-1,obsGISfp(trainsub)'-1);

        if doneglob
            doglobs=0;
        else
            doglobs=[0 1];
        end
        for doglob=doglobs
            [thetTGG(subnotfixed),minima] = SLGPOptimize(Y(trainsub),@(x) traincvTGG(meantime(trainsub),meantime(trainsub),dt1t1,x(1:end-1)*Mfixed+fixedvect,Ycv0(trainsub,trainsub)+diag(x(end)*compactcorr(trainsub)).^2,dy1y1,bedmsk(trainsub,trainsub),fp1fp1),thetTGG(subnotfixed),lbTGG(subnotfixed),ubTGG(subnotfixed),doglob);
            disp(sprintf('%0.3f ',thetTGG));
        end
        
        if optimizesteps(nnn)>=2.1
            % now include geochronological uncertainty, one iteration

             wcvfunc = @(x1,x2,thet) cvfuncTGG(x1,x2,dYears(x1,x2),thet,dy1y1,bedmsk(trainsub,trainsub),fp1fp1);

            dt = abs(time2-time1)/4;

            [dK,df,d2f,yoffset] = GPRdx(meantime(trainsub),Y(trainsub),dt(trainsub),sqrt(dY(trainsub).^2+(thetTGG(end)*compactcorr(trainsub)).^2),@(x1,x2) wcvfunc(x1,x2,thetTGG),2);

            for doglob=[0]
                [thetTGG(subnotfixed),minima] = SLGPOptimize(Y(trainsub),@(x) traincvTGG(meantime(trainsub),meantime(trainsub),dt1t1,x(1:end-1)*Mfixed+fixedvect,Ycv0(trainsub,trainsub)+diag(x(end)*compactcorr(trainsub)).^2+diag(dK),dy1y1,bedmsk(trainsub,trainsub),fp1fp1),thetTGG(subnotfixed),lbTGG(subnotfixed),ubTGG(subnotfixed),doglob);
                disp(sprintf('%0.3f ',thetTGG));

            end
        end
    end

end

% thetTGG2 = [138.548 37.565 59.492 44.370 360.619 1.195 3.435 0.482 0.178 2.020 0.009 1.000 1.000 0.713 ]
