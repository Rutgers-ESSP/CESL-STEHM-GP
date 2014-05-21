function [hp,hl,hl2,df,dsd,dV,outtable,difftimes,diffreg]=PlotPSLOverlay(years,regions,selregions,slf,slV,colrs,starttimes,endtimes,do2s,difftimestep,regionlabels)

% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed May 21 01:09:42 EDT 2014

defval('regions',ones(size(years)));
defval('selregions',unique(regions));
defval('starttimes',repmat(-Inf,1,length(selregions)));
defval('endtimes',repmat(Inf,1,length(selregions)));
defval('colrs',repmat({'b'},size(slf,2)));
defval('do2s',0);
defval('difftimestep',0);
defval('havecov',0);
defval('regionlabels',[]);

defval('hl',[]); defval('df',[]); defval('dV',[]);
defval('dsd',[]); defval('hl2',[]);

if length(endtimes)==1
    endtimes=repmat(endtimes,1,length(selregions));
end

outtable=[];

if prod(size(slV))==length(slV)
    slsd=slV;
    havecov=0;
else
    slsd=sqrt(diag(slV));
    havecov=1;
end

hp(1)=subplot(2,1,1);

for i=1:length(selregions)
    if length(regionlabels)==0
        outtable = [outtable sprintf('\nRegion %0.0f',selregions(i)) sprintf('\nYear\tmm\t1s')];
    else
        outtable = [outtable '\n' regionlabels{i} sprintf('\nYear\tmm\t1s')];
    end
    sub=find((years>=starttimes(i)).*(years<=endtimes(i)).*(regions==selregions(i)));
    if length(sub)>0
        hl(i)=plot(years(sub),slf(sub),[colrs{i}],'linew',2); hold on;
        plot(years(sub),slf(sub)+slsd(sub),[colrs{i} '--']);
        plot(years(sub),slf(sub)-slsd(sub),[colrs{i} '--']);
        if do2s>-1
            plot(years(sub),slf(sub)+slsd(sub),[colrs{i} '--']);
            plot(years(sub),slf(sub)-slsd(sub),[colrs{i} '--']);
        end
        if do2s>0
            plot(years(sub),slf(sub)+2*slsd(sub),[colrs{i} ':']);
            plot(years(sub),slf(sub)-2*slsd(sub),[colrs{i} ':']);
        end
        
        for jj=1:length(sub)
            outtable=[outtable sprintf('\n%0.0f',years(sub(jj))) sprintf('\t%0.2f', [slf(sub(jj)) slsd(sub(jj))])];
        end
    end
    
end

xlabel('Year');
ylabel('mm');

if difftimestep>0

    hp(2)=subplot(2,1,2);
    
    Mdiff = bsxfun(@eq,years,years')-bsxfun(@eq,years,years'+difftimestep);
    Mdiff = Mdiff .* bsxfun(@eq,regions,regions');
    sub=find(sum(Mdiff,2)==0);
    Mdiff=Mdiff(sub,:);
    difftimes=bsxfun(@rdivide,abs(Mdiff)*years,sum(abs(Mdiff),2));;
    diffreg=bsxfun(@rdivide,abs(Mdiff)*regions,sum(abs(Mdiff),2));;
    Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*years);
    
    df = Mdiff*slf;
    if havecov
        dV = Mdiff*slV*Mdiff';
    else
        dV = Mdiff*diag(slsd.^2)*Mdiff';
    end
    dV=.5*(dV+dV');
    dsd=sqrt(diag(dV));

    for i=1:length(selregions)
        if length(regionlabels)==0
            outtable = [outtable sprintf('\nRegion %0.0f',selregions(i)) sprintf('\nYear\tmm/y\t1s\tP>0')];
    	else
            outtable = [outtable '\n' regionlabels{i} sprintf('\nYear\tmm/y\t1s\tP>0')];
     	end
        if length(sub)>0
            sub=find((difftimes>=starttimes(i)).*(difftimes<=endtimes(i)).*(diffreg==selregions(i)));
            if length(sub)>0
                hl2(i)=plot(difftimes(sub),df(sub),[colrs{i}],'linew',2); hold on;
                if do2s>-1
                    plot(difftimes(sub),df(sub)+dsd(sub),[colrs{i} '--']);
                    plot(difftimes(sub),df(sub)-dsd(sub),[colrs{i} '--']);
                end
                if do2s>0
                    plot(difftimes(sub),df(sub)+2*dsd(sub),[colrs{i} ':']);
                    plot(difftimes(sub),df(sub)-2*dsd(sub),[colrs{i} ':']);
                end
                
                for jj=1:length(sub)
                    outtable=[outtable sprintf('\n%0.0f',difftimes(sub(jj))) sprintf('\t%0.3f', [df(sub(jj)) dsd(sub(jj))]) sprintf('\t%0.3f',normcdf(df(sub(jj))./dsd(sub(jj))))];
                end
            end
            
        end
    end
    
    xlabel('Year');
    ylabel(['mm/y (' num2str(difftimestep) '-y avg)']); 


end