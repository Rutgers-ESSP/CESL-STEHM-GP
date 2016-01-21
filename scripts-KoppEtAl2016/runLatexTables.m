% Generate assorted latex tables.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Jan 04 17:00:00 EST 2016

% first do a table with rates for each of the models

regresssetorder = [2 1 3 5 4];

firstyears=[0    0   300 700  1000 1400 1800 1860 1900 0 700 1400 1600 1860];
lastyears= [1700 300 700 1000 1400 1800 1900 1900 2000 700 2000 1600 1800 1960];
datsub=find(testreg==0); sitesub=find(testsites==0);


clear mxtomnq2;
qlevsmm=[.5 .05 .95];
for iii=1:length(regresssets)
   %calculate distribution minimum-to-maximum value between 0 and 1900 CE
    sub=find(testreg==0);
    samps=lhsnorm(f2s{iii}(sub,1),V2s{iii}(sub,sub,1),1000);
    sub2=find((testt<1900).*(testt>=0));
    mxtomn=max(samps(:,sub2),[],2)-min(samps(:,sub2),[],2);
    mxtomnq2(iii,:) = quantile(mxtomn,qlevsmm);
end
    
fid=fopen('bymodel_GSLrates.tex','w');
fprintf(fid,'Model');

wfly=[firstyears ; lastyears];
fprintf(fid,' & %0.0f--%0.0f',wfly);
fprintf(fid,' \\\\ \n');

for iii=regresssetorder
    [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(f2s{iii}(datsub,1),V2s{iii}(datsub,datsub,1),testsites(sitesub),testreg(datsub),testX(datsub,3),firstyears,lastyears);

    fprintf(fid,trainlabels{regressparams(iii)});
    for qqq=1:length(firstyears)
        fprintf(fid,' & $%0.2f \\pm %0.2f$',[fslopeavg(qqq) 2*sdslopeavg(qqq)]);
        if qqq==1
        P=normcdf([fslopeavg(qqq)./sdslopeavg(qqq)]);
            else
        P=normcdf([fslopeavgdiff(qqq-1)./sdslopeavgdiff(qqq-1)]);
        end
        
        if abs(P-.5)>=.49
            fprintf(fid,'***');
        elseif abs(P-.5)>=.45
            fprintf(fid,'**');
        elseif abs(P-.5)>=.4
            fprintf(fid,'*');
        elseif abs(P-.5)>=.133
            fprintf(fid,'$^\\dagger$');
        end
    end
    
    
% $$$     sub=find(diffless==1);
% $$$     for qqq=sub(:)'
% $$$         fprintf(fid,' & $%0.1f \\pm %0.1f$',[fslopeavgdiff(qqq) ; 2*sdslopeavgdiff(qqq)]);
% $$$         P=normcdf([fslopeavgdiff(qqq)./sdslopeavgdiff(qqq)]);
% $$$         if abs(P-.5)>=.49
% $$$             fprintf(fid,'***');
% $$$         elseif abs(P-.5)>=.45
% $$$             fprintf(fid,'**');
% $$$         elseif abs(P-.5)>=.4
% $$$             fprintf(fid,'*');
% $$$         elseif abs(P-.5)>=.167
% $$$             fprintf(fid,'$^\\dagger$');
% $$$         end
% $$$     end
    
    fprintf(fid,' \\\\ \n');
end
%fprintf(fid,'All rates except 0--1700 CE are detrended with respect to 0--1700 CE.');

fprintf(fid,'\n\n\n');

for iii=regresssetorder
    fprintf(fid,['& ' trainlabels{regressparams(iii)}]);
end
for qqq=1:length(firstyears)
    fprintf(fid,'\\\\ \n %0.0f--%0.0f',wfly(:,qqq));
    for iii=regresssetorder
        [fslopeavg,sdslopeavg]=SLRateCompare(f2s{iii}(datsub,1),V2s{iii}(datsub,datsub,1),testsites(sitesub),testreg(datsub),testX(datsub,3),firstyears(qqq),lastyears(qqq));
        P=normcdf([fslopeavg./sdslopeavg]);
        if qqq>1
            [~,~,fslopeavgdiff,sdslopeavgdiff]=SLRateCompare(f2s{iii}(datsub,1),V2s{iii}(datsub,datsub,1),testsites(sitesub),testreg(datsub),testX(datsub,3),firstyears([1 qqq]),lastyears([1 qqq]));
            P=normcdf([fslopeavgdiff./sdslopeavgdiff]);        
        end
        
        fprintf(fid,' & $%0.2f \\pm %0.2f$ (%0.2f)',[fslopeavg 2*sdslopeavg P]);

    end
end
fprintf(fid,'\\\\ \n Amplitude (0--1900)');

for iii=regresssetorder
     fprintf(fid,' & $\\pm %0.0f$ (n%0.0f--%0.0f)',mxtomnq2(iii,:)/2/10);
end

fprintf(fid,'\nAll probabilities except 0--1700 CE are with respect to differences from 0--1700 CE.');


%%% no flattener
% $$$ fprintf(fid,'\n\n\n No flattener \n\n');
% $$$ 
% $$$ regresssetorder2 = 6:10;
% $$$ for iii=regresssetorder2
% $$$     fprintf(fid,['& ' trainlabels{regressparams(iii)}]);
% $$$ end
% $$$ for qqq=1:length(firstyears)
% $$$     fprintf(fid,'\\\\ \n %0.0f--%0.0f',wfly(:,qqq));
% $$$     for iii=regresssetorder2
% $$$         [fslopeavg,sdslopeavg]=SLRateCompare(f2s{iii}(datsub,1),V2s{iii}(datsub,datsub,1),testsites(sitesub),testreg(datsub),testX(datsub,3),firstyears(qqq),lastyears(qqq));       
% $$$         P=normcdf([fslopeavg./sdslopeavg]);
% $$$         if qqq>1
% $$$             [~,~,fslopeavgdiff,sdslopeavgdiff]=SLRateCompare(f2s{iii}(datsub,1),V2s{iii}(datsub,datsub,1),testsites(sitesub),testreg(datsub),testX(datsub,3),firstyears([1 qqq]),lastyears([1 qqq]));
% $$$             P=normcdf([fslopeavgdiff./sdslopeavgdiff]);        
% $$$         end
% $$$         
% $$$         fprintf(fid,' & $%0.2f \\pm %0.2f$ (%0.2f)',[fslopeavg 2*sdslopeavg P]);
% $$$ 
% $$$     end
% $$$ end
fprintf(fid,'\nAll probabilities except 0--1700 CE are with respect to differences from 0--1700 CE.');

%fprintf(fid,'All rates except 0--1700 CE are detrended with respect to 0--1700 CE.');




fclose(fid);

% now do a table with the hyperparameters
fid=fopen('bymodel_theta.tex','w');

fprintf(fid,'Model & log L & $\\sigma_g$ & $\\tau_g$ & $\\sigma_l$ & $\\lambda_l$ & $\\sigma_m$ & $\\tau_m$ & $\\lambda_m$ & $\\sigma_w$ & $\\sigma_0$ & $\\sigma_{g0}$ & $\\sigma_c$ \\\\ \n');
for qqq=1:length(trainspecs)
    fprintf(fid,trainlabels{qqq});
    fprintf(fid,' & %0.0f',logp(qqq));
    fprintf(fid,' & %0.1f',thetTGG{qqq}(1:end-1));
    fprintf(fid,' & %0.2f',thetTGG{qqq}(end));
    fprintf(fid,' \\\\ \n');
end
fclose(fid);

% now do a table with site-specific rates

for iii=regresssetorder

    firstyears=[0    0   700 1400  1800 1900];
    lastyears= [1700 700 1400 1800 1900 2000];

    testreg=testlocs{iii}.reg;
    testsites=testlocs{iii}.sites;
    testX=testlocs{iii}.X;
    testnames2=testlocs{iii}.names2;

    [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(f2s{iii}(:,1),V2s{iii}(:,:,1),testsites,testreg,testX(:,3),firstyears,lastyears);

    fid=fopen(['linrates' labls{iii} '.tex'],'w');
    fprintf(fid,['Rates (mm/y), ' labls{iii} '\n \n']);
    fprintf(fid,'Site ');
    %fprintf(fid,'& Ages');
    %fprintf(fid,' & Lat & Long & Oldest & Youngest');
    for pp=1:length(firstyears)
        fprintf(fid,' & %0.0f--%0.0f',[firstyears(pp) lastyears(pp)]);
    end
    fprintf(fid,' \\\\ \n');

    for kk=1:size(testsites,1)
        fprintf(fid,testnames2{kk});
        %    if testsitedef.sites(kk,2)<360
        %    fprintf(fid,' & %0.2f',testsitedef.sites(kk,2:3));
        %else
        %    fprintf(fid,' & & ');
        %end
        
        %fprintf(fid,' & %0.0f',testsitedef.firstage(kk));
        %fprintf(fid,'--%0.0f',testsitedef.youngest(kk));
        for qqq=1:size(fslopeavg,2)

            if ~isnan(fslopeavg(kk,qqq))
                fprintf(fid,' & $%0.2f \\pm %0.2f$',[fslopeavg(kk,qqq) 2*sdslopeavg(kk,qqq)]);
                if qqq>1
                    P=normcdf([fslopeavgdiff(kk,qqq-1)./sdslopeavgdiff(kk,qqq-1)]);
                else
                    P=normcdf([fslopeavg(kk,qqq)./sdslopeavg(kk,qqq)]);
                end
                
                if abs(P-.5)>=.49
                    fprintf(fid,'***');
                elseif abs(P-.5)>=.45
                    fprintf(fid,'**');
                elseif abs(P-.5)>=.4
                    fprintf(fid,'*');
                elseif abs(P-.5)>=.167
                    fprintf(fid,'$^\\dagger$');
                end
            else
                fprintf(fid,' & ');
            end
        end
        
        fprintf(fid,'\\\\ \n');
    end

    fclose(fid);
end

% list of tide gauges

fid=fopen('tidegauges.tex','w');
fprintf(fid,'Name & Latitude & Longitude & Ages & PSMSL ID');

for kk=1:length(TGold.siteid)
    sub=find(TGold.datid==TGold.siteid(kk));
    fprintf(fid,'\\\\ \n');
    fprintf(fid,TGold.sitenames{kk});
    fprintf(fid,' & %0.2f',[TGold.sitecoords(kk,:)]);
    fprintf(fid,' & %0.0f',[min(TGold.meantime(sub)) max(TGold.meantime(sub))]);
    fprintf(fid,' & %0.0f',TGold.siteid(kk));
end

[u,ui]=sort(TG.siteid);
for kk=ui(:)'
    sub=find(TG0.datid==TG.siteid(kk));
    if length(sub)>0
        fprintf(fid,'\\\\ \n');
        fprintf(fid,TG.sitenames{kk});
        fprintf(fid,' & %0.2f',[TG.sitecoords(kk,:)]);
        fprintf(fid,' & %0.0f',[min(TG0.meantime(sub)) max(TG0.meantime(sub))]);
        fprintf(fid,' & %0.0f',TG.siteid(kk));
    end
    
end
fclose(fid);

% now proxy sites


fid=fopen('proxystudies.tex','w');
fprintf(fid,'Location & Mean Latitude & Mean Longitude & Median Age Range (CE) & N & Reference');

[u,ui]=unique(datPX.textdata(:,3));
[u2,u2i]=unique(datPX.textdata(:,1));
u3i=union(ui,u2i);

for kk=1:length(u3i)
    wname = datPX.textdata(u3i(kk),1); wcite = datPX.textdata(u3i(kk),3);

    sub=find(strcmpi(wcite,datPX.textdata(:,3)));
    sub=intersect(sub,find(strcmpi(wname,datPX.textdata(:,1))));
    wmediantime=datPX.data(sub,6);

    fprintf(fid,'\\\\ \n');
    fprintf(fid,wname{:});
    fprintf(fid,' & %0.1f',mean(datPX.data(sub,[1 2]),1));
    fprintf(fid,' & %0.0f to %0.0f',[min(wmediantime) max(wmediantime)]);
    fprintf(fid,' & %0.0f',length(sub));
    fprintf(fid,[' & ' wcite{:}]);
end

if exist('datLH')
    [u,ui]=unique(datLH.textdata(:,3));
    [u2,u2i]=unique(datLH.textdata(:,1));
    u3i=union(ui,u2i);

    for kk=1:length(u3i)
        wname = datLH.textdata(u3i(kk),1); wcite = datLH.textdata(u3i(kk),3);

        sub=find(strcmpi(wcite,datLH.textdata(:,3)));
        sub=intersect(sub,find(strcmpi(wname,datLH.textdata(:,1))));
        wmediantime=datLH.data(sub,6);

        fprintf(fid,'\\\\ \n');
        fprintf(fid,wname{:});
        fprintf(fid,' & %0.1f',mean(datLH.data(sub,[1 2]),1));
        fprintf(fid,' & %0.0f to %0.0f',[min(wmediantime) max(wmediantime)]);
        fprintf(fid,' & %0.0f',length(sub));
        fprintf(fid,[' & ' wcite{:}]);
    end
end

    fclose(fid);

    % now do a table with half max to min rates
    

clear mxtomnq2;
qlevsmm=[.5 .05 .95];
for iii=1:length(regresssets)
   %calculate distribution minimum-to-maximum value between 0 and 1900 CE
    sub=find(testreg==0);
    samps=lhsnorm(f2s{iii}(sub,1),V2s{iii}(sub,sub,1),1000);
    sub2=find((testt<1900).*(testt>=0));
    mxtomn=max(samps(:,sub2),[],2)-min(samps(:,sub2),[],2);
    mxtomnq2(iii,:) = quantile(mxtomn,qlevsmm);
end
    
fid=fopen('variabilityamplitude_0_1900.tex','w');

fprintf(fid,'Model & \\\\ \n');
for qqq=1:length(regresssets)
    fprintf(fid,labls{qqq});
    fprintf(fid,' & %0.0f (%0.0f--%0.0f)',mxtomnq2(qqq,:)/2/10);
    fprintf(fid,' \\\\ \n');
end
fclose(fid);