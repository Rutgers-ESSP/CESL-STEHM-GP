% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Mar 09 13:47:49 EDT 2015

% first do a table with rates for each of the models

regresssetorder = [2 1 3 5 4];

firstyears=[0 0   400 800  1200 1200 1600 1800 1860 1900];
lastyears= [1700 400 800 1200 1800 1600 1800 1900 1900 2000];
datsub=find(testreg==0); sitesub=find(testsites==0);

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
        P=normcdf([fslopeavg(qqq)./sdslopeavg(qqq)]);
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

    firstyears=[0    0   400 800  1200 1600 1800 1900];
    lastyears= [1700 400 800 1200 1600 1800 1900 2000];

    testreg=testlocs{iii}.reg;
    testsites=testlocs{iii}.sites;
    testX=testlocs{iii}.X;
    testnames2=testlocs{iii}.names2;

    [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(f2s{iii}(:,1),V2s{iii}(:,:,1),testsites,testreg,testX(:,3),firstyears,lastyears);

    fid=fopen(['linrates' labls{iii} '.tex'],'w');
    fprintf(fid,['Rates (mm/y), ' labls{iii} '\n \n']);
    fprintf(fid,'Site ');
    fprintf(fid,'& Ages');
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
        
        fprintf(fid,' & %0.0f',testsitedef.oldest(kk));
        fprintf(fid,'--%0.0f',testsitedef.youngest(kk));
        for qqq=1:size(fslopeavg,2)

            if ~isnan(fslopeavg(kk,qqq))
                fprintf(fid,' & $%0.2f \\pm %0.2f$',[fslopeavg(kk,qqq) 2*sdslopeavg(kk,qqq)]);
                P=normcdf([fslopeavg(kk,qqq)./sdslopeavg(kk,qqq)]);
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
[u,ui]=unique(datPX.textdata(:,3));
[u2,u2i]=unique(datPX.textdata(:,1));
u3i=union(ui,u2i);


fid=fopen('proxystudies.tex','w');
fprintf(fid,'Location & Mean Latitude & Mean Longitude & Median Age Range (CE) & N & Reference');

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
fclose(fid);