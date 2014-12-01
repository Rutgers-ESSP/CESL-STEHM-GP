% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Nov 21 15:13:27 EST 2014

% first do a table with rates for each of the models

firstyears=[0    0   400 800  1200 1600 1800 1900];
lastyears= [1800 400 800 1200 1600 1800 1900 2000];

fid=fopen('bymodel_GSLrates.tex','w');
fprintf(fid,'Model');

wfly=[firstyears ; lastyears];
fprintf(fid,' & %0.0f--%0.0f',wfly);
fprintf(fid,' \\\\ \n');

for iii=1:length(regresssets)
    datsub=find(testreg==0); sitesub=find(testsites==0);
    [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(f2s{iii}(datsub,1),V2s{iii}(datsub,datsub,1),testsites(sitesub),testreg(datsub),testX(datsub,3),firstyears,lastyears);

    fprintf(fid,trainlabels{regressparams(iii)});
    fprintf(fid,' & $%0.2f \\pm %0.2f$',[fslopeavg(1) 2*sdslopeavg(1)]);
    P=normcdf([fslopeavg(1)./sdslopeavg(1)]);
    if abs(P-.5)>=.49
        fprintf(fid,'***');
    elseif abs(P-.5)>=.45
        fprintf(fid,'**');
    elseif abs(P-.5)>=.4
        fprintf(fid,'*');
    elseif abs(P-.5)>=.133
        fprintf(fid,'$^\\dagger$');
    end
    
    sub=find(diffless==1);
    for qqq=sub(:)'
        fprintf(fid,' & $%0.1f \\pm %0.1f$',[fslopeavgdiff(qqq) ; 2*sdslopeavgdiff(qqq)]);
        P=normcdf([fslopeavgdiff(qqq)./sdslopeavgdiff(qqq)]);
        if abs(P-.5)>=.49
            fprintf(fid,'***');
        elseif abs(P-.5)>=.45
            fprintf(fid,'**');
        elseif abs(P-.5)>=.4
            fprintf(fid,'*');
        elseif abs(P-.5)>=.167
            fprintf(fid,'$^\\dagger$');
        end
    end
    
    fprintf(fid,' \\\\ \n');
end
fprintf(fid,'All rates except 0--1800 CE are detrended with respect to 0--1800 CE.');

fclose(fid);

% now do a table with the hyperparameters
fid=fopen('bymodel_theta.tex','w');

fprintf(fid,'Model & log L & $\\sigma_g$ & $\\tau_g$ & $\\sigma_l$ & $\\lambda_l$ & $\\sigma_m$ & $\\tau_m$ & $\\lambda_m$ & $\\sigma_w$ & $\\sigma_0$ & $\\sigma_c$ \\\\ \n');
for qqq=1:length(trainspecs)
    fprintf(fid,trainlabels{qqq});
    fprintf(fid,' & %0.0f',logp(qqq));
    fprintf(fid,' & %0.1f',thetTGG{qqq}(1:end-1));
    fprintf(fid,' & %0.2f',thetTGG{qqq}(end));
    fprintf(fid,' \\\\ \n');
end
fclose(fid);

% now do a table with site-specific rates

iii=1;

firstyears=[0    0   400 800  1200 1600 1800 1900];
lastyears= [1800 400 800 1200 1600 1800 1900 2000];

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
    qqq=1;
linra
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
    
    
    sub=find(diffless==1);
    for qqq=sub(:)'
        if ~isnan(fslopeavgdiff(kk,qqq))
            fprintf(fid,' & $%0.1f \\pm %0.1f$',[fslopeavgdiff(kk,qqq) ; 2*sdslopeavgdiff(kk,qqq)]);
            P=normcdf([fslopeavgdiff(kk,qqq)./sdslopeavgdiff(kk,qqq)]);
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
fprintf(fid,'All rates except 0--1800 CE are detrended with respect to 0--1800 CE.');

fclose(fid);


% now do a table with site sensitivities

iii=1;
firstyears=sitesensfirstyears;
lastyears=sitesenslastyears;

fid=fopen('sitesens_GSL.tex','w');
fid2=fopen('sitesens_theta.tex','w');

fprintf(fid,'Subset');

wfly=[firstyears ; lastyears];
fprintf(fid,' & %0.0f--%0.0f',wfly);
fprintf(fid,' \\\\ \n');

fprintf(fid2,'Subset & $\\sigma_g$ & $\\tau_g$ & $\\sigma_l$ & $\\lambda_l$ & $\\sigma_m$ & $\\tau_m$ & $\\lambda_m$ & $\\sigma_w$ & $\\sigma_0$ & $\\sigma_c$ \\\\ \n');

for reoptimize=0:1
    for qqq=1:size(sitesets,1)
        for rrr=1:size(sitesets,2)
            fprintf(fid1,sitesets{qqq,rrr});
            fprintf(fid2,sitesets{qqq,rrr});
            if reoptimize
                fprintf(fid,'*');
                fprintf(fid2,'*');
            end

            wfslope=sitesensfslope{iii,qqq,rrr,reoptimize+1};
            wsdslope=sitesenssdslope{iii,qqq,rrr,reoptimize+1};
            wfslopediff=sitesensfslopediff{iii,qqq,rrr,reoptimize+1};
            wsdslopediff=sitesenssdslopediff{iii,qqq,rrr,reoptimize+1};
            
            fprintf(fid,' & $%0.2f \\pm %0.2f$',[wfslope(1) 2*wsdslope(1)]);
            P=normcdf([wfslope(1)./wsdslope(1)]);
            if abs(P-.5)>=.49
                fprintf(fid,'***');
            elseif abs(P-.5)>=.45
                fprintf(fid,'**');
            elseif abs(P-.5)>=.4
                fprintf(fid,'*');
            elseif abs(P-.5)>=.133
                fprintf(fid,'$^\\dagger$');
            end
            
            sub=find(diffless==1);
            for uuu=sub(:)'
                fprintf(fid,' & $%0.1f \\pm %0.1f$',[wfslopediff(uuu) ; 2*wsdslopediff(uuu)]);
                P=normcdf([wfslopediff(uuu)./wsdslopediff(uuu)]);
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
            
            fprintf(fid,' \\\\ \n');
            
            dothet=sitesensthet{iii,qqq,rrr,reoptimize+1};
            
            fprintf(fid2,' & %0.1f',dothet(1:end-1));
            fprintf(fid2,' & %0.2f',dothet(end));
            fprintf(fid2,' \\\\ \n');     

        end
    end
end

fprintf(fid,'All rates except 0--1800 CE are detrended with respect to 0--1800 CE.');
fclose(fid);
fclose(fid2);