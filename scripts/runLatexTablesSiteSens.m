% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Mar 10 13:22:27 EDT 2015


% now do a table with site sensitivities

firstyears=sitesensfirstyears;
lastyears=sitesenslastyears;

fid=fopen(['sitesens_GSL' labl '.tex'],'w');
fid2=fopen(['sitesens_theta' labl '.tex'],'w');

fprintf(fid,'Subset');

wfly=[firstyears(2:end) ; lastyears(2:end)];
fprintf(fid,' & %0.0f--%0.0f',wfly);
fprintf(fid,' \\\\ \n');

fprintf(fid2,'Subset & $\\sigma_g$ & $\\tau_g$ & $\\sigma_l$ & $\\lambda_l$ & $\\sigma_m$ & $\\tau_m$ & $\\lambda_m$ & $\\sigma_w$ & $\\sigma_0$ & $\\sigma_{g0}$ & $\\sigma_c$ \\\\ \n');

for reoptimize=0
    for qqq=1:size(sitesets,1)
        for rrr=1:size(sitesets,2)
            fprintf(fid,sitesets{qqq,rrr});
            fprintf(fid2,sitesets{qqq,rrr});
            if reoptimize
                fprintf(fid,'*');
                fprintf(fid2,'*');
            end

            wfslope=sitesensfslope{iii,qqq,rrr,reoptimize+1};
            wsdslope=sitesenssdslope{iii,qqq,rrr,reoptimize+1};
            wfslopediff=sitesensfslopediff{iii,qqq,rrr,reoptimize+1};
            wsdslopediff=sitesenssdslopediff{iii,qqq,rrr,reoptimize+1};
            
           for uuu=2:length(firstyears)
             fprintf(fid,' & $%0.2f \\pm %0.2f$',[wfslope(uuu) 2*wsdslope(uuu)]);
            P=normcdf([wfslope(uuu)./wsdslope(uuu)]);
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

fclose(fid);
fclose(fid2);