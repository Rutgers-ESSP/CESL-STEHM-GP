% let's check whether all 60-year periods beginning 1000 BCE and running through 1700 CE
% are slower than each 60-year period beginning after 1700 CE
% 2019-04-10 14:02:44 -0400

timestep=60;
timestep2=20;
backgroundyrs=[0 1700];
firstyears0=[backgroundyrs(1):timestep2:(2000-timestep)];
lastyears0= [(backgroundyrs(1)+timestep):timestep2:2000];

for kk=1:length(testsites(:,1))

    subSite=find((testreg==testsites(kk,1))); 
    sub2=find(firstyears0>=min(testX(subSite,3)));

    firstyears=firstyears0(sub2);
    lastyears=lastyears0(sub2);

    if length(firstyears)>1
        [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless,Vslopeavg,Vslopeavgdiff]=SLRateCompare(f2s{iii}(subSite,1),V2s{iii}(subSite,subSite,1),testsites(kk),testreg(subSite),testX(subSite,3),firstyears,lastyears);

        startyearplus=firstyears(diffplus);
        endyearplus=lastyears(diffplus);

        startyearminus=firstyears(diffless);
        endyearminus=lastyears(diffless);

        if min(startyearminus)<backgroundyrs(end)
            doyrsub=find(firstyears>=backgroundyrs(2));
            doyrstart=firstyears(doyrsub);
            doyrend=lastyears(doyrsub);

            clear Pmatrix avgP Pmindiffsamps Pallpos mindiffsamps;

            for kkk=1:length(doyrsub)
                doyr=firstyears(doyrsub(kkk));

                sub=find((startyearplus==doyr).*(endyearplus==(doyr+timestep)).*(startyearminus<backgroundyrs(2)));
                U = [fslopeavgdiff ; sdslopeavgdiff]';
                p = 1-normcdf(0,fslopeavgdiff,sdslopeavgdiff);
                Pmatrix(kkk,:)=p(sub);
                avgP(kkk) = mean(p(sub));

                % check prob minimum rate of subsequent periods is faster than random preceeding period
            
                sub=find((startyearplus>=doyr).*(startyearminus<backgroundyrs(2)));
                N = 10000;
                try
                    samps=mvnrnd(fslopeavgdiff(sub),Vslopeavgdiff(sub,sub),N); % samples of differnces 
                catch
                    samps=mvnrnd(fslopeavgdiff(sub),Vslopeavgdiff(sub,sub)+eye(length(sub))*1e-2,N); % samples of differnces 
                end

                u=unique(startyearminus(sub));
                clear mindiffsamps;
                for vvvv=1:length(u)
                    sub2=find(startyearminus(sub)==u(vvvv));
                    mindiffsamps(:,vvvv)=min(samps(:,sub2),[],2);
                end
                Pmindiffsamps(kkk)=mean(mean(mindiffsamps>0,2));
                Pallpos(kkk) = mean(min(mindiffsamps,[],2)>0);

            end

            compyearstart = startyearminus(sub);
            compyearend = endyearminus(sub);

            clf;
            plot(doyrstart,avgP); 
            hold on;
            plot(doyrstart,Pmindiffsamps,'g');
            plot(doyrstart,Pallpos,'r');
            title({[num2str(timestep) '-yr steps: ' testnames2{kk} ],sprintf('Bkgd = %0.0f-%0.0f',firstyears(1),backgroundyrs(2))});

            legend('P > rand bkgd','P all subs > rand bkgd','P all subs > all bkgd','Location','Southeast');
            xlabel('start year');
            fn=testnames2{kk}(setdiff(1:length(testnames2{kk}),strfind(testnames2{kk},' ')));
            pdfwrite(['P-' num2str(timestep) 'y-' fn]);

            fid=fopen(['P-' num2str(timestep) 'y-' fn '.tsv'],'w');
            fprintf(fid,[testnames2{kk} '\n']);
            fprintf(fid,'start year\tend year\trate (mm/yr)\trate (1s)\tP > rand bkgd\tP all subs > rand bkgd\tP all subs > all bkgd\n');
            M=[doyrstart ; doyrend ; fslopeavg(doyrsub) ; sdslopeavg(doyrsub) ; avgP ; Pmindiffsamps ; Pallpos ];
            fprintf(fid,'%0.0f\t%0.0f\t%0.2f\t%0.2f\t%0.3f\t%0.3f\t%0.3f\n',M);
            fclose(fid);
        end
    end
end
