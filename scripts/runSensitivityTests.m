% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Nov 06 14:59:29 EST 2014

        %% do data sensitivity tests
     
        
        wtestsitedef.sites=[0 1e6 1e6];
        wtestsitedef.names={'GSL'};
        wtestsitedef.names2={'GSL'};
        wtestsitedef.firstage=min(oldest);

        wtestt=[1000 1800];
        Mdiff = [-1 1]/diff(wtestt);

        sitesub=(find(wdataset.siteid>9999));
        dosites=wdataset.siteid(sitesub);
        dositenames=wdataset.sitenames(sitesub);
        docoords=wdataset.sitecoords(sitesub,:);
        trainsub0=find((wdataset.limiting==0));
        
        clear sensf senssd sensnames senscoords;
        for qqq=0:length(dosites)
            for doexcl=0:1
                if qqq==0
                    if doexcl
                        trainsub=1;
                        wdataset2=wdataset; wdataset2.dY=ones(size(wdataset2.Y))*1e6; wdataset2.meantime=ones(size(wdataset2.Y))*1e5;
                        wdataset2.time1=wdataset2.meantime; wdataset2.time2=wdataset2.meantime;
                    else
                        wdataset2=wdataset; trainsub=trainsub0; Nincl(qqq+1)=length(trainsub);                   
                    end
                    sensnames{qqq+1} = 'All';
                    senscoords(qqq+1,:)=[NaN NaN];
                else
                    wdataset2=wdataset;
                    
                    if doexcl
                        trainsub=intersect(trainsub0,find(wdataset2.datid~=dosites(qqq)));
                    else                   
                        trainsub=intersect(trainsub0,find(wdataset2.datid==dosites(qqq))); Nincl(qqq+1)=length(trainsub);
                    end
                    sensnames{qqq+1}=dositenames{qqq};
                    senscoords(qqq+1,:)=docoords(qqq,:);
                end
                disp(sensnames{qqq+1});
                
                if length(trainsub)>0
                    
                     [wf,wsd,wV]=RegressHoloceneDataSets(wdataset2,wtestsitedef,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks(1,:),wtestt,refyear,collinear);
                    wf2=Mdiff*wf; wsd2=sqrt(diag(Mdiff*wV*Mdiff'));
                else
                    wf2=NaN; wsd2=NaN;
                end
                
                
                sensf(qqq+1,doexcl+1) = wf2;
                senssd(qqq+1,doexcl+1) = wsd2;
                
            end
            
        end
        
        varnormalizer = senssd(1,2)^2;
        sensvarnormed = senssd.^2/varnormalizer;
        
        fid=fopen(['sitesensitivity' labl '.tsv'],'w');
        fprintf(fid,'Site sensitivity tests\n');
        fprintf(fid,'GSL rate, %0.0f to %0.0f\n\n',wtestt);
        fprintf(fid,'Site\tLat\tLong\tN\tInclusive f (mm/y)\t1s\tpercent var reduction\tpercent var reduction/N\tExclusive f (mm/y)\t1s\tpercent var increase\tpercent var increase/N\n');
        for qqq=1:size(sensf,1)
            fprintf(fid,sensnames{qqq});
            fprintf(fid,'\t%0.3f',senscoords(qqq,:));
            fprintf(fid,'\t%0.0f',Nincl(qqq));
            fprintf(fid,'\t%0.3f',[sensf(qqq,1) 1*senssd(qqq,1) 100*(1-sensvarnormed(qqq,1)) 100*(1-sensvarnormed(qqq,1))/Nincl(qqq) sensf(qqq,2) 1*senssd(qqq,2) 100*(sensvarnormed(qqq,2)-sensvarnormed(1,1)) 100*(sensvarnormed(qqq,2)-sensvarnormed(1,1))/Nincl(qqq)]);
            fprintf(fid,'\n');
        end        
        fclose(fid);
        
         
        %% do author sensitivity tests
        
        wtestsitedef.sites=[0 1e6 1e6];
        wtestsitedef.names={'GSL'};
        wtestsitedef.names2={'GSL'};
        wtestsitedef.firstage=min(oldest);

        wtestt=[1000 1800];
        Mdiff = [-1 1]/diff(wtestt);
        
        
        authors={'AUTHORTEAM','Gehrels','Gehrels-FirstAuthor','Horton','Engelhart','Long','Tornqvist','Kemp','Sivan'};
        studies={
                 {'EH12','Massachusetts','Spain','France','North Carolina','New Jersey','Florida','Connecticut'},
                 {'Isle of Wight','Maine','New Zealand','Tasmania','Nova Scotia','Scotland','Iceland'},
                 {'Maine','New Zealand','Tasmania','Nova Scotia','Iceland'},
                 {'EH12','Massachusetts-Wood Island','Spain','France','North Carolina','New Jersey','Florida'},
                 {'EH12'},
                 {'Scotland','Isle of Wight','Greenland'},
            {'Louisiana'},
            {'North Carolina','Florida','New Jersey','Massachusetts-Wood Island'},
            {'Israel'},
                };
        
        [s,si]=sort(authors);
        authors=authors(si); studies=studies(si);
        
        clear sitesubs dositess docoords;
        for qqq=1:length(authors)
            cursub=[];
            for rrr=1:length(studies{qqq})
                matches=strfind(wdataset.sitenames,studies{qqq}{rrr});
                for sss=1:length(matches)
                    if length(matches{sss})==1
                        cursub=union(cursub,sss);
                    end
                end
                
            end
            sitesubs{qqq}=cursub;
        end

        authors{end+1}='TIDEGAUGES';
        
        sitesubs{end+1}=find(ismember(wdataset.siteid,unique(wdataset.datid(find(wdataset.istg)))));
        
        for qqq=1:length(authors)
            dositess{qqq}=wdataset.siteid(sitesubs{qqq});
            docoords{qqq}=wdataset.sitecoords(sitesubs{qqq},:);
        end
        
        trainsub0=find((wdataset.limiting==0));
        
        clear sensf senssd sensnames senscoords;

        for qqq=0:length(dositess)
            for doexcl=0:1
                wdataset2=wdataset;
                if qqq==0
                    if doexcl
                        trainsub=1;
                        wdataset2=wdataset; wdataset2.dY=ones(size(wdataset2.Y))*1e6; wdataset2.meantime=ones(size(wdataset2.Y))*1e5;
                        wdataset2.time1=wdataset2.meantime; wdataset2.time2=wdataset2.meantime;
                    else
                        trainsub=trainsub0; Nincl(qqq+1)=length(trainsub);                   
                    end
                    sensnames{qqq+1} = 'ALL';
                else
                    
                    if doexcl
                        trainsub=intersect(trainsub0,find(~ismember(wdataset2.datid,dositess{qqq})));
                    else                   
                        trainsub=intersect(trainsub0,find(ismember(wdataset2.datid,dositess{qqq}))); Nincl(qqq+1)=length(trainsub);
                        trainsub=union(trainsub,find(wdataset2.istg));
                    end
                    sensnames{qqq+1}=authors{qqq};
                end
                disp(sensnames{qqq+1});
                
                % add in tide gauges
                
                if length(trainsub)>0
                    
                     [wf,wsd,wV]=RegressHoloceneDataSets(wdataset2,wtestsitedef,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks(1,:),wtestt,refyear,collinear);
                    wf2=Mdiff*wf; wsd2=sqrt(diag(Mdiff*wV*Mdiff'));
                else
                    wf2=NaN; wsd2=NaN;
                end
                
                
                sensf(qqq+1,doexcl+1) = wf2;
                senssd(qqq+1,doexcl+1) = wsd2;
                
            end
            
        end
        
        varnormalizer = senssd(1,2)^2;
        sensvarnormed = senssd.^2/varnormalizer;
        
        fid=fopen(['authorsensitivity' labl '.tsv'],'w');
        fprintf(fid,'Author sensitivity tests\n');
        fprintf(fid,'GSL rate, %0.0f to %0.0f\n\n',wtestt);
        fprintf(fid,'Author\tN\tInclusive f (mm/y)\t1s\tpercent var reduction\tpercent var reduction/N\tExclusive f (mm/y)\t1s\tpercent var increase\tpercent var increase/N\n');
        for qqq=1:size(sensf,1)
            fprintf(fid,sensnames{qqq});
            fprintf(fid,'\t%0.0f',Nincl(qqq));
            fprintf(fid,'\t%0.3f',[sensf(qqq,1) 1*senssd(qqq,1) 100*(1-sensvarnormed(qqq,1)) 100*(1-sensvarnormed(qqq,1))/Nincl(qqq) sensf(qqq,2) 1*senssd(qqq,2) 100*(sensvarnormed(qqq,2)-sensvarnormed(1,1)) 100*(sensvarnormed(qqq,2)-sensvarnormed(1,1))/Nincl(qqq)]);
            fprintf(fid,'\n');
        end        
        fclose(fid);
        