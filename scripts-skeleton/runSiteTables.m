% For each site, make a plot
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2017-07-25 14:17:13 -0400


subsite=find((wdataset.siteid>10000));
maxdistfrom=0.1;
maxerror=1000;
wtestlocs=testlocs{iii};
dosites=wtestlocs.sites(:)';
doNoiseMask=1;

for kkk=dosites

       
    disp(wtestlocs.names{kkk});
    
    datsub=find(wtestlocs.reg==wtestlocs.sites(kkk,1));
    
    for timesteps=[100]
        [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(wtestlocs.X(datsub,3),wtestlocs.reg(datsub),wtestlocs.sites(kkk,1),f2s{iii}(datsub,doNoiseMask),V2s{iii}(datsub,datsub,doNoiseMask),[],testt(1),testt(end),0,timesteps,wtestlocs.names(kkk));


        fid=fopen(['siteplot-' labl '-' wtestlocs.names{kkk} '-' num2str(timesteps) 'y' '.tsv'],'w');
        fprintf(fid,outtable);
        fclose(fid);
        
        
    end    
end
