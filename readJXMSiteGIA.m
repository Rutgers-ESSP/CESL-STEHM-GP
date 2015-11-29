function [testt,testsl,sites,icehist,solidearth] = readJXMSiteGIA(files,testt)

% [testt,testsl,sites,icehist,solidearth] = readJXMSiteGIA(files,testt)
%
% Reads table of GIA projections from Jerry Mitrovica.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun May 25 22:04:07 EDT 2014

defval('testt',[]);
defval('files','*.out.gz');
q=strfind(files,'/');
if length(q)>0
    filepath=files(1:q(end));
    files=files((q(end)+1):end);
else
    filepath='.';
end

files=dir(fullfile(filepath,files));

sites={}; clear age sl siteid;
for i=1:length(files)
    dogz = strfind(files(i).name,'.gz');
    disp(files(i).name);
    if dogz
        system(['gzip -d ' fullfile(filepath,files(i).name)]);
        files(i).name=files(i).name(1:end-3);
    end
    tokn = regexp(files(i).name,'rsl_(\w+)_(\w+).out','tokens'); tokn=tokn{1};
    icehist{i} = tokn{1}; solidearth{i} = tokn{2};

    fid = fopen(fullfile(filepath,files(i).name));
    S=textscan(fid,'%s'); S=S{1};
    fclose(fid);
    if dogz
        system(['gzip ' fullfile(filepath,files(i).name)]);
    end

    readingseries=0; curlabel=''; age{i}=[]; sl{i}=[]; siteid{i}=[]; readinglat=0;
    for jj=1:length(S)
        if readingseries == 1
            if length(str2num(S{jj})) == 0
                readingseries = 0;
                curlabel='';
            end
        end
        if length(str2num(S{jj})) == 0
            readingseries = 0;
            curlabel=[curlabel S{jj}];
        elseif readingseries == 0
            %disp(curlabel);
            indx=strmatch(curlabel,sites,'exact');
            if length(indx)==0
                sites={sites{:},curlabel};
                indx = length(sites);
                disp(curlabel);
            end
            readingseries = 1;
            readingage = 1;
        end
        if readingseries == 1
            if readingage
                age{i} = [age{i} str2num(S{jj})];
                readingage = 0;
                siteid{i} = [siteid{i} indx];
            else
                sl{i} = [sl{i} str2num(S{jj})];
                readingage = 1;
            end
        end
    end
    
end

if length(testt)==0
    testt=unique([age{:}]);
end

clear testsl;
for i=1:length(sl)
    for jj=1:length(sites)
        sub=find(siteid{i}==jj);
        testsl(:,jj,i) = interp1(age{i}(sub),sl{i}(sub),testt,'linear');
    end
    [m,mi]=min(testt);
    testsl(:,:,i) = bsxfun(@minus,testsl(:,:,i),testsl(mi,:,i));
end