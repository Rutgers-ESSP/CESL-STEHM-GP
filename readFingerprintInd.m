function [fp,lo,la] = readFingerprintInd(fpname,subdir)

% [fp,lo,la] = readFingerprintInd(fpname,[subdir])
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Nov 10 17:08:32 EST 2014

defval('fpname','gis');
defval('subdir',[pwd '/IFILES']);

L3=512;

pd=pwd;
cd(subdir);
files=dir(['*_' fpname '.mn.gz']);

for i=1:length(files)
	system(['cp ' files(i).name ' /tmp']);
end
cd('/tmp');
clear fp fpname;
for i=1:length(files)
	system(['gzip -d ' files(i).name ]);
        if exist('ssht_inverse')
            [fp0,lo,la] = readjxms_ssht([files(i).name(1:end-3)],0,0,NaN,L3);
        else
            [fp0,lo,la] = readjxms([files(i).name(1:end-3)],0,0,NaN,L3);
        end
        
	fp(:,:,i)=fp0;
	system(['rm /tmp/' files(i).name(1:end-3)]);	
end
cd(pd);
