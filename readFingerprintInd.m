function [fp,lo,la] = readFingerprintInd(fpname,subdir)

% [fp,lo,la] = readFingerprintInd(fpname,[subdir])
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Dec 21 09:07:12 EST 2013

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
	[fp0,lo,la] = readjxms_ssht([files(i).name(1:end-3)],0,0,NaN,L3);
	fp(:,:,i)=fp0;
	system(['rm /tmp/' files(i).name(1:end-3)]);	
end
cd(pd);
