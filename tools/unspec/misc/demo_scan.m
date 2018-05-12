% Demo how to use unspec in Matlab or Octave:
% extract a scan from MCA data.

% SETUP:
fname = 'manip.spec'	% experimental file
scan_nb = 12		% which scan to extract
% detector_name = 'det' % name of the detector(s)
detector_name = 'corr'  % name of the detector(s)
%detector_name = 'det,corr'  % name of the detector(s)
% END OF SETUP

% note: "=" is the scanned motor name (autodetected)
cmd = sprintf('unspec %s tmp.dat -2 -# -w =,%s -s %i', fname, detector_name, scan_nb)
system(cmd);

a=load('tmp.dat');
plot(a(:,1), a(:,2));
title('demo using unspec: draw scan in a graf')

% end
