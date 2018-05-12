% Demo how to use unspec in Matlab or Octave:
% extract a map.

% SETUP:
fname = 'mcascans.spec'	% experimental file
scan_nb = 5		% which scan to extract
% detector_name = 'det' % name of the detector(s)
detector_name = 'corr'  % name of the detector(s)
%detector_name = 'det,corr'  % name of the detector(s)
% END OF SETUP

% note: "=" is the scanned motor name (autodetected)
cmd = sprintf('unspec %s tmp.dat -2 -# -w =,%s -s %i -Py', fname, detector_name, scan_nb)
system(cmd);

a=load('tmp.dat');
motor = a(:,1);
det   = a(:,2);
map   = a(:,3:end);

figure(1);
plot(motor,det);
title('demo using unspec: draw scan in a graf')
xlabel('motor (deg)');

figure(2);
imagesc(motor, 1:size(map,2), map);
axis xy; colorbar
xlabel('motor (deg)'); ylabel('pixel');
title('demo using unspec: draw MCA spectra as a map')

% end
