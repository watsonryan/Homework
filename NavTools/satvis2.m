function vismat = satvis2(usrxyz,tstart,tstop,tinc,daynum,filename,plotflg,maskang)
%SATVIS2		Generate satellite visibility plots with Rinex ephemeris data
%
%	vismat = SATVIS2(usrxyz,tstart,tstop,tinc,daynum,filename,plotflg,maskang)
%
%   INPUTS
%	usrxyz = user position in cartesian ECEF coordinates
%       tstart,tstop,tinc = start, stop and increment times for satellite
%                           visibility assessment (GPS time of day in seconds)
%       daynum = day number of the week (0=Sunday, 1=Monday, ..., 6=Saturday)
%       filename = name of the ephemeris file; note the ephemeris must be in 
%                  RINEX2 format; the filename must be in quotes 
%                  (e.g., 'stkr2581.o2n')
%       plotflg = 1 for total number of visible satellites versus time
%                 2 for bar plot of visibility of each satellite versus time
%       maskang = Visibility mask angle.  Optional argument.  Default is 5 degrees.
%
%    OUTPUT
%       vismat = matrix of satellite visibility parameters.
%                If plotflg = 1, vismat(:,1) is the total number of visible
%                satellites for the time points specified in vismat(:,2)
%                If plotflg = 2, the last column of VISMAT is the vector of time
%                points at which visibility was evaluated.  Columns 1 through 32
%                (or greater if geostationary satellites are visible) correspond
%                satellite prn's and contain a zero if that satellite is not
%                visible at the given time point.
%
%     See also:   SATVIS3, SATVIS

%	M. & S. Braasch 12-2002; Revised 04-2003
%	Copyright (c) 2002-2003 by GPSoft
%	All Rights Reserved.
%

if nargin<8, maskang=5; end
if nargin<7, error('insufficient number of input arguments'),end

global SQRTSMA
loadrinexn(filename);

id = find(SQRTSMA);

i=0;
h = waitbar(0,'Satvis2: Determining visible satellites versus time');
for t = tstart:tinc:tstop,

   clear svid
   gpstow = daynum*86400 + t;   % GPS time of week in seconds
   k = 0;
   for j = 1:length(id),
      svxyz = svposeph(id(j),gpstow);
      svenu = xyz2enu(svxyz,usrxyz);
      el = (180/pi)*atan(svenu(3)/norm(svenu(1:2)));
      if el >= maskang,
          k = k + 1;
          svid(k) = id(j);
      end
   end
    i = i + 1;
    if plotflg == 1,
       vismat(i,1) = t/3600;
       vismat(i,2) = length(svid);
    elseif plotflg == 2,
       vismat(i,:) = zeros(1,40);  
       vismat(i,40) = t;
       for isv = 1:length(svid),
           vismat(i,svid(isv)) = svid(isv);
       end
    end
    waitbar((t-tstart)/(tstop-tstart),h)
end
close(h)
numpts=i;

close
if plotflg == 1,
	bar(vismat(:,1),vismat(:,2))
        axis([tstart/3600 tstop/3600 0 max(vismat(:,2))+1])
           title('Total Number of Visible SVs - RINEX Ephemeris')
        xlabel('GPS Time of Day (hours)')
        ylabel('number of satellites')
elseif plotflg == 2,
    [nrow,ncol] = size(vismat);
	time=vismat(:,ncol);
    highsv = max(max(vismat(:,1:ncol-1)));
        plot((time(1)/3600)*ones(1,highsv),vismat(1,1:highsv),'*')
        axis([tstart/3600 tstop/3600 1 highsv+1])
           title('Satellite Visibility - RINEX Ephemeris')
        xlabel('GPS Time of Day (hours)')
        ylabel('satellite id number')        
        hold on
        for j = 2:numpts
            plot((time(j)/3600)*ones(1,highsv),vismat(j,1:highsv),'*')
        end
end
hold off
