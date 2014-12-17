function vismat = satvis3(usrxyz,tstart,tstop,tinc,daynum,filename,plotflg,currweek,maskang)
%SATVIS3		Generate satellite visibility plots with Yuma almanac data
%
%	vismat = SATVIS3(usrxyz,tstart,tstop,tinc,daynum,filename,plotflg,currweek,maskang)
%
%   INPUTS
%	usrxyz = user position in cartesian ECEF coordinates
%       tstart,tstop,tinc = start, stop and increment times for satellite
%                           visibility assessment (GPS time of day in seconds)
%       daynum = day number of the week (0=Sunday, 1=Monday, ..., 6=Saturday)
%       filename = name of the almanac file; note the almanac must be in 
%                  Yuma format; the filename must be in quotes 
%                  (e.g., 'yuma19.txt')
%       plotflg = 1 for total number of visible satellites versus time
%                 2 for bar plot of visibility of each satellite versus time
%       currweek = GPS week number for the day (daynum) and time (tstart, etc) specified
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
%     See also:   SATVIS2, SATVIS

%	M. & S. Braasch 12-2002; Revised: 07-2003; 06-2004
%	Copyright (c) 2002-2004 by GPSoft
%	All Rights Reserved.
%

if nargin<9, maskang=5; end
if nargin<8, error('insufficient number of input arguments'),end

global SVIDV MV OMGV RV INCLV TOEV HEALTHV ECCENV
global OMGDOTV ARGPERIV AF0V AF1V WEEKV
loadyuma(filename);

id = find(SVIDV);

i=0;
h = waitbar(0,'Satvis3: Determining visible satellites versus time');
for t = tstart:tinc:tstop,
   clear svid
   gpstow = daynum*86400 + t;   % GPS time of week in seconds
   k = 0;
   for j = 1:length(id),
       N = id(j);
      svxyz = svposalm(RV(N),TOEV(N),MV(N),OMGV(N),INCLV(N),...
           gpstow,ECCENV(N),ARGPERIV(N),OMGDOTV(N),WEEKV(N),currweek);
      svenu = xyz2enu(svxyz,usrxyz);
      el = (180/pi)*atan(svenu(3)/norm(svenu(1:2)));
      if el >= maskang,
          k = k + 1;
          svid(k) = SVIDV(id(j));
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
           title('Total Number of Visible SVs - Yuma Almanac')
        xlabel('Time of Day (hours)')
        ylabel('number of satellites')
elseif plotflg == 2,
    [nrow,ncol] = size(vismat);
	time=vismat(:,ncol);
    highsv = max(max(vismat(:,1:ncol-1)));
        plot((time(1)/3600)*ones(1,highsv),vismat(1,1:highsv),'*')
        axis([tstart/3600 tstop/3600 1 highsv+1])
           title('Satellite Visibility - Yuma Almanac')
        xlabel('Time of Day (hours)')
        ylabel('satellite id number')        
        hold on
        for j = 2:numpts
            plot((time(j)/3600)*ones(1,highsv),vismat(j,1:highsv),'*')
        end
end
hold off
