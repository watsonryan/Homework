function vismat = satvis(usrxyz,tstart,tstop,tinc,sysflg,plotflg,maskang)
%SATVIS		Generate satellite visibility plots with simulated constellations
%
%	vismat = SATVIS(usrxyz,tstart,tstop,tinc,sysflg,plotflg,maskang)
%
%   INPUTS
%	usrxyz = user position in cartesian ECEF coordinates
%       tstart,tstop,tinc = start, stop and increment times for satellite
%                           visibility assessment (GPS time of day in seconds)
%       sysflg = vector specifying which satellite systems are to be included:
%                1 = GPS
%                2 = Glonass
%                3 = GEO
%                4 = Galileo
%                For example, if you want both GPS and Galileo: sysflg = [1 4].
%       plotflg = 1 for total number of visible satellites versus time
%                 2 for bar plot of visibility of each satellite versus time
%       maskang = Visibility mask angle.  Optional argument.  Default is 5 degrees.
%
%    OUTPUT
%       vismat = matrix of satellite visibility parameters.
%                If plotflg = 1, vismat(:,1) is the total number of visible
%                satellites for the time points specified in vismat(:,2)
%                If plotflg = 2, the last column of vismat is the vector of time points
%                at which visibility was evaluated.  The previous columns
%                correspond to the satellite ID numbers in the constellation (either
%                contain a zero if that satellite is not visible at the given time point.
%
%    See also:  SATVIS2, SATVIS3

%	References: 
%                   Understanding GPS: Principles and Applications,
%	            Elliott D. Kaplan, Editor, Artech House Publishers,
%	            Boston, 1996.
%
%	M. & S. Braasch   Revised 07-2003
%	Copyright (c) 1996-2003 by GPSoft
%	All Rights Reserved.
%

global SVIDV

if nargin<7, maskang=5; end
if nargin<6, error('insufficient number of input arguments'),end

for i = 1:length(sysflg),
    if sysflg(i) == 1,
        loadgps
    elseif sysflg(i) == 2,
        loadglo
    elseif sysflg(i) == 3,
        loadgeo
    elseif sysflg(i) == 4,
        loadgalileo
    else
        error('SYSFLG elements must be in the range 1 - 4')
    end
end
totrng = max(SVIDV);

i=0;
for t = tstart:tinc:tstop,
    [svxyzmat,svid] = gensv(usrxyz,t,maskang);
    i = i + 1;
    if plotflg == 1,
       vismat(i,1) = t/3600;
       vismat(i,2) = max(size(svid));
    elseif plotflg == 2,
       vismat(i,:) = zeros(1,totrng+1);  vismat(i,totrng+1) = t;
       for isv = 1:max(size(svid)),
           vismat(i,svid(isv)) = svid(isv);
       end
    end
end
numpts=i;

close
if plotflg == 1,
	stairs(vismat(:,1),vismat(:,2))
    axis([tstart/3600 tstop/3600 0 max(vismat(:,2))+1])
    title('Satellite Visibility')
    xlabel('Time of Day (hours)')
    ylabel('number of satellites')
elseif plotflg == 2,
	time=vismat(:,totrng+1);
        plot((time(1)/3600)*ones(1,totrng),vismat(1,1:totrng),'*')
        axis([tstart/3600 tstop/3600 1 totrng+1])
	   title('Satellite Visibility')
        xlabel('Time of Day (hours)')
        ylabel('satellite id number')        
        hold on
        for j = 2:numpts
            plot((time(j)/3600)*ones(1,totrng),vismat(j,1:totrng),'*')
        end
end
hold off
