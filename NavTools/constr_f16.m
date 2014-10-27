function x = constr_f16(x,gamma,rollrate,pitchrate,turnrate,coord,stab)
%CONSTR_F16     Function used by COSTF16.M to apply constraints to the
%               flight trajectory
%
%     x = constr_f16(x,gamma,rollrate,pitchrate,turnrate,coord,stab)
%
%  INPUTS
%     x = aircraft state vector; see F16_6DOF.M for details
%     gamma = flight path angle (with respect to horizontal) in radians
%     rollrate = aircraft roll rate in radians/sec
%     pitchrate = aircraft pitch rate in radians/sec
%     turnrate = aircraft turn rate in radians/sec
%     coord = coordinated turn flag; 1=coordinated; 2=skidding turn
%      NOTE: Skidding turns are not currently supported by the software
%     stab = stability axis roll flag; 1=stability-axis roll;
%           2=body-axis roll
%
%  OUTPUT
%     x = aircraft state vector; see F16_6DOF.M for details
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%
% 
eps = 1e-10;
alpha = x(2);
beta = x(3);
phi = x(4);
if coord == 1,
    % coordinated turn logic here
    Gscript = turnrate*x(1)/32.2;
    a = 1 - Gscript*tan(alpha)*sin(beta);
    b = sin(gamma)/cos(beta);
    c = 1 + ( Gscript*cos(beta) )^2;
    num1 = Gscript*cos(beta);
    num2 = (a - b*b) + b*tan(alpha)*sqrt(c*(1-b*b)+(Gscript*sin(beta))^2);
    dum = num1*num2/( cos(alpha)*(a*a - b*b*(1+c*(tan(alpha)^2))) );
    phi = atan(dum);
    
    a = cos(alpha)*cos(beta);
    b = sin(phi)*sin(beta) + cos(phi)*sin(alpha)*cos(beta);
    num = a*b + sin(gamma)*sqrt(a*a - (sin(gamma)^2) + b*b);
    den = a*a - (sin(gamma)^2);
    theta = atan(num/den);
    
    x(4) = phi;
    x(5) = theta;

elseif abs(turnrate) > eps,
    % skidding turn logic here
    error('Skidding turns not currently supported; Set coordinated turn flag to 1')
else
    % non-turning flight logic here
    a = cos(alpha)*cos(beta);
    b = sin(phi)*sin(beta) + cos(phi)*sin(alpha)*cos(beta);
    x(4) = phi;
    d = x(2);
    
    if abs(phi) > eps,
        d = -x(2);     % inverted flight
    end
    
    if abs(gamma) > eps,
        sgocb = sin(gamma)/cos(beta);
        x(5) = d + atan( sgocb/sqrt(1-sgocb^2) );  % roc constraint
    else
        %x(5) = d;    % level
        den = a*a - (sin(gamma))^2;
        num = a*b + sin(gamma)*sqrt(den + b*b);
        if den == 0,
            x(5) = d;
        else
            x(5) = atan(num/den);
        end
    end

    x(7) = rollrate;
    x(8) = pitchrate;
    
    if stab == 1,               % stability axis roll
        x(9) = rollrate*sin(alpha)/cos(alpha);
    else
        x(9) = 0;            % body axis roll
    end
    
end    % end 'if coord' loop
