function [Clock_Bias_Estimate,XYZ_Estimate] = GPS_Trop_Iono(nomXYZ,nomCLOCK,pr_OBS,sats_XYZ,numSat,Trop,Iono)
% This function will take the nominal XYZ location, the nominal clock, the
% observed pseudorange, and the satellite information for each satellite
% and will output the estimate position and estimate clock bias.

% The inputs for this function are listed below with a description. 

% nomXYZ - rough guess of GPS location [m]
% nomCLOCK - rough guess of clock bias [m]** 
% pr_OBS - observes pseudoranges set forth by the receiver
% satsXYZ - satellite coordinate information from broadcast navigation message
% nSat - number of satellites used
format long g
c=299792458; % speed of light [m/s]

% pre allocate sizes
prComputed=zeros(numSat,1);
uNom2Sat=zeros(numSat,3);

for i=1:numSat;
    % calculate computed or expected psuedorange
    prComputed(i)=norm(sats_XYZ(i,:)-nomXYZ)+nomCLOCK*c+Trop(i)+Iono(i);
    
    % calculate the unit vectors from nominal position
    % to satellite location 
    % which make up the partials 
    % of the observation model
    uNom2Sat(i,:)=(sats_XYZ(i,:)-nomXYZ)/norm(sats_XYZ(i,:)-nomXYZ);
end

% add the partial of the clock bias, 
% observation model
G=horzcat(uNom2Sat,ones(numSat,1));
% form the innovation vector
deltaRho=prComputed-pr_OBS;
% find the delta solution
dX=inv(G'*G)*G'*deltaRho;

% the clock bias in has a negative sign that must be handled
XYZ_Estimate=nomXYZ+dX(1:3)';
Clock_Bias_Estimate=nomCLOCK+-1*dX(4);
end
