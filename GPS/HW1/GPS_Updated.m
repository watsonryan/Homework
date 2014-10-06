function [Clock_Bias_Estimate,XYZ_Estimate,pdop,tdop,gdop] = GPS_Updated(nomXYZ,nomCLOCK,pr_OBS,sats_XYZ,numSat,sigma_URE,W)
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
    prComputed(i)=norm(sats_XYZ(i,:)-nomXYZ)+nomCLOCK*c;
    
    % calculate the unit vectors from nominal position
    % to satellite location 
    % which make up the partials 
    % of the observation model
    uNom2Sat(i,:)=(sats_XYZ(i,:)-nomXYZ)/norm(sats_XYZ(i,:)-nomXYZ);
end
% Add the weighted matrix
% W = diag(numSat);
% add the partial of the clock bias, observation model
G=horzcat(uNom2Sat,ones(numSat,1));

% form the innovation vector
deltaRho=prComputed-pr_OBS;

% find the delta solution
dX=inv(G'*W*G)*G'*W*deltaRho;

H = inv(G'*G);

% calculate and solve for covariance matrix
P_error = sigma_URE^2 * eye(numSat);
cov_dX = sigma_URE^2 * H;
sigmaX = cov_dX(1,1);
sigmaY = cov_dX(2,2);
sigmaZ = cov_dX(3,3);
sigmat = cov_dX(4,4);

% Calculate different dilutions of precision
pdop = sqrt(sigmaX^2+sigmaY^2+sigmaZ^2);
tdop = sqrt(sigmat^2);
gdop = sqrt(sigmaX^2+sigmaY^2+sigmaZ^2+sigmat^2);
% the clock bias has a negative sign that must be handled
XYZ_Estimate=nomXYZ+dX(1:3)';
Clock_Bias_Estimate=nomCLOCK+-1*dX(4);
end
