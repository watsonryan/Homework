function [idx] = findHighestElevationSat(satsXYZ, nomXYZ)
%[idx] = findHighestElevationSat(satsXYZ)

for i=1:length(satsXYZ)
    for j=1:32
        tempENU=xyz2enu([satsXYZ(i,j,1),satsXYZ(i,j,2),satsXYZ(i,j,3)], nomXYZ);
        el_R(j)=atan2(tempENU(3),sqrt(tempENU(1)^2+tempENU(2)^2));
    end
    idx(i)=find(el_R==max(el_R));
end
end

