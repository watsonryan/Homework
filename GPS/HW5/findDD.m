function [dd_range1, dd_range2, dd_phase1, dd_phase2] = findDD(refStation,usrStation,refSat)
%[ddrange1, ddrange2, ddphase1, ddpahse2] = findDD(range1, range2, phase1, phase2)

for i=1:length(refSat)
    %SD
    sd_range1(i,:)=usrStation.C1(:,i)-refStation.C1(:,i);
    sd_range2(i,:)=usrStation.P2(:,i)-refStation.P2(:,i);
    sd_phase1(i,:)=usrStation.PHASE1(:,i)-refStation.PHASE1(:,i);
    sd_phase2(i,:)=usrStation.PHASE2(:,i)-refStation.PHASE2(:,i);

    %DD
    for j=1:32
        if(j~=refSat(i))
            dd_range1(i,j)=sd_range1(j)-sd_range1(refSat(i));
            dd_range2(i,j)=sd_range2(j)-sd_range2(refSat(i));
            dd_phase1(i,j)=sd_phase1(j)-sd_phase1(refSat(i));
            dd_phase2(i,j)=sd_phase2(j)-sd_phase2(refSat(i));
        end
    end
end

end

