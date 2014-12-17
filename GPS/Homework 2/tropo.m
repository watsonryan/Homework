function [Tropo] = tropo(z,Tot_pressure, Par_pressure,T, XYZ_Estimate);


    llh(z,:) = xyz2llh(XYZ_Estimate(:,z)');

    Td = .002277*(1+.0026*(cos(llh(z,1)))+(.0028*(llh(z,3))))*Tot_pressure;
    Tw = .002277*((1255/T)+.05)*Par_pressure;
    Tropo(:,z) = Td+Tw;

    Tropo(:,z) = Td+Tw;