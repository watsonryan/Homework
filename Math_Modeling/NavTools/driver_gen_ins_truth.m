%
%   driver_gen_ins_truth

earth_shape_flag = 1;    % 1 for ellipsoidal
vert_mech = 0;   % 0 for north-pointing mechanization
lat_start = 39;   % latitude (in deg) of starting point of vehicle trajectory
lon_start = -81;  % longitude (in deg) of starting point of vehicle trajectory

gen_ins_truth('flt_prof_1_dat',earth_shape_flag,vert_mech,...
    lat_start,lon_start,'prof_1_ins_truth')
