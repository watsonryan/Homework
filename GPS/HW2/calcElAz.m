%This Function Compute Azimuth and Elevation of satellite from reciever 
%************************************************************************           
%    ==================================================================
%    Input :                                                            *
%        Pos_Rcv       : XYZ position of reciever               (Meter) *
%        Pos_SV        : XYZ matrix position of GPS satellites  (Meter) *
%    Output:                                                            *
%        E             :Elevation (Rad)                                 *
%        A             :Azimuth   (Rad)                                 *
%************************************************************************           


function [E,A]=calcElAz(Pos_Rcv,Pos_SV)
ENU=xyz2enu(Pos_SV(:),Pos_Rcv(:));
Elevation=atan2(ENU(3),norm(ENU));
Azimuth=atan2(ENU(1)/norm(ENU),ENU(2)/norm(ENU));
E=Elevation;
A=Azimuth;