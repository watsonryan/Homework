function storerinnav(rinexfilename,matfilename)
%STORERINNAV   Program to load a RINEX2 navigation (ephemeris)
%              file and store the data in a MATLAB file
%              format (MAT-file)
%
%   storerinnav(rinexfilename,matfilename)
%
%   INPUTS
%  rinexfilename = Name of the ASCII text file containing the
%             RINEX2-formatted Navigation message data 
%  matfilename = filename for the MAT-file which will be created
%                If matfilename='day101nav' then a MAT-file will
%                be created called:  day101nav.mat
%
%  NOTE: make sure to put the names in single 
%  quotation marks (e.g.,  storerinnav('stkr2581.02n','day258nav')

%
%   Copyright (c) 2002-2003    Michael S. Braasch / GPSoft LLC
%
global ALPHA BETA UTC_A0 UTC_A1 UTC_TOT UTC_WN LEAP_SEC
global SV_ID_VEC TOC_YEAR TOC_MONTH TOC_DAY TOC_HOUR
global TOC_MINUTE TOC_SEC AF0 AF1 AF2
global IODE CRS DELTAN MZERO CUC ECCEN CUS SQRTSMA
global TOE CIC OMEGAZERO CIS IZERO CRC ARGPERI OMEGADOT
global IDOT CODES_ON_L2 TOE_WN L2_P_FLAG
global URA SV_HEALTH TGD IODC TRANS_TIME_OF_MESSAGE
global FIT_INTERVAL SPARE1 SPARE2
%
loadrinexn(rinexfilename)
%
save(matfilename)
