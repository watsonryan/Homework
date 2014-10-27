function storerinob(rinexfilename,matfilename)
%STORERINOB   Program to load a RINEX2 observation (measurement data)
%             file and store the data in a MATLAB file format (MAT-file)
%
%   storerinob(rinexfilename,matfilename)
%
%   INPUTS
%  rinexfilename = Name of the ASCII text file containing the
%             RINEX2-formatted observation (measurement) data 
%  matfilename = filename for the MAT-file which will be created
%                If matfilename='day258ob' then a MAT-file will
%                be created called:  day258ob.mat
%
%  NOTE: make sure to put the names in single 
%  quotation marks (e.g.,  storerinob('stkr2581.02o','day258ob')

%
%   Copyright (c) 2002-2003    Michael S. Braasch / GPSoft LLC
%
global SVID_MAT TOWSEC PHASE1 PHASE2 C1 P1 P2 D1 D2 S1 S2
global PHASE1LLI PHASE1SS PHASE2LLI PHASE2SS
global C1LLI C1SS P1LLI P1SS P2LLI P2SS
global MARKER_XYZ ANTDELTA OBSINT TIMESTART TIMESTOP CLOCKOFFSET
%
loadrinexob(rinexfilename)
%
save(matfilename)
