function Nmat=nmatgen
%NMATGEN	Generate matrix of candidate ambiguity sets for
%               the case of 5 satellites (4 double-differences)
%               with uncertainty less than +/- 1 meter 
%               (+/- 5 wavelengths)
%
%	nmatgen
%
%    Output
%        Nmat = candidate ambiguity set matrix
%        Rows in Nmat range from [-5 -5 -5 -5] up to [5 5 5 5]
%
%    Note: Since this routine takes awhile to run, this has been
%          done for you and the result may be loaded directly
%          from nmat.mat

%    Copyright (c) 1996 by GPSoft
%
      kk = 0;
      for N1 = -5:5,
        for N2 = -5:5,
          for N3 = -5:5,
            for N4 = -5:5,
                kk = kk + 1;
                Nmat(kk,:) = [N1 N2 N3 N4];
            end
          end
        end
      end
