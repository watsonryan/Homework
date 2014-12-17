%  sndemo27.m
%
%  DGPS Example 6D  -  Using 7 sv's
%
%  Double-Difference Carrier DGPS - Full Ambiguity Resolution
%
%  To simplify the program, we are taking advantage of the fact that
%  we know certain satellites are visible for the entire run
%
clear all;
close all;
%
%
load nmat               % nmat is a matrix of candidate ambiguities
                        % for a one meter (+/- 5 wavelength) uncertainty
                        % for four double-differences (5 satellites)
%
threshold = 0.1;        % threshold to reject/retain ambiguity sets
%
lambda=299792458/1575.42e6;
samat=sagen(24,3600,3,54321);
mpmatusr=mpgen(24,3600,1,22222);
mpmatref=mpgen(24,3600,1,33333);
refllh = [0*pi/180 0.1*pi/180 0];
refxyz = llh2xyz(refllh);
usrllh = [0.01*pi/180 0.04*pi/180 0];
usrxyz = llh2xyz(usrllh);

loadgps
initref=[0 0 0 0];
initusr=[0 0 0 0];
svflg1=0; svflg2=0; ambflg1=0; ambflg2=0;        % ambiguity resol. flags
prrefmat=[]; adrrefmat=[]; prusrmat=[]; adrusrmat=[];
i=0;  ii=0; parmag=[]; jpar=0; persistance=0; N_previous=[];

bar1 = waitbar(0,'Initializing (200 s) ...  ');

randn('state',74347098);
for t = 1:1:3600,

    i=i+1;

%%%Reference Station  -  Raw measurements and carrier-smoothing
    clear svmatref svmatusr svidref svidusr prref prusr
    [svmatref,svidref] = gensv(refxyz,t,15);
    [trprref,tradrref] = genrng(1,refxyz,svmatref,svidref,t,0); % Error free
    [prref,adrref] = genrng(1,refxyz,svmatref,svidref,t,...
                            [1 1 1 0.1 1],samat,mpmatref);
    [prsmref,prrefmat,adrrefmat]=...
                   hatch(prref,adrref,svidref,100,prrefmat,adrrefmat);

%%%User  -  Raw measurements and carrier-smoothing
    [svmatusr,svidusr] = gensv(usrxyz,t,15);
    [trprusr,tradrusr] = genrng(2,usrxyz,svmatusr,svidusr,t,0);  % Error free
    [prusr,adrusr] = genrng(2,usrxyz,svmatusr,svidusr,t,...
                            [1 1 1 0.1 1],samat,mpmatusr);
    [prsmusr,prusrmat,adrusrmat]=...
                   hatch(prusr,adrusr,svidusr,100,prusrmat,adrusrmat);
%
    if i > 199,        % Wait for 200 seconds to allow carrier-smoothing
                       % to produce good initial estimate of ambiguities

		if i == 200, 
			close(bar1);	  % end initialization bar
      end
      
		fprintf(1,'time is = %i \n',t)
      ii = ii + 1;
%                      % Estimate center of baseline between ref and user
      if ii == 1,
         estusr = olspos(prusr,svmatusr,initusr);  initusr=estusr;
         baspos = mean([refxyz(1:3); estusr(1:3)]);
      else
         baspos = mean([refxyz; ddestusr]);
      end
%
      j=0; maxel=0;
%%%      For this example, arbitrarily choose sv's 3,6,7,9,16,18,19
      for k = 1:max(size(svidref)),
           id=svidref(k);
                                                 % arb loop start
           if ( (id==3)|(id==6)|(id==7)|(id==9)|(id==16)|(id==18)|(id==19) ),
           if find(svidusr==id),
              j = j + 1;
                     % This section is storing the ranges (pseudorange and
                     % integrated Doppler; measured and true) which are
                     % common both to the reference and the user
              compr(j,1) = prref(k); 
              compr(j,2) = prusr(find(svidusr==id));
              compr(j,3) = trprref(k);
              compr(j,4) = trprusr(find(svidusr==id));
              compr(j,5) = prsmref(k);
              compr(j,6) = prsmusr(find(svidusr==id));
              comadr(j,1) = adrref(k);
              comadr(j,2) = adrusr(find(svidusr==id));
              comadr(j,3) = tradrref(k);
              comadr(j,4) = tradrusr(find(svidusr==id));
              comsv(j,:) = svmatref(k,:);
              vec = comsv(j,:) - baspos;        % Compute unit vector from
              unit_vec(j,:) = vec / norm(vec);  % baseline to satellite

              if ii == 1,                   % choose key sv at first epoch
                 enu = xyz2enu(comsv(j,:),baspos);
                 el(j) = atan2(enu(3),norm(enu(1:2)));
                 if el(j) > maxel,
                    maxel = el(j);
                    keyid = j;
                 end
              end
           end
           end                                              % arb loop end
           if ( (svflg1==0)&(j==5) ), break; end  % Here we are limiting
           if ( (svflg2==0)&(j==6) ), break; end  % the number of satellites
                                                  % used depending upon
                                                  % how many candidate
                                                  % ambiguity sets remain
      end  
      numcomsv = j;
      fprintf('Number of Satellites Used Now = %i \n',numcomsv)

%
      sd_pr = compr(:,2) - compr(:,1);         % Compute single-differences
      trsd_pr = compr(:,4) - compr(:,3);
      smsd_pr = compr(:,6) - compr(:,5);
      sd_adr = comadr(:,2) - comadr(:,1);
      trsd_adr = comadr(:,4) - comadr(:,3);
      j = 0;
      for k = 1:numcomsv,
          if k ~= keyid,                      % Compute double-differences
                                              % and data matrix
             j = j + 1;
             dd_pr(j,1) = sd_pr(keyid) - sd_pr(k);
             trdd_pr(j,1) = trsd_pr(keyid) - trsd_pr(k);
             smdd_pr(j,1) = smsd_pr(keyid) - smsd_pr(k);
             dd_adr(j,1) = sd_adr(keyid) - sd_adr(k);
             trdd_adr(j,1) = trsd_adr(keyid) - trsd_adr(k);
             H(j,:) = unit_vec(keyid,:) - unit_vec(k,:);
          end
      end
%
%
      if ii == 1,      % The 200 seconds of carrier-smoothing is done
                       % and now we start with just 5 satellites 
         tr_N_lambda = trdd_adr - trdd_pr;   % True ambiguities (meters)
         true_NL = tr_N_lambda';
         true_N = (1/lambda)*true_NL        % True amb. in wavelengths
         smest_N = round( (1/lambda)*(dd_adr - smdd_pr)' )
         smest_NL = lambda*smest_N;   % Carrier-smoothed code estimate
                                      % of ambiguities
%%%         rawest_NL = (dd_adr - dd_pr)';
      end
%
      if ambflg1 == 1,  % The number of candidate ambiguities for the
                        % 5 satellite case has just been reduced below 1000
                        % and now we add the 6th satellite
         tr_N_lambda = [tr_N_lambda; trdd_adr(5)-trdd_pr(5)];
         true_NL = tr_N_lambda';
         true_N = (1/lambda)*true_NL
         smest_N = [smest_N round((1/lambda)*(dd_adr(5)-smdd_pr(5)))]
         smest_NL = lambda*smest_N;
         ambflg1=0;
      end
%
      if ambflg2 == 1,  % The number of candidate ambiguities for the
                        % 6 satellite case has just been reduced below 1000
                        % and now we add the 7th satellite
         tr_N_lambda = [tr_N_lambda; trdd_adr(6)-trdd_pr(6)];
         true_NL = tr_N_lambda';
         true_N = (1/lambda)*true_NL
         smest_N = [smest_N round((1/lambda)*(dd_adr(6)-smdd_pr(6)))]
         smest_NL = lambda*smest_N;
         ambflg2=0;
      end
%
      [Q,R] = qr(H);   % Compute QR decomposition to aid in determining
                       % the parity vector for each candidate 
                       % ambiguity set
      Qt = Q';
      Q_beta = Qt(1:3,:);
      U = R(1:3,:);
%
%%      base_est = inv(U)*Q_beta*(dd_adr - tr_N_lambda);   % Steady-state
%%      ddestusr = refxyz - base_est';                     % carrier-phase
%%      err_adr(ii,:) = ddestusr - usrxyz;                 % solution
%
%%      base_est = inv(U)*Q_beta*dd_pr;          % Double-difference code
%%      ddestusr = refxyz - base_est';           % position solution
%%      err_pr(ii,:) = ddestusr - usrxyz;
%
      base_est = inv(U)*Q_beta*smdd_pr;         % Double-difference
      ddestusr = refxyz - base_est';            % carrier-smoothed code
      err_smpr(ii,:) = ddestusr - usrxyz;       % position solution
%
      [r,c] = size(Qt);
      Q_parity = Qt(4:r,:);          % The portion of Q transpose used
                                     % to compute the parity vector
%
      par_mag0(ii,1)=norm(Q_parity*(dd_adr-tr_N_lambda));  % Magnitude of
                                                           % the parity
                                                           % vector for the
                                                           % true ambiguity
                                                           % set
%
      Ncnttmp(ii) = max(size(Nmat));    % Number of candidate 
	                                     % ambiguity sets
      ipar = 0; Ntmp=[]; 
      for kk = 1:max(size(Nmat)),     % Loop through all remaining
                                      % candidate ambiguity sets
         tmp=norm(Q_parity*(dd_adr-(smest_NL'+lambda*Nmat(kk,:)')));
         if Ncnttmp(ii) > 199,
            if tmp < threshold,             % Retain this ambiguity set only
               ipar = ipar + 1;             % if the magnitude of its
               Ntmp(ipar,:) = Nmat(kk,:);   % parity vector is less than
            end                             % the threshold
         else
            if kk == 1,                 % We now have less than 200 
              jpar = jpar + 1;          % candidate ambiguity sets.
            end                         % We will store the magnitude of
            ipar = ipar + 1;            % the parity vector for each
            parmag(jpar,ipar) = tmp;    % ambiguity set for each epoch
            Ntmp(ipar,:) = Nmat(kk,:);
         end
      end
      Nmat = Ntmp;
      if jpar > 25,               % We've stored up enough values that
                                  % we can compute some statistics
         parstd = std(parmag);
         parmean = mean(parmag);
         [valmean,ordmean] = sort(parmean);   % Sort the means and standard
         [valstd,ordstd] = sort(parstd);      % deviations for the parity
                                              % magnitudes in ascending
                                              % order
%      
            N_current = Nmat(ordmean(1),:)       % It would be nice if the
                                                 % true ambiguity set was
            if persistance > 0,                  % the one whose
               if N_current == N_previous,       % parity vector magnitude
                  persistance = persistance + 1  % over time had the lowest
               else                              % mean.  This is true only
                  persistance = 1                % if you are willing to
                  N_previous = N_current;        % wait long enough.  Here
               end                               % we have arbitrarily 
               if persistance > 19,              % decided that 20 epochs,
                  time_est = t                   % IN A ROW, are sufficient.
                  tru_N_diff = true_N - smest_N  % Once this occurs, we
                  est_N_diff = N_current         % declare victory and quit.
                  return
               end
            else,
               persistance = 1
               N_previous = N_current;
            end,
      end                        %  end IF JPAR>25 loop
%
		Ncnttmp(ii) = max(size(Nmat));
      fprintf('Remaining Number of Ambiguity Sets = %i \n\n',Ncnttmp(ii))

%
      if ( (Ncnttmp(ii)<1000)&(svflg1==1)&(svflg2==0) ),   % See the 
         svflg2=1;                                         % comments for
         ambflg2=1;                                        % the next IF
      end                                                  % loop.  Now we
                                                           % go from 6
                                                           % satellites to 7
%
      if ( (Ncnttmp(ii)<1000)&(svflg1==0) ),   % If we have less than 1000
         svflg1=1;                             % candidate ambiguity sets
         ambflg1=1;                            % and we have only been using
      end                                      % 5 satellites so far, set 
                                               % svflg1=1 so that we will
                                               % start using 6 satellites
                                               % and set ambflg1=1 so that 
                                               % we expand our candidate
                                               % ambiguity sets to include
                                               % the additional double-
                                               % difference
%
      if ( (ambflg1==1) | (ambflg2==1) ),  % The flag has been set so we
         Ntmp=[];                          % can expand the ambiguity sets
         for kk = 1:max(size(Nmat)),       % to include one additional
            Nz=Nmat(kk,:);                 % double-difference
         Ntemp=[Nz -5;Nz -4;Nz -3;Nz -2;Nz -1;Nz 0;Nz 1;Nz 2;Nz 3;Nz 4;Nz 5];
            Ntmp = [Ntmp; Ntemp];
         end
         Nmat = Ntmp;
      end
%
      clear sd_pr dd_pr H dd_adr trsd_pr sd_adr trsd_adr
      clear trdd_pr trdd_adr Q R Qt Q_beta Q_parity base_est
      clear svidref svidusr svmatref svmatusr t tradrref tradrusr
      clear trprref trprusr unit_vec vec U adrref adrusr c comadr compr 
      clear comsv el enu estref estusr id j k maxel numcomsv 
      clear prref prusr r smdd_pr smsd_pr
   else
      waitbar(i/200)
   end  % End i>199 loop
end  % End t loop
%
%