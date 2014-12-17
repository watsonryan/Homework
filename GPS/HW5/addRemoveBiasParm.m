function [xo,Po]=addRemoveBiasParm(xi,Pi,nNB,ddIDin,ddIDout, resetSig)
% function that handles changing state vector due to satellites 
% either being picked up or dropped from tracking.  It adds new ambiguties
% to state vector by initializing the diagonal of the P matrix for that
% paramter.  It pulls out row/column of ambiguties that are no longer
% needed  due a satellite dropping.
%
% 
% note THIS FUNCTION ASSUMES 
% 1) DUAL FREQUENCY DATA (i.e. 2 ambiguities for
%    each double difference satellite pair)
% 2) The reference satellite does not changes for the previous epoch to the
%    current epoch


% Function inputs:
% xi -- state vector (non ambiguity and ambiguity estimates from previous epoch)
% Pi -- error-covariance vector (from previous epoch, corresponds to xi)
% nNB -- number of non-ambiguity parameters in the state vector 
%        (this is typically = 3, for the delta X, delta Y, delta Z)
% ddIDin -- list of PRN # correponding to the Double Differences from last epoch 
%           (that is, PRN of non-Reference satellite used)
% ddIDout -- list of PRN # correponding to the Double Differences from
%            current epoch (that is, PRN of non-Reference satellite used)
% resetSig -- If a new DD is added, this is the uncertainty used for the
%             new ambiguity term added to the state vector

% Author: Jason Gross



nNewDD=length(ddIDout);
nOldDD=length(ddIDin);

n=2*nNewDD+nNB;
no=nNewDD+nNB;
Po=zeros(n,n);
xo=zeros(n,1);


% if # DD grows, add parm
% and use resetSig to initialize undertainty
% Note this portion of the code could be updated
% to use DD pseudorange to provide an initial guess of
% ambiguity, currently it is just initialized to zero.
if nNewDD>nOldDD
    
    % start with an identity of reset sig
    Po=resetSig*eye(n);
   
   
    % just pass through the non-ambiguty portion of
    % the state and error-covariance matrix
    Po(1:nNB,1:nNB)=Pi(1:nNB,1:nNB); 
    xo(1:nNB)=xi(1:nNB);
    
    
    for id=1:length(ddIDout);
        % find where new PRN is and add that portion of state vector
        if(isempty(find(ddIDin==ddIDout(id))))
            
            xo(nNB+id)=0;
            xo(nNB+id+nNewDD)=0;
            Po(nNB+id,nNB+id)=resetSig;
            Po(nNB+id+nNewDD,nNB+id+nNewDD)=resetSig;
            
        else
         % otherwise make sure to shuffle old values correctly 
            kk=find(ddIDin==ddIDout(id));
            % for L1
            xo(nNB+id)=xi(nNB+kk);
            % for L2
            xo(nNB+id+nNewDD)=xi(nNB+kk+nOldDD);
   
            % for L1            
            Po(nNB+id,nNB+id)=Pi(nNB+kk,nNB+kk);
            % for L2
            Po(nNB+id+nNewDD,nNB+id+nNewDD)=Pi(nNB+kk+nOldDD,nNB+kk+nOldDD);
        end
    end
% in this case we need to remove a part of xi and pi     
elseif nNewDD<nOldDD
    
    % determine number of dropped satellites 
    % and save which ones have been dropped
    drop=0;
    for id=1:length(ddIDin);
        if(isempty(find(ddIDout==ddIDin(id))))
            drop=drop+1;
            notincluded(drop)=ddIDin(id);
        end
    end
    
    % now store indices of dropped sats
    rmIndex=zeros(1,2*drop);
    for i=1:drop
        rmIndex(i)=find(ddIDin==notincluded(drop))+nNB;
        rmIndex(i+drop)=find(ddIDin==notincluded(drop))+nNB+nOldDD;
    end
    
    
    % effeciently remove the rows/columns that correspond
    % to the ambiguity of the satellite dropped
    Po=Pi;
    Po(rmIndex,:)=[];
    Po(:,rmIndex)=[];
    xo=xi;
    xo(rmIndex)=[];
    Po=Po;
    xo=xo;
else
    % if nDD didnt change, don't do anything but pass through.
    Po=Pi;
    xo=xi;
end
