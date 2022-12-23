function [So So2 T2 T2s rmse] = SAGE_fast(imgs,TEs,te2,mask,displayON)
% ZZjing: dsage fitting (can be multi-groups imgs, every group contains one ge, one mix, one se)
% the size of AA depend on the number of groups


% MKM  - fit T2 and T2* using \ method for SAGE sequence
% uses 4 parameter fit from equation in Schmiedeskamp et al 2012
%
%
% Inputs
% imgs:          input images at each echo time
% TEs:            echo times
% te2:             time of 180
% mask:          mask of voxels to fit
% displayON;    enables display
%
% Outputs
% So So2 T2 T2s rmse : four parameter fit and rmse

sz=size(imgs);

for i = 1:size(TEs,2)
    
    ind = find(TEs(:,i)>te2(:,i),1);
    ind2 = find(TEs(:,i)>te2(:,i)*2,1);
    if (isempty(ind2))
        % jaejin 
        % from  : ind2 = size(TEs,1)+1;
        % to    : ind2 = size(TEs,1);
        ind2 = size(TEs,1);
    end
    if ind2 == ind
        ind2 = ind2+1;
    end
    indA=ind;ind2A=ind2;
    
    A1(:,:,i) = [ones(ind-1,1), zeros(ind-1,1), zeros(ind-1,1),-TEs(1:ind-1,i)];
    A2(:,:,i) = [zeros(ind2-ind,1), ones(ind2-ind,1), (2*te2(:,i)-2*TEs(ind:ind2-1,i)), (TEs(ind:ind2-1,i)-2*te2(:,i))];
    A3(:,:,i) = [zeros(size(TEs,1)-ind2+1,1), ones(size(TEs,1)-ind2+1,1), -2*te2(:,i).*ones(size(TEs,1)-ind2+1,1), (2*te2(:,i)-TEs(ind2:end,i))];
    
    A = [A1; A2; A3];
end

if size(A,3) == 1
    AA = cat(1,A(:,:,1));
elseif size(A,3) == 2
    AA = cat(1,A(:,:,1),A(:,:,2));
elseif size(A,3) == 3
    AA = cat(1,A(:,:,1),A(:,:,2),A(:,:,3)); % depend on the number of groups
elseif size(A,3) == 4
    AA = cat(1,A(:,:,1),A(:,:,2),A(:,:,3),A(:,:,4));
elseif size(A,3) == 5
    AA = cat(1,A(:,:,1),A(:,:,2),A(:,:,3),A(:,:,4),A(:,:,5));
elseif size(A,3) == 6
    AA = cat(1,A(:,:,1),A(:,:,2),A(:,:,3),A(:,:,4),A(:,:,5),A(:,:,6));
end


for ii=1:sz(1)
    for jj=1:sz(2)
        %        for kk = 1:sz(3)
        if mask(ii,jj,1) == 0 || imgs(ii,jj,1) ==0
            So(ii,jj)=0;
            So2(ii,jj) = 0;
            T2(ii,jj)=0;
            T2s(ii,jj) = 0;
            rmse(ii,jj)=0;
            continue
        end
        % fit voxel
        signal=squeeze(imgs(ii,jj,:));
        param = AA\log(signal);
        
        %  result
        So(ii,jj) = exp(param(1));
        So2(ii,jj) = exp(param(2));
        T2(ii,jj) = 1/param(3);
        T2s(ii,jj) = 1/param(4);
        rmse(ii,jj) = 100*(sqrt(mean(SAGE_cost([So(ii,jj) So2(ii,jj) T2(ii,jj) T2s(ii,jj)],signal,TEs,te2).^2))./mean(signal));
        
    end % end jj
end % end ii


So(isnan(So))   = 0;
So2(isnan(So2)) = 0;
T2(isnan(T2))   = 0;
T2s(isnan(T2s)) = 0;
rmse(isnan(rmse)) = 0;

