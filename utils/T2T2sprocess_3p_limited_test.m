function [So T2 T2s rmse] = T2T2sprocess_3p_limited_test(imgs,TEs,te2,delta,beta0,mask)

displayON = 0;
%options=optimset('Display','off','MaxIter',100,'TolFun',1e-4);
options=optimset('Display','off','MaxIter',100);
sz=size(imgs);

%ind = find(TEs>te2,1);
%ind2 = find(tout>te2*2,1);                   
%if (isempty(ind2)) 
%  ind2 = length(te);
%end

%A1 = [ones(ind-1,1), zeros(ind-1,1), zeros(ind-1,1),-TEs(1:ind-1)];
%A2 = [zeros(length(TEs)-ind+1,1), ones(length(TEs)-ind+1,1), (te2-2*TEs(ind:end)), (TEs(ind:end)-te2)];
%A = [A1; A2];

for ii=1:sz(1)
    %ii
    for jj=1:sz(2)
%        for kk = 1:sz(3)
            if mask(ii,jj,1) == 0 || imgs(ii,jj,1) ==0 
                So(ii,jj)=0;
                T2(ii,jj)=0;
                T2s(ii,jj) = 0;
                rmse(ii,jj)=0;
                continue
            end
           % fit voxel
           signal=squeeze(imgs(ii,jj,:));
           beta0_vox = squeeze(beta0(ii,jj,:));
           maxlim_vox = beta0_vox*500;
           minlim_vox = beta0_vox*.002;
          % [param] = lsqnonlin(@T2T2Scost_3p,beta0_vox,minlim_vox,maxlim_vox,options, signal, TEs, te2, delta(ii,jj));
           [param, resnorm, residual, exitflag, output] = lsqnonlin(@T2T2Scost_3p_test,beta0_vox,minlim_vox,maxlim_vox,options, signal, TEs, te2, delta(ii,jj));
          % length(residual)
           %param = A\log(signal);
           %  result  
            %So(ii,jj) = exp(param(1));
            %So2(ii,jj) = exp(param(2));
            %T2(ii,jj) = 1/param(3);
            %T2s(ii,jj) = 1/param(4);
            So(ii,jj) = param(1);
            T2(ii,jj) = param(2);
            T2s(ii,jj) = param(3);
           % rmse(ii,jj) = resnorm(1);
           % rmse(ii,jj) = sqrt(mean(residual.^2));
          %  rmse(ii,jj) =1;
            rmse(ii,jj) = 100*(sqrt(mean(T2T2Scost_3p_test([So(ii,jj)  T2(ii,jj) T2s(ii,jj)],signal,TEs,te2,delta(ii,jj)).^2))./mean(signal));

%             if displayON && ii==floor(sz(1)/2) && jj==floor(sz(2)/2) %&& kk==25
%                     tout = TEs(1):1:TEs(end);
%                     yout = zeros(1,length(tout));
%                     ind = find(tout>te2,1);
%                     yout(1:ind-1) = So(ii,jj).*exp(-tout(1:ind-1)./T2s(ii,jj));
%                     %yout(ind:end) = So2(ii,jj).*exp(-te2.*2.*(1./T2s(ii,jj)-1./T2(ii,jj)))...
%                     %                .*exp(-tout(ind:end).*(2./T2(ii,jj)-1./T2s(ii,jj)));
%                     ind2 = find(tout>te2*2,1);
%                     if (isempty(ind2)) 
%                         ind2 = length(tout);
%                     end
%                     yout(ind:ind2-1) = So(ii,jj)./delta(ii,jj).*exp(-te2.*2.*(1./T2s(ii,jj)-1./T2(ii,jj)))...
%                                     .*exp(-tout(ind:ind2-1).*(2./T2(ii,jj)-1./T2s(ii,jj)));
%                     yout(ind2:end) = So(ii,jj)./delta(ii,jj).*exp(-te2.*2.*(1./T2(ii,jj)-1./T2s(ii,jj)))...
%                                     .*exp(-tout(ind2:end).*(1./T2s(ii,jj)));
% 
%                                 
%                     figure; plot(TEs,signal,'o','LineWidth',2); 
%                     hold on; plot(tout,yout,'LineWidth',2); 
%                     title(['T2 and T2* fit, RMSE = ', num2str(rmse(ii,jj),4)]);
%                     text(TEs(2),max(yout(:))*.9,['So=',num2str(So(ii,jj),2),  ...
%                           '    T2=',num2str(T2(ii,jj),4), '    T2*=',num2str(T2s(ii,jj),4)]); 
%                     drawnow;
%             end % end display 
%        end % end kk 
    end % end jj 
end % end ii