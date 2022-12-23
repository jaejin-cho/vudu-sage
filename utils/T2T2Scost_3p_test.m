function cost = T2T2Scost_3p_test(params,S,TEs,te2,delta)
S0 = params(1);
S02 = S0./delta;
T2 = params(2);
T2s = params(3);
Sfit = zeros(1,size(TEs,1),size(TEs,2));


for i = 1:size(TEs,2)
ind = find(TEs(:,i)>te2(:,i),1);
Sfit(1,1:ind-1,i) = S0.*exp(-TEs(1:ind-1,i)./T2s);
%Sfit(ind:end) = S02.*exp(-te2.*2.*(1./T2s-1./T2))...
%                .*exp(-te(ind:end).*(2./T2-1./T2s));


ind2 = find(TEs(:,i)>te2(:,i)*2,1); 
if (isempty(ind2)) 
    ind2 = size(TEs,1)+1;
end
if ind2 == ind
    ind2 = ind2+1;
end
Sfit(1,ind:ind2-1,i) = S02.*exp(-te2(:,i).*2.*(1./T2s-1./T2))...
                .*exp(-TEs(ind:ind2-1,i).*(2./T2-1./T2s));
Sfit(1,ind2:end,i) = S02.*exp(-te2(:,i).*2.*(1./T2-1./T2s))...
                .*exp(-TEs(ind2:end,i).*(1./T2s));
            
end


if size(Sfit,3) == 1
    Sfitting = cat(2,Sfit(:,:,1));
elseif size(Sfit,3) == 2
    Sfitting = cat(2,Sfit(:,:,1),Sfit(:,:,2));
elseif size(Sfit,3) == 3
    Sfitting = cat(2,Sfit(:,:,1),Sfit(:,:,2),Sfit(:,:,3)); % depend on the number of groups 
elseif size(Sfit,3) == 4
    Sfitting = cat(2,Sfit(:,:,1),Sfit(:,:,2),Sfit(:,:,3),Sfit(:,:,4));
elseif size(Sfit,3) == 5
    Sfitting = cat(2,Sfit(:,:,1),Sfit(:,:,2),Sfit(:,:,3),Sfit(:,:,4),Sfit(:,:,5));
end

cost = abs(Sfitting')-S;
cost = abs(abs(Sfitting')- abs(S));





%cost = sum(cost.^2); % lsqnonlin wants vector
% only use sos for fminsearch
