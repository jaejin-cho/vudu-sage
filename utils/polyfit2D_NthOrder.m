function I_fitted = polyfit2D_NthOrder(I,mask,Order)
%polyfit the part of the image 'I' that is inside the mask

Ivec = I(mask==1);

matrix(:, 1) = ones(sum(mask(:)), 1);
FullMatrix(:, 1) = double(ones(size(mask,1)*size(mask,2),1));

x = (0 : size(mask,2) - 1)/size(mask,2) - 1/2;
Xmat = repmat(x, size(mask, 1), 1);
y = ((0 : size(mask,1) - 1)/size(mask,1) - 1/2).';
Ymat = repmat(y, 1, size(mask, 2));

for CurrentOrder = 1:Order
    for Xpower = 0:CurrentOrder
        matrix = [matrix, Xmat(mask == 1).^Xpower .*Ymat(mask == 1).^(CurrentOrder-Xpower)];
        FullMatrix = [FullMatrix, Xmat(:).^Xpower .*Ymat(:).^(CurrentOrder-Xpower)];
    end
end

coefs = matrix \ Ivec; % determine coefficients based on data in W2
I_fittedVec = FullMatrix * coefs;
%I_fitted = reshape(I_fittedVec,size(I,1),size(I,2)).*mask;
I_fitted = reshape(I_fittedVec,size(I,1),size(I,2));

