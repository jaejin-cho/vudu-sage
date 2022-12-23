%--------------------------------------------------------------------------
%% load data
%--------------------------------------------------------------------------
load("data.mat");
addpath ('utils');
[nx,ny,nz,nc,ne] = size(kspace);

%--------------------------------------------------------------------------
%% Echo spacing time
%--------------------------------------------------------------------------
esp     =   prot.iEffectiveEpiEchoSpacing * 1e-6;          
time_ap	=   vec(linspace(0,ny-1,ny) * esp);
time_pa	=   flip(time_ap,1);

%--------------------------------------------------------------------------
%% LORAKS param
%--------------------------------------------------------------------------
FTy         =   fftshift(fft(eye(ny,ny),[],2)) ;
LORAKS_type =   1;
R           =   5;
rank        =   250;
tol         =   5e-3;
lambda      =   1e-2;
max_iter    =   50;
[in1,in2]   =   meshgrid(-R:R,-R:R);
idx         =   find(in1.^2+in2.^2<=R^2);
patchSize   =   numel(idx);

%--------------------------------------------------------------------------
%% buda reconstruction
%--------------------------------------------------------------------------
h           =   mifft(kspace,1) .* sqrt(nx);   
im_recon	=   zeross([nx,ny,nz,ne]);

for zz  =   1:nz    
    %--------------------------------------------------------------------------
    %% data
    %--------------------------------------------------------------------------
    hslice      =   h(:,:,zz,:,:);
    hmask       =   sq(abs(hslice) > 0);
    cslice      =   sq(csm(:,:,zz,:));
    
    %--------------------------------------------------------------------------
    %% FT matrix
    %--------------------------------------------------------------------------
    FTappa      =   zeross([ny,ny,nx,ne]);
    for xx = 1:nx
        fe  =   vec(field(xx,:,zz));
        for sh = 1:ne
            if mod(sh,2)
                phs =   exp( 1i .* 2 .* pi .* (time_ap * fe.'));
            else
                phs =   exp( 1i .* 2 .* pi .* (time_pa * fe.'));
            end
            m       =   find(sum(hmask(xx,:,:,sh),3));
            FTappa(m,:,xx,sh)	=	FTy(m,:) .* phs(m,:);
        end
    end
    
    %--------------------------------------------------------------------------
    %% define function
    %--------------------------------------------------------------------------
    A1          =   @(x) vec(hmask.*ft2_ms_epi_b0( repmat(reshape(x,[nx,ny,1,ne]),[1,1,nc,1]) .* repmat(cslice,[1,1,1,ne]), FTappa));
    A2          =   @(x) vec(sum(ift2_ms_epi_b0(hmask.*reshape(x,[nx,ny,nc,ne]), FTappa) .* repmat(conj(cslice),[1,1,1,ne]),3));
    A3          =   @(x) vec(A2(A1(x)));
    
    Aty         =   reshape(A2(sq(hslice)),[nx,ny,ne]);
    [x_out,~]   =   pcg(@(x,flag) A3(x), vec(Aty) );
    im_recon(:,:,zz,:)     =   reshape(x_out,[nx,ny,ne]);
    
    %--------------------------------------------------------------------------
    %% loraks functions
    %--------------------------------------------------------------------------
    B1  =   @(x) vec(mfft2(reshape(x, [nx ny ne])));
    B2  =   @(x) vec(mifft2(reshape(x,[nx ny ne])));
    P1  =   @(x) LORAKS_operators(x,nx,ny,ne,R,LORAKS_type,[]);
    P2  =   @(x) LORAKS_operators(x,nx,ny,ne,R,-LORAKS_type,[]);
    Z1  =   @(x) padarray(reshape(x,[nx,ny,ne]),[2*R, 2*R], 'post');
    Z2  =   @(x) x(1:nx,1:ny,:,:);
    
    %--------------------------------------------------------------------------
    %% loop - iteration
    %--------------------------------------------------------------------------
    z   =   vec(im_recon(:,:,zz,:));
    
    for iter = 1:max_iter
        z_cur   =   z;
        pz      =   B1(z_cur);
        MM      =   P1(pz);
        
        Um      =   svd_left(MM);
        nmm     =   Um(:,rank+1:end)'; % null space
        Bhr     =   0;
        
        if LORAKS_type == 1 % S
            nf      = size(nmm,1);
            nmm     = reshape(nmm,[nf, patchSize, 2*ne]);
            nss_h   = reshape(nmm(:,:,1:2:end)+1j*nmm(:,:,2:2:end),[nf, patchSize*ne]);
            Nis     = filtfilt(nss_h,'C',nx,ny,ne,R);
            Nis2    = filtfilt(nss_h,'S',nx,ny,ne,R);
            LhL     = @(x) 2*B2((Z2(ifft2(squeeze(sum(Nis.*repmat(fft2(Z1(B1(x))),[1 1 1 ne]),3))))) ...
                -(Z2(ifft2(squeeze(sum(Nis2.*repmat(conj(fft2(Z1(B1(x)))),[1 1 1 ne]),3))))));
        end
        
        % data fitting
        M       =   @(x) A3(x) + lambda*LhL(x);
        [z,~]   =   pcg(M, vec(Aty) + lambda*Bhr, [], 10, [], [], z_cur);
        
        t       =   (norm(z_cur-z)/norm(z));
        % display the status
        if ~rem(iter,1)
            disp(['iter ' num2str(iter) ', relative change in solution: ' num2str(t)]);
        end
        im_recon(:,:,zz,:) = reshape(z,[nx,ny,1,ne]);
        % check for convergence
        if t < tol
            disp('Convergence tolerance met: change in solution is small');
            break;
        end
    end    
end

mosaic(abs(sq(rot90(cat(4,im_recon(:,:,:,1:2:end),im_recon(:,:,:,2:2:end))))),2,5,11,'vudu sage recon AP/PA', [0,1e-3]);

save('result.mat','im_recon','-v7.3');

%--------------------------------------------------------------------------
%% Fitting
%--------------------------------------------------------------------------

T2_4p_fit   =   zeross([nx,ny,nz,2]);
T2_3p_fit   =   T2_4p_fit;
fit_iter    =   10;
net         =   denoisingNetwork('dncnn');
zz          =   1;

%--------------------------------------------------------------------------
% fitting param
%--------------------------------------------------------------------------
imgs_use = sq(rsos(reshape(im_recon(:,:,zz,:),[nx,ny,2,5]),3))*5e2;

% mask
mask = rsos(imgs_use,3);
mask = mask / max(mask(:));
mask = bwareaopen( imfill( mask > 0.1, 'holes' ), 1500, 4);

%--------------------------------------------------------------------------
% 4-parameter fit
%--------------------------------------------------------------------------

[So_4p, So2_4p, T2_4p, T2s_4p,  rmse] = SAGE_fast(imgs_use, TEs, te2, mask, 0);

T2_4p_fit(:,:,zz,1) = T2s_4p;
T2_4p_fit(:,:,zz,2) = T2_4p;


%--------------------------------------------------------------------------
% denoise
%--------------------------------------------------------------------------

for ec = 1:5
    imgs_use(:,:,ec) = denoiseImage(imgs_use(:,:,ec), net);
end

%--------------------------------------------------------------------------
% 3-parameter fit
%--------------------------------------------------------------------------

So_3p   = So_4p;
So2_3p  = So2_4p;
T2_3p   = T2_4p;
T2s_3p  = T2s_4p;

for n_iter = 1:fit_iter
    beta0 = cat(3, So_3p, abs(T2_3p), abs(T2s_3p));
    delta = polyfit2D_NthOrder( abs(So_3p ./  So2_3p), mask, 2 );
    [So_3p, T2_3p, T2s_3p] = T2T2sprocess_3p_limited_test( double(imgs_use), TEs, te2, delta, beta0, mask);
    So2_3p = So_3p ./ delta;
end

T2_3p_fit(:,:,zz,1) = T2s_3p;
T2_3p_fit(:,:,zz,2) = T2_3p;

%--------------------------------------------------------------------------
%% delta analysis (mask region, should select one slice)
%--------------------------------------------------------------------------
delta_max = max(max(delta(mask==1)));
delta_min = min(min(delta(mask==1)));
mean(delta(mask==1));               

num_dict	= 100;                  % how many kinds of dict you want to build
delta_msk   = delta .* mask;
delta_value = linspace(delta_min, delta_max, num_dict); 
step        = ( delta_value(2) - delta_value(1) ) / 2;

for i = 1: size(delta_value,2)
    delta_msk(delta_msk <= (delta_value(i)+step)  &  delta_msk > (delta_value(i)-step) ) = delta_value(i);    
end
mask_roi = zeross([size(delta_msk,1), size(delta_msk,2), size(delta_value,2)]);
for i = 1: size(delta_value,2)
    mask_roi(:,:,i) = delta_msk == delta_value(i);
end

%--------------------------------------------------------------------------
%% dictionary matching
%--------------------------------------------------------------------------
N       =   0;
T2      =   [1:1:50 52:2:150 155:5:500]; 
T2s     =   [1:1:50 52:2:150 155:5:300]; 
T2_sz   =   size(T2,2);
T2s_sz  =   size(T2s,2);
delta_value_sz  = size(delta_value,2);
TEs_sz          = size(TEs,2);
for ii = 1:T2_sz
    for jj = 1:T2s_sz
        if (T2(ii) > T2s(jj)) && (T2(ii)/T2s(jj) < 4)
            N = N+1;
        end            
    end
end 
signal          =   zeross([N, length(TEs), delta_value_sz]); 
lookup_table    =   zeross([N, 2, delta_value_sz]);
N               =   0;

for ii = 1:T2_sz
    for jj = 1:T2s_sz            
        if (T2(ii) > T2s(jj)) && (T2(ii)/T2s(jj) < 4)
            N=N+1;                
            for kk = 1:delta_value_sz                
                signal1 = delta_value(kk) * exp(-TEs(1)./T2s(jj));
                signal2 = delta_value(kk) * exp(-TEs(2)./T2s(jj));                              
                signal3 = exp(-TEs(5).*(1./T2s(jj)-1./T2(ii))) .* exp(-TEs(3).*(2./T2(ii)-1./T2s(jj)));
                signal4 = exp(-TEs(5).*(1./T2s(jj)-1./T2(ii))) .* exp(-TEs(4).*(2./T2(ii)-1./T2s(jj)));
                signal5 = exp(-TEs(5).*(1./T2s(jj)-1./T2(ii))) .* exp(-TEs(5).*(2./T2(ii)-1./T2s(jj)));       
                signal(N, :, kk) = cat(2,signal1,signal2,signal3,signal4,signal5);
                lookup_table(N, :, kk) = [T2(ii) T2s(jj)];                
            end
        end        
    end
end

% normalize dict.signal
for m = 1:size(signal,1)
    for n = 1:size(signal,3)
        signal(m,:,n) = signal(m,:,n)/sum(abs(signal(m,:,n)).^2)^0.5;
    end
end

%--------------------------------------------------------------------------
% Fast Match
%--------------------------------------------------------------------------

sz          =   size(imgs_use);
T2_2par     =   zeross(sz(1:2));
T2s_2par    =   zeross(sz(1:2)); 
num_voxel   =   zeross([1, delta_value_sz]); 

for i = 1: delta_value_sz    
    num_voxel(i)    =   sum(sum(mask_roi(:,:,i)), 2);    
    imgs_voxel      =   zeross([length(TEs), num_voxel(i)]);    
    for t = 1:length(TEs)
        temp = imgs_use(:,:,t);
        imgs_voxel(t,:) = temp(mask_roi(:,:,i)==1);
    end
    res     =   signal(:,:,i) * imgs_voxel;
    [~, max_idx] = max(abs(res), [], 1);
    res_map =   lookup_table(max_idx,:,i);    
    T2_2par(mask_roi(:,:,i)==1)     =   res_map(:,1);
    T2s_2par(mask_roi(:,:,i)==1)    =   res_map(:,2);
end

mosaic(abs(sq(rot90(T2s_2par))),1,1,12,'est T_2^*', [0,150]); colormap hot;  
mosaic(abs(sq(rot90( T2_2par))),1,1,13,'est T_2 ',  [0,200]); colormap hot;  

save('result.mat','T2s_2par','T2_2par','-append');
