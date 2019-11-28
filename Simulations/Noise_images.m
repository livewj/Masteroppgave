% Generate your own MK (and other) map with various levels of noise

% Load data
% Subject 006
dwi = niftiread('/Users/live/Documents/006/eddy_corrected_data.nii');
dwi = double(dwi);
bvecs = importdata('/Users/live/Documents/006/eddy_corrected_data.eddy_rotated_bvecs');
bvals = importdata('/Users/live/Documents/006/bvals'); %17 b vals
grad = [bvecs' bvals'];
slice = 41;
mask = niftiread('/Users/live/Documents/006/brain_mask.nii'); % GM + WM + CSF
mask = boolean(mask(:,:,slice)); % Mask for selected slice
[nx, ny, nz, np] = size(dwi);
dwi_re = reshape(dwi, [nx*ny nz np]);
S_brain = squeeze(dwi_re(:, slice, :)); % Signal from selected slice

%% DKI parameter estimation of real signal

[~, dt_all] = dki_fit(dwi(:,:,slice,:), grad, mask); % b < 3000
[fa_all, md_all, rd_all, ad_all, ~, mk_all, rk_all, ak_all] = dki_parameters(dt_all);

% Real signal MK MAP
% Remove outliers
mk_all( mk_all > 2 ) = NaN; mk_all( mk_all < 0 ) = NaN;
figure, imagesc(mk_all), colormap gray, axis equal
w = 7; h = 5; set(gcf, 'PaperSize', [w h]); set(gcf, 'PaperPosition', [0 0 w h]);


%% Generate "noise-free" MK map

nvox = size(S_brain,1);

% Generate signal for each voxel
fa = NaN(106,106); md = NaN(106,106); rd = NaN(106,106); ad = NaN(106,106); 
mk = NaN(106,106); rk = NaN(106,106); ak = NaN(106,106); 

S_DKI = NaN(nvox, 138);

for v = 1:nvox
    if mask(v) == 0
        continue
    else
        S = S_brain(v,:);
        S = reshape(S, [1 1 1 length(S)]);
        [~, dt_real, w] = dki_fit(S, grad);
        S_DKI(v,:) = w; % Estimated signal

        % Estimate parameters
        S = S_DKI(v,:);
        S = reshape(S, [1 1 1 length(S)]);
        [~, dt] = dki_fit(S, grad(:,:));
        [fa(v), md(v), rd(v), ad(v), ~, mk(v), rk(v), ak(v)] = dki_parameters(dt);
    end
end

% Lag MK kart av generert "støyfritt" signal:
% Remove outliers
mk( mk > 2 ) = NaN; mk( mk < 0 ) = NaN;
figure, imagesc(mk), colormap gray, axis equal
figure, imagesc(fa), colormap gray, axis equal
w = 7; h = 5; set(gcf, 'PaperSize', [w h]); set(gcf, 'PaperPosition', [0 0 w h]);


% Normalize signals
inds0 = find(bvals < 10);
S_DKI_norm = S_DKI./nanmean(S_DKI(:,inds0),2);
S_brain_norm = S_brain./nanmean(S_brain(:,inds0),2);

% Plot signal from some arbitrary voxel (all dirs)
inds_voxel = 3433;
plot(bvals(inds3000),S_DKI_norm(inds_voxel,:),'b*'), hold on, ...
plot(bvals(inds3000),S_brain_norm(inds_voxel,inds3000), 'r*')

beep
%% ADD NOISE

n = 1000;
SNR = [5 10 20 40 60 80 100];


%Preallocate
fa_data = zeros(nvox, length(SNR), n);
md_data = zeros(nvox, length(SNR), n);
rd_data = zeros(nvox, length(SNR), n);
ad_data = zeros(nvox, length(SNR), n);
mk_data = zeros(nvox, length(SNR), n);
rk_data = zeros(nvox, length(SNR), n);
ak_data = zeros(nvox, length(SNR), n);

% MAIN LOOP
for v = 1:nvox % ~ 60 sec per voxel
    if mask(v) == 0
        continue
    else
        S = S_DKI_norm(v,:);
        for i = 1:n
            for j = 1:length(SNR)
                %Add rician noise
                N_real = randn(1,size(S,2))*(1/SNR(j));
                N_imag = randn(1,size(S,2))*(1/SNR(j));
                x = N_real + S;
                y = N_imag;
                S_noise = sqrt(x.^2 + y.^2);

                % Calculate tensor from signal with noise
                S_noise_reshaped = reshape(S_noise, [1, 1, 1, length(S_noise)]);
                [~, dt_app] = dki_fit(S_noise_reshaped, grad);
                dt_app = double(dt_app); 

                [fa_n, md_n, rd_n, ad_n, ~, mk_n, rk_n, ak_n] = dki_parameters(dt_app);

                fa_data(v,j,i) = fa_n;
                md_data(v,j,i) = md_n;
                rd_data(v,j,i) = rd_n;
                ad_data(v,j,i) = ad_n;
                mk_data(v,j,i) = mk_n;
                rk_data(v,j,i) = rk_n;
                ak_data(v,j,i) = ak_n;

            end
        end
    end
end

beep on; beep


% Draw random MK values for each voxel
mk_final = NaN(nvox, length(SNR));
fa_final = NaN(nvox, length(SNR));
for v = 1:nvox
    if mask(v) == 0
        continue
    else
        for j = 1:length(SNR)
            x = randi(n); % normalfordelt
            fa_final(v,j) = fa_data(v,j,x);
            mk_final(v,j) = mk_data(v,j,x);
        end
    end
end

% Reshape
fa_resh = reshape(fa_final, [106 106 length(SNR)]);
mk_resh = reshape(mk_final, [106 106 length(SNR)]);

% Parameter map with added noise
figure, imagesc(fa_resh(:,:,1)), colormap gray, axis equal
w = 7; h = 5; set(gcf, 'PaperSize', [w h]); set(gcf, 'PaperPosition', [0 0 w h]);

% Remove outliers
plot(mk_resh(:,:,1),'*')
mk_resh( mk_resh > 2 ) = NaN; mk_resh( mk_resh < 0 ) = NaN;
figure, imagesc(mk_resh(:,:,1)), colormap gray, axis equal
w = 7; h = 5; set(gcf, 'PaperSize', [w h]); set(gcf, 'PaperPosition', [0 0 w h]);


%% Save data

save('fa_final.mat','fa_final');
save('mk_final.mat','mk_final');

%%  Load data
A = load('fa_final.mat');
B = load('mk_final.mat');

fa_final = A.fa_final;
mk_final = B.mk_final;

%% Make histograms of FA and MK values in WM mask

nvox = 11236;

fa_final_wm = NaN(11236,7); mk_final_wm = NaN(11236,7);
for v = 1:nvox    
    for j = 1:length(SNR)
            if white_matter_mask(v) == 0
               continue
            else
                fa_final_wm(v,j) = fa_final(v,j);
                mk_final_wm(v,j) = mk_final(v,j);
        end
    end
end


% histogram(fa_final_wm(:,1),80), hold on, histogram(fa_final_wm(:,2),80), hold on, histogram(fa_final_wm(:,3),80)
% hold on, histogram(fa_final_wm(:,4),80), hold on, histogram(fa_final_wm(:,5),80)
% 
% histogram(mk_final_wm(:,1),80), hold on, histogram(mk_final_wm(:,2),80), hold on, histogram(mk_final_wm(:,3),80)
% hold on, histogram(mk_final_wm(:,4),80), hold on, histogram(mk_final_wm(:,5),80)


fontsize = 20;

%histogram(fa_final(:,1),80), hold on, histogram(fa_final(:,2),80), hold on, histogram(fa_final(:,3),80)
histogram(fa_final(:,4),80), hold on, histogram(fa_final(:,5),80), hold on, histogram(fa,80)
set(gca, 'fontsize', 20)
title('Histogram of FA values', 'fontsize', fontsize)
legend(['SNR = 40'], ['SNR = 60'], ['No added noise'])
xlabel('FA', 'fontsize', fontsize)
ylabel('Number of voxels', 'fontsize', fontsize)
w = 7; h = 5; set(gcf, 'PaperSize', [w h]); set(gcf, 'PaperPosition', [0 0 w h]);


%histogram(mk_final(:,1),80), hold on, histogram(mk_final(:,2),80), hold on, histogram(mk_final(:,3),80)
histogram(mk_final(:,4),80), hold on, histogram(mk_final(:,5),80), hold on, histogram(mk,80)
set(gca, 'fontsize', 20)
title('Histogram of MK values', 'fontsize', fontsize)
legend(['SNR = 40'], ['SNR = 60'], ['No added noise'])
xlabel('MK', 'fontsize', fontsize)
ylabel('Number of voxels', 'fontsize', fontsize)
w = 7; h = 5; set(gcf, 'PaperSize', [w h]); set(gcf, 'PaperPosition', [0 0 w h]);

