% Simulation experiment 1:
% Estimation of diffuison parameters at different levels of SNR
%% Load data

subject = '006'; % 002, 003 or 006

if subject == '002'
    dwi = niftiread('/Users/live/Documents/DTI_AP_002_corrected/eddy_corrected_data.nii');
    dwi = double(dwi);
    bvecs = importdata('/Users/live/Documents/DTI_AP_002_corrected/eddy_corrected_data.eddy_rotated_bvecs');
    bvals = importdata('/Users/live/Documents/DTI_AP/002/1_006_DTI_AP_20180131/002_bvals');
    grad = [bvecs' bvals'];
    slice = 41;
    mask_whole = niftiread('/Users/live/Documents/DTI_AP/002/1_006_DTI_AP_20180131/brain_b0_mask.nii'); % GM + WM + CSF
    mask_whole = boolean(mask_whole(:,:,slice));

    [nx, ny, nz, np] = size(dwi);
    dwi_re = reshape(dwi, [nx*ny nz np]);
    S_brain = squeeze(dwi_re(:, slice, :));
end

if subject == '003'
    dwi = niftiread('/Users/live/Documents/DTI_AP_003_corrected/eddy_corrected_data.nii'); 
    dwi = double(dwi);
    bvecs = importdata('/Users/live/Documents/DTI_AP_003_corrected/eddy_corrected_data.eddy_rotated_bvecs');
    bvals = importdata('/Users/live/Documents/DTI_AP_003/003/003.bval');
    grad = [bvecs' bvals'];
    
    % sort bvals in ascending order
    %[~,idx] = sort(grad(:,4)); % sort just the first column
    %sortedgrad = grad(idx,:);
    
    slice = 41;
    mask_whole = niftiread('/Users/live/Documents/DTI_AP_003/003/brain_b0_mask.nii'); % GM + WM + CSF
    mask_whole = boolean(mask_whole(:,:,slice));
    
    [nx, ny, nz, np] = size(dwi);
    dwi_re = reshape(dwi, [nx*ny nz np]);
    S_brain = squeeze(dwi_re(:, slice, :));
end

if subject == '006'
    dwi = niftiread('/Users/live/Documents/006/eddy_corrected_data.nii');
    dwi = double(dwi);
    bvecs = importdata('/Users/live/Documents/006/eddy_corrected_data.eddy_rotated_bvecs');
    bvals = importdata('/Users/live/Documents/006/bvals');
    grad = [bvecs' bvals'];
    slice = 44;
    mask_whole = niftiread('/Users/live/Documents/006/brain_mask.nii'); % GM + WM + CSF
    mask_whole = boolean(mask_whole(:,:,slice));

    [nx, ny, nz, np] = size(dwi);
    dwi_re = reshape(dwi, [nx*ny nz np]);
    S_brain = squeeze(dwi_re(:, slice, :));
end

%% ----- Noise-free parameter estimation based on real signal ------


[~, dt_all] = dki_fit(dwi(:,:,slice,:), grad, mask_whole);
[fa_all, md_all, rd_all, ad_all, ~, mk_all, rk_all, ak_all] = dki_parameters(dt_all);

% Generate parametric maps
%figure, imagesc(fa_all), colormap gray, axis equal
%figure, imagesc(md_all), colormap gray, axis equal
%figure, imagesc(rd_all), colormap gray, axis equal
%figure, imagesc(ad_all), colormap gray, axis equal

% mk_all( mk_all > 2 ) = NaN; % Remove outliers
% figure, imagesc(mk_all), colormap gray, axis equal
% figure, imagesc(rk_all), colormap gray, axis equal
% figure, imagesc(ak_all), colormap gray, axis equal

% Backround noise
%S_brain_b0 = dwi(:,:,slice,1);
%roi = roipoly(S_brain_b0/max(S_brain_b0(:))); 
%roi_std = std(S_brain_b0(roi))
%roi_mean = mean(S_brain_b0(roi))

beep on
beep
%% Experiment setup

NSA = 2;    % 1, 2

folder = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/subject' subject];

vox_nr = 3; % 1 = White matter
            % 2 = Gray matter
            % 3 = CSF   
           
            
if subject == '002'
    if vox_nr == 1 
        inds_voxel = 4189;
        tissue = 'white matter';
    end  
    if vox_nr == 2
        inds_voxel = 3657;
        tissue = 'gray matter';
    end  
    if vox_nr == 3
       inds_voxel = 4407;
       tissue = 'CSF';
    end    
    S_voxel = S_brain(inds_voxel,:);
end

if subject == '003'
    if vox_nr == 1 
        inds_voxel = 3983;
        tissue = 'white matter';
    end
    if vox_nr == 2
        inds_voxel = 3339;
        tissue = 'gray matter';
    end    
    if vox_nr == 3
        inds_voxel = 6730;
        tissue = 'CSF';
    end     
    S_voxel = S_brain(inds_voxel,:);
end

if subject == '006'
    if vox_nr == 1    
        inds_voxel = 3656;
        tissue = 'white matter';
    end   
    if vox_nr == 2
        inds_voxel = 2172;
        tissue = 'gray matter';
    end
    if vox_nr == 3
        inds_voxel = 4190;
        tissue = 'CSF';
    end    
    S_voxel = S_brain(inds_voxel,:);
end

disp('ok')

%% Generate signal and establish ground truth from selected voxel 

% Estimate signal from voxel
clear S_DKI
fa = []; md = []; rd = []; ad = []; mk = []; rk = []; ak = [];
for v = 1:size(S_voxel,1)
    S = S_voxel(v,:);
    S = reshape(S, [1 1 1 length(S)]);
    [~, dt_real, w] = dki_fit(S, grad); %full gradient set
    S_DKI(:,v) = w; % Estimated signal
    
    % Estimate parameters
    S = S_DKI(:,v)'; 
    S = reshape(S, [1 1 1 length(S)]);
    [~, dt] = dki_fit(S, grad); % b < 3000 
    [fa(v), md(v), rd(v), ad(v), ~, mk(v), rk(v), ak(v)] = dki_parameters(dt); % Ground truth
end


%Plot tensors
% DT = squeeze(dt(1:6));
% KT = squeeze(dt(7:21));
% addpath('/Users/live/Documents/Masterthesis/MATLAB/DKI/DKI_Estimation_Method/')
% figure(1), plotTensors(DT);
% figure(2), plotTensors(KT);

% Compare ground truth parameters with real-signal parameters
[fa_all(inds_voxel) fa; md_all(inds_voxel) md; rd_all(inds_voxel) rd; ad_all(inds_voxel) ad; ...
    mk_all(inds_voxel) mk; rk_all(inds_voxel) rk; ak_all(inds_voxel) ak]

% Print ground truth
[fa; md; rd; ad; mk; rk; ak]

% Normalize signal
bvals = grad(:,4);
%bvals = sortedgrad(:,4);
inds0 = find(bvals < 10);
S_DKI_norm = S_DKI./nanmean(S_DKI(inds0,:));
S_voxel_norm = S_voxel/mean(S_voxel(:,inds0));

% Compare real and estimated signal (all directions)
fontsize = 22;
inds3000 = find(bvals < 3010);
figure(3), plot(bvals(inds3000),S_voxel_norm(inds3000), 'r*'), hold on, plot(bvals(inds3000), S_DKI_norm,'b*')
set(gca,'fontsize', fontsize)
xlabel('b-value [s/mm^2]','fontsize', fontsize)
ylabel('Scaled signal intensity, S/S_0','fontsize', fontsize)
title(['Actual vs fitted signal, ' tissue],'fontsize', fontsize)
legend(['Actual signal'], ['Fitted signal'],'fontsize', fontsize)


% Compare real and estimated signal (one direction)
% inds_grad_dir = [1 7 25 50 89]; % 002
% inds_grad_dir = [1 7 25 51 89]; % 006
% inds_grad_dir = [2 3 16 47 88]; % 003
%figure(4), plot(bvals(inds_grad_dir), S_voxel_norm(inds_grad_dir), '-*'), hold on, ...
    %plot(bvals(inds_grad_dir), S_DKI_norm(inds_grad_dir),'-*')
% figure(4), plot(bvals(inds_grad_dir), S_voxel_norm_WM(inds_grad_dir), '-*'), hold on, ...
%     plot(bvals(inds_grad_dir), S_voxel_norm_GM(inds_grad_dir),'-*'), hold on, ...
%     plot(bvals(inds_grad_dir), S_voxel_norm_CSF(inds_grad_dir),'-*')


%% ADD NOISE

S = S_DKI_norm;
nvox=size(S,2);

n = 1000;
SNR = [10 15 20 25 30 40 50 60 80 100];

%Preallocate
fa_data = zeros(n, length(SNR), nvox);
md_data = zeros(n, length(SNR), nvox);
rd_data = zeros(n, length(SNR), nvox);
ad_data = zeros(n, length(SNR), nvox);
mk_data = zeros(n, length(SNR), nvox);
rk_data = zeros(n, length(SNR), nvox);
ak_data = zeros(n, length(SNR), nvox);

S_noise_data = zeros(length(S), n, length(SNR), NSA);

% MAIN LOOP
for v = 1:nvox % ~ 60 sec per voxel
    for k = 1:NSA 
        for i = 1:n
            for j = 1:length(SNR)
                %Add rician noise
                N_real = randn(1,size(S,1))*(1/SNR(j)); % randn()*sigma, sigma = noise = S0/SNR = 1/SNR
                N_imag = randn(1,size(S,1))*(1/SNR(j));
                x = N_real + S(:,v)';
                y = N_imag;
                S_noise = sqrt(x.^2 + y.^2);
                
                S_noise_data(:,i,j,k) = S_noise;
                
            end
        end
    end
end

% Average signals
S_noise_mean = mean(S_noise_data,4);


% Calculate parameters from averaged signal
for i = 1:n
    for j = 1:length(SNR)

        S_noise = S_noise_mean(:,i,j);

        % Calculate tensor from signal with noise
        S_noise_reshaped = reshape(S_noise, [1, 1, 1, length(S_noise)]);
        [~, dt_app] = dki_fit(S_noise_reshaped, grad(:,:));
        dt_app = double(dt_app); 

        [fa_n, md_n, rd_n, ad_n, ~, mk_n, rk_n, ak_n] = dki_parameters(dt_app);

        fa_data(i,j) = fa_n;
        md_data(i,j) = md_n;
        rd_data(i,j) = rd_n;
        ad_data(i,j) = ad_n;
        mk_data(i,j) = mk_n;
        rk_data(i,j) = rk_n;
        ak_data(i,j) = ak_n;
    end
end


% Remove outliers: Tukey's method
k = 1.5; % Outliers , k = 3; % Extreme outliers
fa_data( fa_data < (prctile(fa_data, 25) - k*iqr(fa_data)) ) = NaN;
fa_data( fa_data > (prctile(fa_data, 75) + k*iqr(fa_data)) ) = NaN;

md_data( md_data < (prctile(md_data, 25) - k*iqr(md_data)) ) = NaN;
md_data( md_data > (prctile(md_data, 75) + k*iqr(md_data)) ) = NaN;

rd_data( rd_data < (prctile(rd_data, 25) - k*iqr(rd_data)) ) = NaN;
rd_data( rd_data > (prctile(rd_data, 75) + k*iqr(rd_data)) ) = NaN;

ad_data( ad_data < (prctile(ad_data, 25) - k*iqr(ad_data)) ) = NaN;
ad_data( ad_data > (prctile(ad_data, 75) + k*iqr(ad_data)) ) = NaN;

mk_data( mk_data < (prctile(mk_data, 25) - k*iqr(mk_data)) ) = NaN;
mk_data( mk_data > (prctile(mk_data, 75) + k*iqr(mk_data)) ) = NaN;

rk_data( rk_data < (prctile(rk_data, 25) - k*iqr(rk_data)) ) = NaN;
rk_data( rk_data > (prctile(rk_data, 75) + k*iqr(rk_data)) ) = NaN;

ak_data( ak_data < (prctile(ak_data, 25) - k*iqr(ak_data)) ) = NaN;
ak_data( ak_data > (prctile(ak_data, 75) + k*iqr(ak_data)) ) = NaN;


% Average of all voxels 
% fa_app = mean(fa_data,4);
% md_app = mean(md_data,4);
% rd_app = mean(rd_data,4);
% ad_app = mean(ad_data,4);
% mk_app = mean(mk_data,4);
% rk_app = mean(rk_data,4); 
% ak_app = mean(ak_data,4);


beep on
beep
disp('SNR finished')

%% Compute relative error [%]

fa_rel = (fa_data - fa)./fa*100;
md_rel = (md_data - md)./md*100;
rd_rel = (rd_data - rd)./rd*100;
ad_rel = (ad_data - ad)./ad*100;
mk_rel = (mk_data - mk)./mk*100;
rk_rel = (rk_data - rk)./rk*100;
ak_rel = (ak_data - ak)./ak*100;

disp('Relative error computed')

%% SAVE RELATIVE ERROR TO FILE

voxel = ['voxel' num2str(vox_nr)];
filename = ['/subject' subject '_' voxel '_DKI_rel_nsa.txt'];
write_to_file_DKI(folder, filename, fa_rel, md_rel, rd_rel, ad_rel, mk_rel, rk_rel, ak_rel);

disp('Data saved')


%% Plot parameters

% ------ Boxplot --------

w = 7;
h = 5;

params = [{fa_rel} {md_rel} {rd_rel} {ad_rel} {mk_rel} {rk_rel} {ak_rel}];
params_string = ['FA'; 'MD'; 'RD'; 'AD'; 'MK'; 'RK'; 'AK'];

for i = 1:length(params)
    figure(i), clf, boxplot(params{i}, SNR)
    set(gca, 'fontsize', fontsize)
    hold on
    plot(linspace(0,105), 0*ones(1,100), '--');
    hold off
    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    %pause;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    %saveas(figure(i), fullfile(folder, ['subject' subject '_' voxel '_DKI_' params_string(i,:)]), 'pdf');
end

%% Pool subject data together

NSA = 2;
voxel = 'voxel3'; % 1, 2, 3

subject002_data_DKI = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/subject002/subject002_' voxel '_DKI_rel_nsa.txt']);
subject003_data_DKI = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/subject003/subject003_' voxel '_DKI_rel_nsa.txt']);
subject006_data_DKI = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/subject006/subject006_' voxel '_DKI_rel_nsa.txt']);

subject002_data_DTI = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/subject002/subject002_' voxel '_DTI2_rel_nsa.txt']);
subject003_data_DTI = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/subject003/subject003_' voxel '_DTI2_rel_nsa.txt']);
subject006_data_DTI = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/subject006/subject006_' voxel '_DTI2_rel_nsa.txt']);

fa_pool_DKI = [subject002_data_DKI(1:1000,:); subject003_data_DKI(1:1000,:); subject006_data_DKI(1:1000,:)];
md_pool_DKI = [subject002_data_DKI(1001:2000,:); subject003_data_DKI(1001:2000,:); subject006_data_DKI(1001:2000,:)];
rd_pool_DKI = [subject002_data_DKI(2001:3000,:); subject003_data_DKI(2001:3000,:); subject006_data_DKI(2001:3000,:)];
ad_pool_DKI = [subject002_data_DKI(3001:4000,:); subject003_data_DKI(3001:4000,:); subject006_data_DKI(3001:4000,:)];
mk_pool_DKI = [subject002_data_DKI(4001:5000,:); subject003_data_DKI(4001:5000,:); subject006_data_DKI(4001:5000,:)];
rk_pool_DKI = [subject002_data_DKI(5001:6000,:); subject003_data_DKI(5001:6000,:); subject006_data_DKI(5001:6000,:)];
ak_pool_DKI = [subject002_data_DKI(6001:7000,:); subject003_data_DKI(6001:7000,:); subject006_data_DKI(6001:7000,:)];

fa_pool_DTI = [subject002_data_DTI(1:1000,:); subject003_data_DTI(1:1000,:); subject006_data_DTI(1:1000,:)];
md_pool_DTI = [subject002_data_DTI(1001:2000,:); subject003_data_DTI(1001:2000,:); subject006_data_DTI(1001:2000,:)];
rd_pool_DTI = [subject002_data_DTI(2001:3000,:); subject003_data_DTI(2001:3000,:); subject006_data_DTI(2001:3000,:)];
ad_pool_DTI = [subject002_data_DTI(3001:4000,:); subject003_data_DTI(3001:4000,:); subject006_data_DTI(3001:4000,:)];


folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/Pool/'];

% SAVE TO FILE
filename_pool_DKI = ['/' voxel '_DKI_22.txt'];
filename_pool_DTI = ['/' voxel '_DTI_22.txt'];

write_to_file_DKI(folder_pool, filename_pool_DKI, fa_pool_DKI, md_pool_DKI, rd_pool_DKI, ad_pool_DKI, mk_pool_DKI, rk_pool_DKI, ak_pool_DKI);
write_to_file_DTI(folder_pool, filename_pool_DTI, fa_pool_DTI, md_pool_DTI, rd_pool_DTI, ad_pool_DTI);

disp('ok')

%% Plot pooled parameters

w = 7;
h = 5;

% params_pool = [{fa_pool_DKI} {md_pool_DKI} {rd_pool_DKI} {ad_pool_DKI} {mk_pool_DKI} {rk_pool_DKI} {ak_pool_DKI}];
% params_string = ['FA'; 'MD'; 'RD'; 'AD'; 'MK'; 'RK'; 'AK'];
% model = 'DKI';

params_pool = [{fa_pool_DTI} {md_pool_DTI} {rd_pool_DTI} {ad_pool_DTI}];
params_string = ['FA'; 'MD'; 'RD'; 'AD'];
model = 'DTI';


for i = 1:length(params_pool)
    figure(i), clf, boxplot(params_pool{i}, SNR)
    set(gca, 'fontsize', fontsize)
    hold on
    plot(linspace(0,105), 0*ones(1,100), '--');
    hold off
    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    %pause;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    %saveas(figure(i), fullfile(folder_pool, [voxel '_' model '_' params_string(i,:)]), 'pdf');
end
