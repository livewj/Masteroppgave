% Simulation experiment 3:
% Evaluating the influence of the number of b-values and gradient
% directions on the parameter estimations
% SNR = constant

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
%     [~,idx] = sort(grad(:,4)); % sort just the first column
%     sortedgrad = grad(idx,:);
    
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

%% Create gradient sets

inds_all = find(bvals < 3000);

inds0 = find(bvals<10);
inds500 = find(400<bvals & bvals<600);
inds1000 = find(800<bvals & bvals<1200);
inds2000 = find(1900<bvals & bvals<2100);
inds3000 = find(2900<bvals & bvals<3100);

% DKI gradient sets
set1 = [inds0 inds500 inds1000]; % 6 + 12 + 30
set2 = [inds0(1) inds500(1:2:end) inds3000(1:2:end)]; % 1 + 6 + 25
set3 = [inds0(1) inds500 inds3000]; % 1 + 12 + 50
set4 = [inds0 inds500 inds3000];    % 6 + 12 + 50
set5 = [inds0(1) inds500(1:2:end) inds1000(1:2:end) inds3000(1:2:end)]; % 1 + 6 + 15 + 25
set6 = [inds0(1) inds500 inds1000 inds3000]; % 1 + 12 + 30 + 50
set7 = [inds0 inds500 inds1000 inds3000];    % 6 + 12 + 30 + 50
set8 = [inds0 inds500 inds1000 inds2000 inds3000]; % 6 + 12 + 30 + 40 + 50


% Plot gradient directions
% inds = inds3000(1:2:end);
% figure, hold on, 
% for i=1:length(inds)
%     plot3([-grad(inds(i),1) grad(inds(i),1)], [-grad(inds(i),2) grad(inds(i),2)], [-grad(inds(i),3) grad(inds(i),3)]);
%     %pause;
%     axis equal;
%      set(gca,'FontSize',14)
%      xlabel('x','FontSize',14)
%      ylabel('y','FontSize',14)
%      zlabel('z','FontSize',14)
%      title('Gradient directions','FontSize',14)
% end


disp('ok')
%% Noise-free parameter estimation based on real signal

[~, dt_all] = dki_fit(dwi(:,:,slice,:), grad, mask_whole);
[fa_all, md_all, rd_all, ad_all, ~, mk_all, rk_all, ak_all] = dki_parameters(dt_all);


%% Experiment setup

NSA = 1;


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

% Estimate signal from voxel (Full set)
clear S_DKI
fa = []; md = []; rd = []; ad = []; mk = []; rk = []; ak = [];
for v = 1:size(S_voxel,1)
    S = S_voxel(v,:);
    S = reshape(S, [1 1 1 length(S)]);
    [~, ~, w] = dki_fit(S, grad(:,:));
    S_DKI(:,v) = w; % Estimated signal
    
    % Estimate parameters
    S = S_DKI(:,v)'; 
    S = reshape(S, [1 1 1 length(S)]);
    [~, dt] = dki_fit(S, grad(:,:));
    [fa(v), md(v), rd(v), ad(v), ~, mk(v), rk(v), ak(v)] = dki_parameters(dt); % Ground truth
end

% Compare estimated (ground truth) parameters with real-signal parameters
[fa_all(inds_voxel) fa; md_all(inds_voxel) md; rd_all(inds_voxel) rd; ad_all(inds_voxel) ad; ...
    mk_all(inds_voxel) mk; rk_all(inds_voxel) rk; ak_all(inds_voxel) ak]

% Print ground truth
[fa; md; rd; ad; mk; rk; ak];


% Normalize signal
bvals = grad(:,4);
inds0 = find(bvals < 10);
S_DKI_norm = S_DKI./nanmean(S_DKI(inds0,:));


%% ADD NOISE

nvox=size(S_DKI_norm,2);

n = 1000;
SNR = 50;

set_nr = [{set1} {set2} {set3} {set4} {set5} {set6} {set7} {set8}];


%Preallocate
fa_data = zeros(n, length(SNR), NSA, nvox);
md_data = zeros(n, length(SNR), NSA, nvox);
rd_data = zeros(n, length(SNR), NSA, nvox);
ad_data = zeros(n, length(SNR), NSA, nvox);
mk_data = zeros(n, length(SNR), NSA, nvox);
rk_data = zeros(n, length(SNR), NSA, nvox);
ak_data = zeros(n, length(SNR), NSA, nvox);

% MAIN LOOP
for s = 1:length(set_nr)
    inds = set_nr{s};
    S = S_DKI_norm(inds);
    s
    for v = 1:nvox % ~ 60 sec per voxel
        for k = 1:NSA 
            for i = 1:n
                for j = 1:length(SNR)
                    %Add rician noise
                    N_real = randn(1,size(S,1))*(1/SNR(j));
                    N_imag = randn(1,size(S,1))*(1/SNR(j));
                    x = N_real + S(:,v)';
                    y = N_imag;
                    S_noise = sqrt(x.^2 + y.^2);

                    % Calculate tensor from signal with noise
                    S_noise_reshaped = reshape(S_noise, [1, 1, 1, length(S_noise)]);
                    [~, dt_app] = dki_fit(S_noise_reshaped, grad(inds,:));
                    dt_app = double(dt_app); 

                    [fa_n, md_n, rd_n, ad_n, ~, mk_n, rk_n, ak_n] = dki_parameters(dt_app);

                    fa_data(i,j,k,v) = fa_n;
                    md_data(i,j,k,v) = md_n;
                    rd_data(i,j,k,v) = rd_n;
                    ad_data(i,j,k,v) = ad_n;
                    mk_data(i,j,k,v) = mk_n;
                    rk_data(i,j,k,v) = rk_n;
                    ak_data(i,j,k,v) = ak_n;
                    
                    fa_set{s} = fa_data;
                    md_set{s} = md_data;
                    rd_set{s} = rd_data;
                    ad_set{s} = ad_data;
                    mk_set{s} = mk_data;
                    rk_set{s} = rk_data;
                    ak_set{s} = ak_data;
                 

                end
            end
        end
    end
end


% FOR ALL SETS:

% Remove outliers: Tukey's method
k = 1.5;
for s = 1:length(set_nr)
    fa_set{s}( fa_set{s} < (prctile(fa_set{s}, 25) - k*iqr(fa_set{s})) ) = NaN;
    fa_set{s}( fa_set{s} > (prctile(fa_set{s}, 75) + k*iqr(fa_set{s})) ) = NaN;

    md_set{s}( md_set{s} < (prctile(md_set{s}, 25) - k*iqr(md_set{s})) ) = NaN;
    md_set{s}( md_set{s} > (prctile(md_set{s}, 75) + k*iqr(md_set{s})) ) = NaN;
    
    rd_set{s}( rd_set{s} < (prctile(rd_set{s}, 25) - k*iqr(rd_set{s})) ) = NaN;
    rd_set{s}( rd_set{s} > (prctile(rd_set{s}, 75) + k*iqr(rd_set{s})) ) = NaN;
    
    ad_set{s}( ad_set{s} < (prctile(ad_set{s}, 25) - k*iqr(ad_set{s})) ) = NaN;
    ad_set{s}( ad_set{s} > (prctile(ad_set{s}, 75) + k*iqr(ad_set{s})) ) = NaN;
    
    mk_set{s}( mk_set{s} < (prctile(mk_set{s}, 25) - k*iqr(mk_set{s})) ) = NaN;
    mk_set{s}( mk_set{s} > (prctile(mk_set{s}, 75) + k*iqr(mk_set{s})) ) = NaN;
    
    rk_set{s}( rk_set{s} < (prctile(rk_set{s}, 25) - k*iqr(rk_set{s})) ) = NaN;
    rk_set{s}( rk_set{s} > (prctile(rk_set{s}, 75) + k*iqr(rk_set{s})) ) = NaN;
    
    ak_set{s}( ak_set{s} < (prctile(ak_set{s}, 25) - k*iqr(ak_set{s})) ) = NaN;
    ak_set{s}( ak_set{s} > (prctile(ak_set{s}, 75) + k*iqr(ak_set{s})) ) = NaN;
end


% Average value from all NSA
clear fa_nsa; clear md_nsa; clear rd_nsa; clear ad_nsa; clear mk_nsa; clear rk_nsa, clear ak_nsa;
for s = 1:length(set_nr)
    fa_nsa{s} = nanmean(fa_set{s},3);
    md_nsa{s} = nanmean(md_set{s},3);
    rd_nsa{s} = nanmean(rd_set{s},3);
    ad_nsa{s} = nanmean(ad_set{s},3);
    mk_nsa{s} = nanmean(mk_set{s},3);
    rk_nsa{s} = nanmean(rk_set{s},3);
    ak_nsa{s} = nanmean(ak_set{s},3);
end

% Average of all voxels 
clear fa_app, clear md_app, clear rd_app; clear ad_app; clear mk_app, clear rk_app, clear ak_app;
for s = 1:length(set_nr)
    fa_app{s} = mean(fa_nsa{s},4);
    md_app{s} = mean(md_nsa{s},4);
    rd_app{s} = mean(rd_nsa{s},4);
    ad_app{s} = mean(ad_nsa{s},4);
    mk_app{s} = mean(mk_nsa{s},4);
    rk_app{s} = mean(rk_nsa{s},4); 
    ak_app{s} = mean(ak_nsa{s},4);
end

disp('SNR finished')

% Compute relative error [%]
clear md_rel, clear rd_rel; clear ad_rel; clear fa_rel, clear mk_rel, clear rk_rel, clear ak_rel
for s = 1:length(set_nr)
    fa_rel{s} = (fa_app{s} - fa)./fa*100;
    md_rel{s} = (md_app{s} - md)./md*100;
    rd_rel{s} = (rd_app{s} - rd)./rd*100;
    ad_rel{s} = (ad_app{s} - ad)./ad*100;
    mk_rel{s} = (mk_app{s} - mk)./mk*100;
    rk_rel{s} = (rk_app{s} - rk)./rk*100;
    ak_rel{s} = (ak_app{s} - ak)./ak*100;
end


disp('Relative error computed')

%% SAVE RELATIVE ERROR TO FILE

voxel = ['voxel' num2str(vox_nr)];
folder = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/subject' subject '/' voxel];


for s = 1:length(set_nr)
    filename = ['/subject' subject '_' voxel '_DKI_set' num2str(s) '.txt'];
    write_to_file_DKI(folder, filename, fa_rel{s}, md_rel{s}, rd_rel{s}, ad_rel{s}, mk_rel{s}, rk_rel{s}, ak_rel{s});
end

disp('Data saved')

%% Plot: All parameters, one set

w = 7;
h = 5;
fontsize = 20;
params = [{fa_rel} {md_rel} {rd_rel} {ad_rel} {mk_rel} {rk_rel} {ak_rel}];

for s = 1:length(set_nr)
    all_params_setX = [params{1}{s} params{2}{s} params{3}{s} params{4}{s} params{5}{s} params{6}{s} params{7}{s}];    
    figure(s), clf, boxplot(all_params_setX, 'labels', {'FA','MD','RD','AD','MK','RK','AK'})
    set(gca, 'fontsize', fontsize)
    hold on
    plot(linspace(0,105), 0*ones(1,100), '--');
    hold off
    xlabel('Diffusion parameters', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title(['Set ' num2str(s)], 'fontsize', fontsize)
    %pause;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    %saveas(figure(s), fullfile(folder, ['subject' subject '_' voxel '_DKI_set' num2str(s)]), 'pdf');
end
    

%% One parameter, all sets
w = 7;
h = 5;
fontsize = 20;
params = [{fa_rel} {md_rel} {rd_rel} {ad_rel} {mk_rel} {rk_rel} {ak_rel}];
params_string = ['FA'; 'MD'; 'RD'; 'AD'; 'MK'; 'RK'; 'AK'];

for p = 1:length(params)
    one_param_all_sets = [params{p}{1} params{p}{2} params{p}{3} params{p}{4} params{p}{5} params{p}{6} params{p}{7} params{p}{8}];
    figure(p), clf, boxplot(one_param_all_sets, 'labels', {'1','2','3','4','5','6','7','8'})
    set(gca, 'fontsize', fontsize)
    hold on
    plot(linspace(0,105), 0*ones(1,100), '--');
    hold off
    xlabel('Gradient set', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(p,:)], 'fontsize', fontsize)
    %pause;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(p), fullfile(folder, ['subject' subject '_' voxel '_DKI_' params_string(p,:)]), 'pdf');
end

%% Pool subject data together

voxel = 'voxel1';
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/Pool/'];

% DKI
% clear fa_pool_DKI md_pool_DKI rd_pool_DKI ad_pool_DKI mk_pool_DKI rk_pool_DKI ak_pool_DKI 
% for s = 1:8
%     subject002_data_DKI = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/subject002/' voxel '/subject002_' voxel '_DKI_set' num2str(s) '.txt']);
%     subject003_data_DKI= importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/subject003/' voxel '/subject003_' voxel '_DKI_set' num2str(s) '.txt']);
%     subject006_data_DKI = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/subject006/' voxel '/subject006_' voxel '_DKI_set' num2str(s) '.txt']);
% 
%     %subject003_data_DKI(7000) = NaN;
%     %subject006_data_DKI(7000) = NaN;
%     
%     fa_pool_DKI(:,s) = [subject002_data_DKI(1:1000); subject003_data_DKI(1:1000); subject006_data_DKI(1:1000)];
%     md_pool_DKI(:,s) = [subject002_data_DKI(1001:2000); subject003_data_DKI(1001:2000); subject006_data_DKI(1001:2000)];
%     rd_pool_DKI(:,s) = [subject002_data_DKI(2001:3000); subject003_data_DKI(2001:3000); subject006_data_DKI(2001:3000)];
%     ad_pool_DKI(:,s) = [subject002_data_DKI(3001:4000); subject003_data_DKI(3001:4000); subject006_data_DKI(3001:4000)];
%     mk_pool_DKI(:,s) = [subject002_data_DKI(4001:5000); subject003_data_DKI(4001:5000); subject006_data_DKI(4001:5000)];
%     rk_pool_DKI(:,s) = [subject002_data_DKI(5001:6000); subject003_data_DKI(5001:6000); subject006_data_DKI(5001:6000)];
%     ak_pool_DKI(:,s) = [subject002_data_DKI(6001:7000); subject003_data_DKI(6001:end); subject006_data_DKI(6001:end)];
%    
% end
% 
% %SAVE TO FILE
% filename_pool_DKI = [voxel '_DKI.txt'];
% write_to_file_DKI(folder_pool, filename_pool_DKI, fa_pool_DKI, md_pool_DKI, rd_pool_DKI, ad_pool_DKI, mk_pool_DKI, rk_pool_DKI, ak_pool_DKI);


% DTI
clear fa_pool_DTI md_pool_DTI rd_pool_DTI ad_pool_DTI
for s = 7:8
    subject002_data_DTI = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/subject002/' voxel '/subject002_' voxel '_DTI3_set' num2str(s) '.txt']);
    %subject003_data_DTI= importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/subject003/' voxel '/subject003_' voxel '_DTI_set' num2str(s) '.txt']);
    %subject006_data_DTI = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/subject006/' voxel '/subject006_' voxel '_DTI3_set' num2str(s) '.txt']);
 
    %subject002_data_DTI(4000) = NaN;
    %subject006_data_DTI(4000) = NaN;
    
    fa_pool_DTI(:,s) = [subject002_data_DTI(1:1000,:);  subject006_data_DTI(1:1000,:)];
    md_pool_DTI(:,s) = [subject002_data_DTI(1001:2000,:);  subject006_data_DTI(1001:2000,:)];
    rd_pool_DTI(:,s) = [subject002_data_DTI(2001:3000,:); subject006_data_DTI(2001:3000,:)];
    ad_pool_DTI(:,s) = [subject002_data_DTI(3001:4000,:);  subject006_data_DTI(3001:4000,:)];
    
%     fa_pool_DTI(:,s) = [subject002_data_DTI(1:1000,:); subject003_data_DTI(1:1000,:); subject006_data_DTI(1:1000,:)];
%     md_pool_DTI(:,s) = [subject002_data_DTI(1001:2000,:); subject003_data_DTI(1001:2000,:); subject006_data_DTI(1001:2000,:)];
%     rd_pool_DTI(:,s) = [subject002_data_DTI(2001:3000,:); subject003_data_DTI(2001:3000,:); subject006_data_DTI(2001:3000,:)];
%     ad_pool_DTI(:,s) = [subject002_data_DTI(3001:4000,:); subject003_data_DTI(3001:4000,:); subject006_data_DTI(3001:4000,:)];
end

%SAVE TO FILE
filename_pool_DTI = [voxel '_DTI3.txt'];
write_to_file_DTI(folder_pool, filename_pool_DTI, fa_pool_DTI, md_pool_DTI, rd_pool_DTI, ad_pool_DTI);    

disp('ok')