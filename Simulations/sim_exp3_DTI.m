% Simulation experiment 3:
% Evaluating the influence of the number of gradient directions and
% b-values on the parameter estimations
% SNR = constant

%% Load data 

subject = '002'; % 002, 003 or 006

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


%% Create gradient sets
%bvals = sortedgrad(:,4);
bvals = grad(:,4);

inds0 = find(bvals<10);
inds500 = find(400<bvals & bvals<600);
inds1000 = find(800<bvals & bvals<1200);
inds2000 = find(1900<bvals & bvals<2100);
inds3000 = find(2900<bvals & bvals<3100);

% DTI gradient sets
% S_DTI
set1 = [inds0(1); inds1000(1:2:end)]; % 1 + 15
set2 = [inds0(1); inds1000]; % 1 + 30
set3 = [inds0; inds1000];    % 6 + 30
set4 = [inds0(1); inds500(1:2:end); inds1000(1:2:end)]; % 1 + 6 + 15
set5 = [inds0(1); inds500; inds1000]; % 1 + 12 + 30
set6 = [inds0; inds500; inds1000]; % 6 + 12 + 30

% S_DTI3
set7 = [inds0(1); inds500; inds1000; inds3000]; % 1 + 12 + 30 + 50
set8 = [inds0; inds500; inds1000; inds3000];    % 6 + 12 + 30 + 50


% ENSURE ISOTROPICALLY DISTRIBUTED GRADIENT DIRECTIONS
% Plot gradient directions
% inds = inds1000(1:5:end)
% figure, hold on, 
% for i=1:length(inds)
%     plot3([-grad_roi(inds(i),1) grad_roi(inds(i),1)], [-grad_roi(inds(i),2) grad_roi(inds(i),2)], [-grad_roi(inds(i),3) grad_roi(inds(i),3)]);
%     %pause;
%     %axis equal;
%      set(gca,'FontSize',14)
%      xlabel('x','FontSize',14)
%      ylabel('y','FontSize',14)
%      zlabel('z','FontSize',14)
%      title('Gradient directions','FontSize',14)
% end


disp('ok')
%% Noise-free parameter estimation based on real signal

% Compute d for each voxel (Full set)
FA_all = NaN(106,106); MD_all = NaN(106,106); RD_all = NaN(106,106); AD_all = NaN(106,106);
for v = 1:length(S_brain)
    if mask_whole(v) == 0
        continue
    else
        S = S_brain(v,:);
        S_resh = reshape(S, [1 1 1 length(S)]);
        [~, d] = dti_fit(S_resh, grad(:,:), [], [], [], 1.5);

        D = [d(1) d(2) d(3); d(2) d(4) d(5); d(3) d(5) d(6)];
        [MD_all(v), FA_all(v), RD_all(v), AD_all(v)] = dti_parameters(D);
    end
end

FA_resh = reshape(FA_all, [nx, ny]);
imagesc(FA_resh), colormap gray, axis equal


%% Experminent setup

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

%% A) Generate signal and establish ground truth from selected voxel 
inds1500 = find(bvals < 1500);

% Estimate signal
clear S_DTI
for v = 1:size(S_voxel,1)
    S = S_voxel(v,:);
    S = reshape(S, [1 1 1 length(S)]);
    [~, d_real, w] = dti_fit(S, grad, [],[],[], 1.5);
    %[~, ~, w] = dti_fit(S, sortedgrad, [],[],[],1.5); % subject 003
    S_DTI(:,v) = w; % Estimated signal
    
    % Estimate parameters
    S = S_DTI(:,v)'; 
    S = reshape(S, [1 1 1 length(S)]);
    [~, d] = dti_fit(S, grad(inds1500,:), [],[],[], 1.5);
   % [~, d] = dti_fit(S, sortedgrad(inds1500,:), [],[],[], 1.5); % subject 003
    D = [d(1) d(2) d(3); d(2) d(4) d(5); d(3) d(5) d(6)];
    [MD(v), FA(v), RD(v), AD(v)] = dti_parameters(D);
end

% Compare estimated (ground truth) parameters with real-signal parameters
[FA_all(inds_voxel) FA; MD_all(inds_voxel) MD; RD_all(inds_voxel) RD; AD_all(inds_voxel) AD];

% Print ground truth
[FA; MD; RD; AD];


% Generate DTI signal with b 3000
clear S_DTI_2
for v = 1:size(S_voxel,1)
    S = S_voxel(v,:);
    S = reshape(S, [1 1 1 length(S)]);
    [~, ~, w] = dti_fit(S, grad, [],[],[], 3);
    %[~, ~, w] = dti_fit(S, sortedgrad, [],[],[],3); % subject 003
    S_DTI_2(:,v) = w; % Estimated signal
end



% Normalize signal
inds0 = find(bvals < 10);
for i = 1:length(inds0)
    if inds0(i) > length(S_DTI)
        inds0(i) = 0;
    end
end
inds0(find(inds0 == 0)) = [];

S_DTI_norm = S_DTI./nanmean(S_DTI(inds0,:));
S_DTI_2_norm = S_DTI_2./nanmean(S_DTI_2(inds0,:));
S_voxel_norm = S_voxel/mean(S_voxel(:,inds0));

% Compare real and estimated signal
inds3000 = find(bvals < 3000);
inds1500 = find(bvals < 1500);
plot(bvals(inds3000),S_voxel_norm(:,inds3000), 'r*'), hold on, plot(bvals(inds1500), S_DTI_norm,'b*'), ...
    hold on, plot(bvals(inds3000), S_DTI_2_norm, 'g*')

% figure(4), plot(bvals(inds_grad_dir_DTI3), S_voxel_dir, '-*'), ...
%     hold on, plot(bvals(inds_grad_dir_DTI), S_DTI_norm(inds_grad_dir_DTI),'g-o'), ...
%     hold on, plot(bvals(inds_grad_dir_DTI3), S_DTI3_norm(inds_grad_dir_DTI3),'r-o');
% set(gca,'fontsize', fontsize)
% xlabel('b-value [s/mm^2]','fontsize', fontsize)
% ylabel('Scaled signal intensity, S/S_0','fontsize', fontsize)
% title(['Actual vs fitted signal, white matter'],'fontsize', fontsize)
% legend(['Actual signal'], ['DTI, b < 1000'], ['DTI, b < 3000'],'fontsize', fontsize)
% w = 7;
% h = 5;
% set(gcf, 'PaperSize', [w h]);
% set(gcf, 'PaperPosition', [0 0 w h]);
% folder_signal = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/Signal_plot'];
% voxel = ['voxel' num2str(vox_nr)];
%saveas(figure(4), fullfile(folder_signal, [subject '_' voxel '_DTIb1000_DTIb3000']), 'pdf');
 
%% B) Generate DTI signal analytically based on tensor from real signal

S0 = 1;
G = grad(inds3000,1:3);
b = grad(inds3000,4);

% Generate signal based on tensor from experimental DW-MRI data
%d_real = [1.6053, 0.1482, -0.0631, 0.1654, -0.0197, 0.1202]*1e-3; % White matter
%d_real = [0.6940, 0.0090, -0.0606, 0.6893, 0.0854, 0.7667]*1e-3; % GM
%d_real = [2.1055, 0.0064, 0.0443, 2.5200, -0.0087, 2.3436]*1e-3; % CSF

% Change order
D = [d_real(1) d_real(2) d_real(3); ...
    d_real(2) d_real(4) d_real(5); ...
    d_real(3) d_real(5) d_real(6)];
d_real = [D(1,1) D(2,2) D(3,3) D(1,2) D(1,3) D(2,3)]*1e-3;


for i = 1:length(G)
    H(i,1) = G(i,1)^2;
    H(i,2) = G(i,2)^2;
    H(i,3) = G(i,3)^2;
    H(i,4) = 2*G(i,1)*G(i,2);
    H(i,5) = 2*G(i,1)*G(i,3);
    H(i,6) = 2*G(i,2)*G(i,3);
end

Y = H*d_real';

% Calculate signal
S_DTI3_2 = S0*exp(-b.*Y);


% Estimate D from S_DTI
Y = -(1./b).*log(S_DTI3_2./S0); 
H_ps = inv(H'*H)*H'; 
d = H_ps*Y;


% Estimate parameters
D = [d(1) d(4) d(5); d(4) d(2) d(6); d(5) d(6) d(3)]*1000;
[MD_2, FA_2, RD_2, AD_2] = dti_parameters(D);

% WLS truth vs analytical ground truth
[FA FA_2; MD MD_2; RD RD_2; AD AD_2]


% Normalize signal
inds0 = find(bvals < 10);
S_DTI3_2_norm = S_DTI3_2./nanmean(S_DTI3_2(inds0,:));


% Plot real signal, DTI fitted signal and analytic DTI signal (one direction)
inds_grad_dir_DTI = [1 7 25];
fontsize = 14;
S0_voxel = mean(S_voxel_norm(inds0));

S_voxel_dir = [S0_voxel S_voxel_norm([7 25])];

figure(4), plot(bvals(inds_grad_dir_DTI), S_voxel_dir, '-*'), ...
    hold on, plot(bvals(inds_grad_dir_DTI), S_DTI_norm(inds_grad_dir_DTI),'g-o'), ...
    hold on, plot(bvals(inds_grad_dir_DTI), S_DTI3_2_norm(inds_grad_dir_DTI),'r-o');
set(gca,'fontsize', fontsize)
xlabel('b-value [s/mm^2]','fontsize', fontsize)
ylabel('Scaled signal intensity, S/S_0','fontsize', fontsize)
%title(['Actual vs fitted signal, white matter'],'fontsize', fontsize)
%title(['Actual vs fitted signal, gray matter'],'fontsize', fontsize)
%title(['Actual vs fitted signal, CSF'],'fontsize', fontsize)
legend(['Actual signal'], ['Fitted signal DTI, WLS'], ['Fitted signal DTI, analytic'],'fontsize', fontsize)
w = 7;
h = 5;
set(gcf, 'PaperSize', [w h]);
set(gcf, 'PaperPosition', [0 0 w h]);

%% ADD NOISE

%nvox=size(S_DTI_norm,2);
%S = S_DTI_2_norm;  % analytic DTI signal b 1000
S = S_DTI3_2_norm; % analytic DTI signal b 3000
nvox = 1;

n = 1000;
SNR = 50;

%set_nr = [{set1} {set2} {set3} {set4} {set5} {set6} {set7} {set8}];
%set_nr = [{set1} {set2} {set3} {set4} {set5} {set6}];
set_nr = [{set7} {set8}];

%Preallocate
FA_data = zeros(n, length(SNR), NSA, nvox);
MD_data = zeros(n, length(SNR), NSA, nvox);
RD_data = zeros(n, length(SNR), NSA, nvox);
AD_data = zeros(n, length(SNR), NSA, nvox);

% MAIN LOOP
for s = 1:length(set_nr)
    inds = set_nr{s};
    %S = S_DTI3_norm(inds);
    %S = S_DTI_2_norm(inds);
    S = S_DTI3_2_norm(inds);
    for v = 1:nvox
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
                    S_noise_reshaped = reshape(S_noise, [1 1 1 length(S_noise)]);
                    [~, d_app] = dti_fit(S_noise_reshaped, grad(inds,:));
                    %[~, d_app] = dti_fit(S_noise_reshaped, sortedgrad(inds,:));
                
                    d = double(d_app);
                    D_app = [d(1) d(2) d(3); d(2) d(4) d(5); d(3) d(5) d(6)];

                    [MD_n, FA_n, RD_n, AD_n] = dti_parameters(D_app);
                    
                    FA_data(i,j,k,v) = FA_n;
                    MD_data(i,j,k,v) = MD_n;
                    RD_data(i,j,k,v) = RD_n;
                    AD_data(i,j,k,v) = AD_n;   
                    
                    FA_set{s} = FA_data;
                    MD_set{s} = MD_data;
                    RD_set{s} = RD_data;
                    AD_set{s} = AD_data;

                end
            end
        end
    end
end


% FOR ALL SETS:
% Remove extreme outliers: Tukey's method
k = 1.5;
for s = 1:length(set_nr)
    FA_set{s}( FA_set{s} < (prctile(FA_set{s}, 25) - k*iqr(FA_set{s})) ) = NaN;
    FA_set{s}( FA_set{s} > (prctile(FA_set{s}, 75) + k*iqr(FA_set{s})) ) = NaN;

    MD_set{s}( MD_set{s} < (prctile(MD_set{s}, 25) - k*iqr(MD_set{s})) ) = NaN;
    MD_set{s}( MD_set{s} > (prctile(MD_set{s}, 75) + k*iqr(MD_set{s})) ) = NaN;
    
    RD_set{s}( RD_set{s} < (prctile(RD_set{s}, 25) - k*iqr(RD_set{s})) ) = NaN;
    RD_set{s}( RD_set{s} > (prctile(RD_set{s}, 75) + k*iqr(RD_set{s})) ) = NaN;
    
    AD_set{s}( AD_set{s} < (prctile(AD_set{s}, 25) - k*iqr(AD_set{s})) ) = NaN;
    AD_set{s}( AD_set{s} > (prctile(AD_set{s}, 75) + k*iqr(AD_set{s})) ) = NaN;
end


% Average value from all NSA
clear FA_nsa; clear MD_nsa; clear RD_nsa; clear AD_nsa;
for s = 1:length(set_nr)
    FA_nsa{s} = nanmean(FA_set{s},3);
    MD_nsa{s} = nanmean(MD_set{s},3);
    RD_nsa{s} = nanmean(RD_set{s},3);
    AD_nsa{s} = nanmean(AD_set{s},3);
end

% Average of all voxels
clear FA_app, clear MD_app, clear RD_app, clear AD_app;
for s = 1:length(set_nr)
    FA_app{s} = mean(FA_nsa{s},4);
    MD_app{s} = mean(MD_nsa{s},4);
    RD_app{s} = mean(RD_nsa{s},4);
    AD_app{s} = mean(AD_nsa{s},4); 
end

disp('SNR finished')
%% Compute relative error [%]

clear MD_rel, clear FA_rel, clear RD_rel, clear AD_rel
% for s = 1:length(set_nr)
%     FA_rel{s} = (FA_app{s} - FA)./FA*100;
%     MD_rel{s} = (MD_app{s} - MD)./MD*100;
%     RD_rel{s} = (RD_app{s} - RD)./RD*100;
%     AD_rel{s} = (AD_app{s} - AD)./AD*100;
% end

for s = 1:length(set_nr)
    FA_rel{s} = (FA_app{s} - FA_2)./FA_2*100;
    MD_rel{s} = (MD_app{s} - MD_2)./MD_2*100;
    RD_rel{s} = (RD_app{s} - RD_2)./RD_2*100;
    AD_rel{s} = (AD_app{s} - AD_2)./AD_2*100;
end

disp('Relative error computed')


%% SAVE RELATIVE ERROR TO FILE

voxel = ['voxel' num2str(vox_nr)];
folder = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_3/subject' subject '/' voxel ];
%set_str = ['1', '2', '3', '4', '5','6']; % b < 1500
set_str = ['7','8']; % b < 3000
%set_str = ['1', '2', '3', '4', '5','6','7','8']; % b < 1500

for s = 1:length(set_nr)
    filename = ['/subject' subject '_' voxel '_DTI3_2_set' set_str(s) '.txt'];
    write_to_file_DTI(folder, filename, FA_rel{s}, MD_rel{s}, RD_rel{s}, AD_rel{s});
end

disp('Data saved')

%% All parameters, one set

w = 7;
h = 5;

params = [{FA_rel} {MD_rel} {RD_rel} {AD_rel}];

    
for s = 1:length(set_nr)
    all_params_setX = [params{1}{s} params{2}{s} params{3}{s} params{4}{s}];    
    figure(s), clf, boxplot(all_params_setX, 'labels', {'FA','MD','RD','AD'})
    set(gca, 'fontsize', 14)
    hold on
    plot(linspace(0,105), 0*ones(1,100), '--');
    hold off
    xlabel('Diffusion parameters', 'fontsize', 14)
    ylabel('Relative error [%]', 'fontsize', 14);
    title(['Set ' set_str(s)], 'fontsize', 14)
    %pause;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    %saveas(figure(s), fullfile(folder, ['subject' subject '_' voxel '_DTI_set' set_str(s)]), 'pdf');
end
    

%% One parameter, all sets

w = 7;
h = 5;

params = [{FA_rel} {MD_rel} {RD_rel} {AD_rel}];
params_string = ['FA'; 'MD'; 'RD'; 'AD'];

for p = 1:length(params)
    one_param_all_sets = [params{p}{1} params{p}{2}];
    %one_param_all_sets = [params{p}{1} params{p}{2} params{p}{3} params{p}{4} params{p}{5} params{p}{6}];
    %figure(p), clf, boxplot(one_param_all_sets, 'labels', {'1','2','3','4','5','6','7','8'})
    %figure(p), clf, boxplot(one_param_all_sets, 'labels', {'1','2','3','4','5','6'})
    figure(p), clf, boxplot(one_param_all_sets, 'labels', {'7','8'})
    set(gca, 'fontsize', 14)
    hold on
    plot(linspace(0,105), 0*ones(1,100), '--');
    hold off
    xlabel('Gradient set', 'fontsize', 14)
    ylabel('Relative error [%]', 'fontsize', 14);
    title([params_string(p,:)], 'fontsize', 14)
    %pause;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    %saveas(figure(p), fullfile(folder, ['subject' subject '_' voxel '_DTI_' params_string(p,:)]), 'pdf');
end


