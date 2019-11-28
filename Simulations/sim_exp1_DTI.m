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

% Compute d for each voxel
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

beep on
beep

%% Experiment setup

NSA = 2;    % 1, 2

folder = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/subject' subject];
%folder = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=' num2str(NSA) '/subject' subject '/DTIb3000'];

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
bvals = grad(:,4);
inds1500 = find(bvals < 1500);

% Estimate signal from voxel
clear S_DTI
for v = 1:size(S_voxel,1)
    S = S_voxel(v,:);
    S = reshape(S, [1 1 1 length(S)]);
    [~, d_real, w] = dti_fit(S, grad, [],[],[], 1.5); % b < 1500
    S_DTI(:,v) = w; % Estimated signal
    
    % Estimate parameters
    S = S_DTI(:,v)'; 
    S = reshape(S, [1 1 1 length(S)]);
    [~, d] = dti_fit(S, grad(inds1500,:), [],[],[], 1.5); % b < 1500 
    D = [d(1) d(2) d(3); d(2) d(4) d(5); d(3) d(5) d(6)];
    [MD(v), FA(v), RD(v), AD(v)] = dti_parameters(D); % Ground truth
end

% Compare
[FA_all(inds_voxel) FA; MD_all(inds_voxel) MD; RD_all(inds_voxel) RD; AD_all(inds_voxel) AD]

% Print ground truth
[FA; MD; RD; AD]

% Normalize signal
inds0 = find(bvals < 10);
for i = 1:length(inds0)
    if inds0(i) > length(S_DTI)
        inds0(i) = 0;
    end
end
inds0(find(inds0 == 0)) = [];
S_DTI_norm = S_DTI./nanmean(S_DTI(inds0,:));
S_voxel_norm = S_voxel/mean(S_voxel(:,inds0));

% Compare real and estimated signal
inds3000 = find(bvals < 3000);
fontsize = 22;

figure(3), plot(bvals(inds1500),S_voxel_norm(:,inds1500), 'r*'), hold on, plot(bvals(inds1500), S_DTI_norm,'b*')


%% B) Generate DTI signal analytically based on tensor from real signal

S0 = 1;
G = grad(inds1500,1:3);
b = grad(inds1500,4);

% Generate signal based on real tensor
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
S_DTI_2 = S0*exp(-b.*Y);


% Estimate D from S_DTI
Y = -(1./b).*log(S_DTI_2./S0); 
H_ps = inv(H'*H)*H'; 
d = H_ps*Y;


% Estimate parameters
D = [d(1) d(4) d(5); d(4) d(2) d(6); d(5) d(6) d(3)]*1000;
[MD_2, FA_2, RD_2, AD_2] = dti_parameters(D);

% WLS truth vs analytical ground truth
[FA FA_2; MD MD_2; RD RD_2; AD AD_2]


% Normalize signal
inds0 = find(bvals < 10);
S_DTI_2_norm = S_DTI_2./nanmean(S_DTI_2(inds0,:));


% Plot real signal, DTI fitted signal and analytic DTI signal (one direction)
inds_grad_dir_DTI = [1 7 25];
fontsize = 14;
S0_voxel = mean(S_voxel_norm(inds0));

S_voxel_dir = [S0_voxel S_voxel_norm([7 25])];

figure(4), plot(bvals(inds_grad_dir_DTI), S_voxel_dir, '-*'), ...
    hold on, plot(bvals(inds_grad_dir_DTI), S_DTI_norm(inds_grad_dir_DTI),'g-o'), ...
    hold on, plot(bvals(inds_grad_dir_DTI), S_DTI_2_norm(inds_grad_dir_DTI),'r-o');
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
folder_signal = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/Signal_plot'];
voxel = ['voxel' num2str(vox_nr)];
%saveas(figure(4), fullfile(folder_signal, [subject '_' voxel '_DTI_one_dir_2']), 'pdf');


%% ADD NOISE

S = S_DTI_norm;
%S = S_DTI_2_norm;
nvox=size(S,2);

n = 1000;
SNR = [10 15 20 25 30 40 50 60 80 100];

%Preallocate
FA_data = zeros(n, length(SNR));
MD_data = zeros(n, length(SNR));
RD_data = zeros(n, length(SNR));
AD_data = zeros(n, length(SNR));


%NSA = 2;
%SNR = [20 sqrt(NSA)*20 NSA*20]; 
S_noise_data = zeros(length(S), n, length(SNR), NSA);

% Add noise
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
        S_noise_reshaped = reshape(S_noise, [1 1 1 length(S_noise)]);
        [b0, d_app] = dti_fit(S_noise_reshaped, grad(inds1500,:), [],[],[], 1.5); % b < 1500
        %[b0, d_app] = dti_fit(S_noise_reshaped, grad, [],[],[], 3); % b < 3000

        d = double(d_app);
        D_app = [d(1) d(2) d(3); d(2) d(4) d(5); d(3) d(5) d(6)];

        [MD_n, FA_n, RD_n, AD_n] = dti_parameters(D_app);

        MD_data(i,j) = MD_n;
        FA_data(i,j) = FA_n;
        RD_data(i,j) = RD_n;
        AD_data(i,j) = AD_n;             

    end
end

% Remove extreme outliers: Tukey's method
k= 1.5;
FA_data( FA_data < (prctile(FA_data, 25) - k*iqr(FA_data)) ) = NaN;
FA_data( FA_data > (prctile(FA_data, 75) + k*iqr(FA_data)) ) = NaN;

MD_data( MD_data < (prctile(MD_data, 25) - k*iqr(MD_data)) ) = NaN;
MD_data( MD_data > (prctile(MD_data, 75) + k*iqr(MD_data)) ) = NaN;

RD_data( RD_data < (prctile(RD_data, 25) - k*iqr(RD_data)) ) = NaN;
RD_data( RD_data > (prctile(RD_data, 75) + k*iqr(RD_data)) ) = NaN;

AD_data( AD_data < (prctile(AD_data, 25) - k*iqr(AD_data)) ) = NaN;
AD_data( AD_data > (prctile(AD_data, 75) + k*iqr(AD_data)) ) = NaN;


% Average value from all NSA
% if NSA > 1
%     FA_data = nanmean(FA_data,3);
%     MD_data = nanmean(MD_data,3);
%     RD_data = nanmean(RD_data,3);
%     AD_data = nanmean(AD_data,3);
% end

% Average of all voxels
FA_app = mean(FA_data,4);
MD_app = mean(MD_data,4);
RD_app = mean(RD_data,4);
AD_app = mean(AD_data,4);


beep on
beep
disp('SNR finised')

%% SNR histogram

% all b vals
% NSA = 1
S_noise_40_NSA1 = S_noise_data(:,:,3,1); 
S_noise_sqrt20_NSA1 = S_noise_data(:,:,2,1); 
S_noise_20_NSA1 = S_noise_data(:,:,1,1); 


% NSA = 2
S_noise_NSA2 = mean(S_noise_data,4); % midler over

S_noise_40_NSA2 = S_noise_NSA2(:,:,3);
S_noise_sqrt20_NSA2 = S_noise_NSA2(:,:,2);
S_noise_20_NSA2 = S_noise_NSA2(:,:,1);

% signal at b = 500
find(bvals < 600 & bvals > 400);
% NSA = 1
S_noise_40_NSA1_b500 = S_noise_40_NSA1(12,:);
S_noise_sqrt20_NSA1_b500 = S_noise_sqrt20_NSA1(12,:);
S_noise_20_NSA1_b500 = S_noise_20_NSA1(12,:);

% NSA = 2
S_noise_40_NSA2_b500 = S_noise_40_NSA2(12,:);
S_noise_sqrt20_NSA2_b500 = S_noise_sqrt20_NSA2(12,:);
S_noise_20_NSA2_b500 = S_noise_20_NSA2(12,:);

% Histogrammet viser at SNR går som sqrt(NSA)
% SNR = 20 ved NSA = 2 er ca det samme som SNR = sqrt(2)*20 ved NSA = 1
histogram(S_noise_20_NSA2_b500), hold on, histogram(S_noise_sqrt20_NSA1_b500), hold on, histogram(S_noise_40_NSA1_b500)
legend(['SNR = 20, nsa = 2'], ['SNR = sqrt(2)*20, nsa = 1'], ['SNR = 2*20, nsa = 1']);



%% Compute relative error [%]

FA_rel = (FA_app - FA)./FA*100;
MD_rel = (MD_app - MD)./MD*100;
RD_rel = (RD_app - RD)./RD*100;
AD_rel = (AD_app - AD)./AD*100;

% FA_rel = (FA_app - FA_2)./FA_2*100;
% MD_rel = (MD_app - MD_2)./MD_2*100;
% RD_rel = (RD_app - RD_2)./RD_2*100;
% AD_rel = (AD_app - AD_2)./AD_2*100;

disp('Relative error computed')

%% SAVE RELATIVE ERROR TO FILE

voxel = ['voxel' num2str(vox_nr)];
filename = ['/subject' subject '_' voxel '_DTI2_rel_nsa.txt'];
write_to_file_DTI(folder, filename, FA_rel, MD_rel, RD_rel, AD_rel);

disp('Data saved')


%% Plot parameteres

% ------ Boxplot --------

w = 7;
h = 5;
fontsize = 20;

params = [{FA_rel} {MD_rel} {RD_rel} {AD_rel}];
params_string = ['FA'; 'MD'; 'RD'; 'AD'];

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
    %saveas(figure(i), fullfile(folder, ['subject' subject '_' voxel '_DTI2_' params_string(i,:)]), 'pdf');
end

