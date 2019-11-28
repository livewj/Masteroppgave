% Plots results of simulation experiment 2 as boxplots

%% Load data: White matter
tissue = 'WM';

% DKI, WLS, NSA=1
data_DKI_WM1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/voxel1_DKI.txt']);
fa_DKI_WM1 = data_DKI_WM1(1:3000,:);
md_DKI_WM1 = data_DKI_WM1(3001:6000,:);
rd_DKI_WM1 = data_DKI_WM1(6001:9000,:);
ad_DKI_WM1 = data_DKI_WM1(9001:12000,:);
mk_DKI_WM1 = data_DKI_WM1(12001:15000,:);
rk_DKI_WM1 = data_DKI_WM1(15001:18000,:);
ak_DKI_WM1 = data_DKI_WM1(18001:21000,:);
% DTI, WLS, b 1000, NSA=1
data_DTI_WM1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/voxel1_DTI.txt']);
fa_DTI_WM1 = data_DTI_WM1(1:3000,:);
md_DTI_WM1 = data_DTI_WM1(3001:6000,:);
rd_DTI_WM1 = data_DTI_WM1(6001:9000,:);
ad_DTI_WM1 = data_DTI_WM1(9001:12000,:);
% DTI, WLS, b 3000, NSA=1
data_DTI_WM3_1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/DTIb3000pool/voxel1_DTI.txt']);
fa_DTI_WM3_1 = data_DTI_WM3_1(1:3000,:);
md_DTI_WM3_1 = data_DTI_WM3_1(3001:6000,:);
rd_DTI_WM3_1 = data_DTI_WM3_1(6001:9000,:);
ad_DTI_WM3_1 = data_DTI_WM3_1(9001:12000,:);


% DKI, WLS, NSA=2
data_DKI_WM2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel1_DKI.txt']);
fa_DKI_WM2 = data_DKI_WM2(1:3000,:);
md_DKI_WM2 = data_DKI_WM2(3001:6000,:);
rd_DKI_WM2 = data_DKI_WM2(6001:9000,:);
ad_DKI_WM2 = data_DKI_WM2(9001:12000,:);
mk_DKI_WM2 = data_DKI_WM2(12001:15000,:);
rk_DKI_WM2 = data_DKI_WM2(15001:18000,:);
ak_DKI_WM2 = data_DKI_WM2(18001:21000,:);
% DTI, WLS, b 1000, NSA=2
data_DTI_WM2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel1_DTI.txt']);
fa_DTI_WM2 = data_DTI_WM2(1:3000,:);
md_DTI_WM2 = data_DTI_WM2(3001:6000,:);
rd_DTI_WM2 = data_DTI_WM2(6001:9000,:);
ad_DTI_WM2 = data_DTI_WM2(9001:12000,:);
% DTI, WLS, b 3000, NSA=2
data_DTI_WM3_2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/DTIb3000pool/voxel1_DTI.txt']);
fa_DTI_WM3_2 = data_DTI_WM3_2(1:3000,:);
md_DTI_WM3_2 = data_DTI_WM3_2(3001:6000,:);
rd_DTI_WM3_2 = data_DTI_WM3_2(6001:9000,:);
ad_DTI_WM3_2 = data_DTI_WM3_2(9001:12000,:);


% DKI and DTI, NSA = 1, 2
fa_pool = [fa_DKI_WM1 fa_DKI_WM2 fa_DTI_WM1 fa_DTI_WM2 fa_DTI_WM3_1 fa_DTI_WM3_2];
md_pool = [md_DKI_WM1 md_DKI_WM2 md_DTI_WM1 md_DTI_WM2 md_DTI_WM3_1 md_DTI_WM3_2];
rd_pool = [rd_DKI_WM1 rd_DKI_WM2 rd_DTI_WM1 rd_DTI_WM2 rd_DTI_WM3_1 rd_DTI_WM3_2];
ad_pool = [ad_DKI_WM1 ad_DKI_WM2 ad_DTI_WM1 ad_DTI_WM2 ad_DTI_WM3_1 ad_DTI_WM3_2];

mk_pool = [mk_DKI_WM1 mk_DKI_WM2];
rk_pool = [rk_DKI_WM1(1:2000,:) rk_DKI_WM2(1:2000,:)];
ak_pool = [ak_DKI_WM1 ak_DKI_WM2];

%% Load data: Gray matter
tissue = 'GM';

% DKI, WLS, NSA=1
data_DKI_GM1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/voxel2_DKI.txt']);
fa_DKI_GM1 = data_DKI_GM1(1:3000,:);
md_DKI_GM1 = data_DKI_GM1(3001:6000,:);
rd_DKI_GM1 = data_DKI_GM1(6001:9000,:);
ad_DKI_GM1 = data_DKI_GM1(9001:12000,:);
mk_DKI_GM1 = data_DKI_GM1(12001:15000,:);
rk_DKI_GM1 = data_DKI_GM1(15001:18000,:);
ak_DKI_GM1 = data_DKI_GM1(18001:21000,:);
% DTI, WLS, b 1000 NSA=1
data_DTI_GM1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/voxel2_DTI.txt']);
fa_DTI_GM1 = data_DTI_GM1(1:3000,:);
md_DTI_GM1 = data_DTI_GM1(3001:6000,:);
rd_DTI_GM1 = data_DTI_GM1(6001:9000,:);
ad_DTI_GM1 = data_DTI_GM1(9001:12000,:);
% DTI, WLS, b 3000 NSA=1
data_DTI_GM3_1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/DTIb3000pool/voxel2_DTI.txt']);
fa_DTI_GM3_1 = data_DTI_GM3_1(1:3000,:);
md_DTI_GM3_1 = data_DTI_GM3_1(3001:6000,:);
rd_DTI_GM3_1 = data_DTI_GM3_1(6001:9000,:);
ad_DTI_GM3_1 = data_DTI_GM3_1(9001:12000,:);


% DKI, WLS, NSA=2
data_DKI_GM2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel2_DKI.txt']);
fa_DKI_GM2 = data_DKI_GM2(1:3000,:);
md_DKI_GM2 = data_DKI_GM2(3001:6000,:);
rd_DKI_GM2 = data_DKI_GM2(6001:9000,:);
ad_DKI_GM2 = data_DKI_GM2(9001:12000,:);
mk_DKI_GM2 = data_DKI_GM2(12001:15000,:);
rk_DKI_GM2 = data_DKI_GM2(15001:18000,:);
ak_DKI_GM2 = data_DKI_GM2(18001:21000,:);
% DTI, WLS, b 1000 NSA=2
data_DTI_GM2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel2_DTI.txt']);
fa_DTI_GM2 = data_DTI_GM2(1:3000,:);
md_DTI_GM2 = data_DTI_GM2(3001:6000,:);
rd_DTI_GM2 = data_DTI_GM2(6001:9000,:);
ad_DTI_GM2 = data_DTI_GM2(9001:12000,:);
% DTI, WLS, b 3000 NSA=2
data_DTI_GM3_2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/DTIb3000pool/voxel2_DTI.txt']);
fa_DTI_GM3_2 = data_DTI_GM3_2(1:3000,:);
md_DTI_GM3_2 = data_DTI_GM3_2(3001:6000,:);
rd_DTI_GM3_2 = data_DTI_GM3_2(6001:9000,:);
ad_DTI_GM3_2 = data_DTI_GM3_2(9001:12000,:);

% DKI and DTI, NSA = 1, 2
fa_pool = [fa_DKI_GM1 fa_DKI_GM2 fa_DTI_GM1 fa_DTI_GM2 fa_DTI_GM3_1 fa_DTI_GM3_2];
md_pool = [md_DKI_GM1 md_DKI_GM2 md_DTI_GM1 md_DTI_GM2 md_DTI_GM3_1 md_DTI_GM3_2];
rd_pool = [rd_DKI_GM1 rd_DKI_GM2 rd_DTI_GM1 rd_DTI_GM2 rd_DTI_GM3_1 rd_DTI_GM3_2];
ad_pool = [ad_DKI_GM1 ad_DKI_GM2 ad_DTI_GM1 ad_DTI_GM2 ad_DTI_GM3_1 ad_DTI_GM3_2];

mk_pool = [mk_DKI_GM1 mk_DKI_GM2];
rk_pool = [rk_DKI_GM1 rk_DKI_GM2];
ak_pool = [ak_DKI_GM1 ak_DKI_GM2];

%% Load data: CSF
tissue = 'CSF';

% DKI, WLS, NSA=1
data_DKI_CSF1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/voxel3_DKI.txt']);
fa_DKI_CSF1 = data_DKI_CSF1(1:3000,:);
md_DKI_CSF1 = data_DKI_CSF1(3001:6000,:);
rd_DKI_CSF1 = data_DKI_CSF1(6001:9000,:);
ad_DKI_CSF1 = data_DKI_CSF1(9001:12000,:);
mk_DKI_CSF1 = data_DKI_CSF1(12001:15000,:);
rk_DKI_CSF1 = data_DKI_CSF1(15001:18000,:);
ak_DKI_CSF1 = data_DKI_CSF1(18001:21000,:);
% DTI, WLS, b 1000 NSA=1
data_DTI_CSF1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/voxel3_DTI.txt']);
fa_DTI_CSF1 = data_DTI_CSF1(1:3000,:);
md_DTI_CSF1 = data_DTI_CSF1(3001:6000,:);
rd_DTI_CSF1 = data_DTI_CSF1(6001:9000,:);
ad_DTI_CSF1 = data_DTI_CSF1(9001:12000,:);
% DTI, WLS, b 3000 NSA=1
data_DTI_CSF3_1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/DTIb3000pool/voxel3_DTI.txt']);
fa_DTI_CSF3_1 = data_DTI_CSF3_1(1:3000,:);
md_DTI_CSF3_1 = data_DTI_CSF3_1(3001:6000,:);
rd_DTI_CSF3_1 = data_DTI_CSF3_1(6001:9000,:);
ad_DTI_CSF3_1 = data_DTI_CSF3_1(9001:12000,:);


% DKI, WLS, NSA=2
data_DKI_CSF2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel3_DKI.txt']);
fa_DKI_CSF2 = data_DKI_CSF2(1:3000,:);
md_DKI_CSF2 = data_DKI_CSF2(3001:6000,:);
rd_DKI_CSF2 = data_DKI_CSF2(6001:9000,:);
ad_DKI_CSF2 = data_DKI_CSF2(9001:12000,:);
mk_DKI_CSF2 = data_DKI_CSF2(12001:15000,:);
rk_DKI_CSF2 = data_DKI_CSF2(15001:18000,:);
ak_DKI_CSF2 = data_DKI_CSF2(18001:21000,:);
% DTI, WLS, b 1000 NSA=2
data_DTI_CSF2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel3_DTI.txt']);
fa_DTI_CSF2 = data_DTI_CSF2(1:3000,:);
md_DTI_CSF2 = data_DTI_CSF2(3001:6000,:);
rd_DTI_CSF2 = data_DTI_CSF2(6001:9000,:);
ad_DTI_CSF2 = data_DTI_CSF2(9001:12000,:);
% DTI, WLS, b 3000 NSA=2
data_DTI_CSF3_2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/DTIb3000pool/voxel3_DTI.txt']);
fa_DTI_CSF3_2 = data_DTI_CSF3_2(1:3000,:);
md_DTI_CSF3_2 = data_DTI_CSF3_2(3001:6000,:);
rd_DTI_CSF3_2 = data_DTI_CSF3_2(6001:9000,:);
ad_DTI_CSF3_2 = data_DTI_CSF3_2(9001:12000,:);



% DKI and DTI, NSA = 1, 2
fa_pool = [fa_DKI_CSF1 fa_DKI_CSF2 fa_DTI_CSF1 fa_DTI_CSF2 fa_DTI_CSF3_1 fa_DTI_CSF3_2];
md_pool = [md_DKI_CSF1 md_DKI_CSF2 md_DTI_CSF1 md_DTI_CSF2 md_DTI_CSF3_1 md_DTI_CSF3_2];
rd_pool = [rd_DKI_CSF1 rd_DKI_CSF2 rd_DTI_CSF1 rd_DTI_CSF2 rd_DTI_CSF3_1 rd_DTI_CSF3_2];
ad_pool = [ad_DKI_CSF1 ad_DKI_CSF2 ad_DTI_CSF1 ad_DTI_CSF2 ad_DTI_CSF3_1 ad_DTI_CSF3_2];

mk_pool = [mk_DKI_CSF1 mk_DKI_CSF2];
rk_pool = [rk_DKI_CSF1 rk_DKI_CSF2];
ak_pool = [ak_DKI_CSF1 ak_DKI_CSF2];



%% Plot DKI, DTI1, DTI3 for NSA = 1, 2 (Groups of 6)

params_pool = [{fa_pool} {md_pool} {rd_pool} {ad_pool}]; % DKI1 DKI2 DTI1 DTI2 DTI3_1 DTI3_2
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/Plots_together_ALL'];
params_string = ['FA'; 'MD'; 'RD'; 'AD'];
group = [1 7 13 19 25 31 37 43 49 55 2 8 14 20 26 32 38 44 50 56 3 9 15 21 27 33 39 45 51 57 4 10 16 22 28 34 40 46 52 58 ...
    5 11 17 23 29 35 41 47 53 59 6 12 18 24 30 36 42 48 54 60];
positions = [1 1.25 1.5 1.75 2 2.25 3 3.25 3.5 3.75 4 4.25 5 5.25 5.5 5.75 6 6.25 7 7.25 7.5 7.75 8 8.25 ...
    9 9.25 9.5 9.75 10 10.25 11 11.25 11.5 11.75 12 12.25 13 13.25 13.5 13.75 14 14.25 15 15.25 15.5 15.75 16 16.25...
    17 17.25 17.5 17.75 18 18.25 19 19.25 19.50 19.75 20 20.25];


color1 = [0, 0.4470, 0.7410];      % DKI color (blue) nsa 1
color2 = [0, 0, 1];                % DKI color (dark blue) nsa 2
color3 = [0.8500, 0.3250, 0.0980]; % DTI color (orange) nsa 1
color4 = [0.9290, 0.6940, 0.1250]; % DTI color (yellow) nsa 2
color5 = [0.55, 0.71, 0];          % DTI3 color (green) nsa 1
color6 = [0, 0.26, 0.15];          % DTI3 color (dark green) nsa 2

color = [{color6},{color5},{color4},{color3},{color2},{color1},{color6},{color5},{color4},{color3},{color2},{color1},{color6},{color5},{color4},{color3},{color2},{color1},{color6},{color5},{color4},{color3},{color2},{color1},{color6},{color5},{color4},{color3},{color2},{color1},...
    {color6},{color5},{color4},{color3},{color2},{color1},{color6},{color5},{color4},{color3},{color2},{color1},{color6},{color5},{color4},{color3},{color2},{color1},{color6},{color5},{color4},{color3},{color2},{color1},{color6},{color5},{color4},{color3},{color2},{color1}];

% x-axis
snr = {'10','15','20','25','30','40','50','60','80','100'};
added_noise = {'1/100','1/80','1/60','1/50','1/40','1/30','1/25','1/20','1/15','1/10'};

for i = 1:length(params_pool)
    x = params_pool{i};
    figure(i), boxplot(x, group, 'positions', positions);
    
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off

    % Group name
    set(gca,'xtick',[mean(positions(1:6)) mean(positions(7:12)) mean(positions(13:18)) mean(positions(19:24)) ...
        mean(positions(25:30)) mean(positions(31:36)) mean(positions(37:42)) mean(positions(43:48)) ...
        mean(positions(49:54)) mean(positions(55:60))])
    set(gca,'xticklabel',added_noise)


    % Change box colors
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color{j},'FaceAlpha',.5);
    end
    
    % Change color of outliers
    outliers = findobj(gcf,'tag', 'Outliers');
    for j = 1:length(outliers)
        outliers(j).MarkerEdgeColor = color{j};
    end
    
    % Legend
    set(gca, 'fontsize', 16)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:6),'DKI, NSA = 1','DKI, NSA = 2',...
        'DTI, b = 1000, NSA = 1','DTI, b = 1000, NSA = 2','DTI, b = 3000, NSA = 1','DTI, b = 3000, NSA = 2');
    legend('location','best');

    xlabel('Added noise [1/SNR]', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 12;
    h = 7;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp2_' tissue '_' params_string(i,:) '_ALL']), 'pdf');

end


