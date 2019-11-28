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
% DTI, WLS, NSA=1
data_DTI_WM1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/voxel1_DTI.txt']);
fa_DTI_WM1 = data_DTI_WM1(1:3000,:);
md_DTI_WM1 = data_DTI_WM1(3001:6000,:);
rd_DTI_WM1 = data_DTI_WM1(6001:9000,:);
ad_DTI_WM1 = data_DTI_WM1(9001:12000,:);

% DKI, WLS, NSA=2
data_DKI_WM2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel1_DKI_22.txt']);
fa_DKI_WM2 = data_DKI_WM2(1:3000,:);
md_DKI_WM2 = data_DKI_WM2(3001:6000,:);
rd_DKI_WM2 = data_DKI_WM2(6001:9000,:);
ad_DKI_WM2 = data_DKI_WM2(9001:12000,:);
mk_DKI_WM2 = data_DKI_WM2(12001:15000,:);
rk_DKI_WM2 = data_DKI_WM2(15001:18000,:);
ak_DKI_WM2 = data_DKI_WM2(18001:21000,:);
% DTI, WLS, NSA=2
data_DTI_WM2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel1_DTI_22.txt']);
fa_DTI_WM2 = data_DTI_WM2(1:3000,:);
md_DTI_WM2 = data_DTI_WM2(3001:6000,:);
rd_DTI_WM2 = data_DTI_WM2(6001:9000,:);
ad_DTI_WM2 = data_DTI_WM2(9001:12000,:);

% DKI and DTI, NSA = 1, 2
fa_pool = [fa_DKI_WM1 fa_DKI_WM2 fa_DTI_WM1 fa_DTI_WM2];
md_pool = [md_DKI_WM1 md_DKI_WM2 md_DTI_WM1 md_DTI_WM2];
rd_pool = [rd_DKI_WM1 rd_DKI_WM2 rd_DTI_WM1 rd_DTI_WM2];
ad_pool = [ad_DKI_WM1 ad_DKI_WM2 ad_DTI_WM1 ad_DTI_WM2];

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
% DTI, WLS, NSA=1
data_DTI_GM1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/voxel2_DTI.txt']);
fa_DTI_GM1 = data_DTI_GM1(1:3000,:);
md_DTI_GM1 = data_DTI_GM1(3001:6000,:);
rd_DTI_GM1 = data_DTI_GM1(6001:9000,:);
ad_DTI_GM1 = data_DTI_GM1(9001:12000,:);

% DKI, WLS, NSA=2
data_DKI_GM2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel2_DKI_22.txt']);
fa_DKI_GM2 = data_DKI_GM2(1:3000,:);
md_DKI_GM2 = data_DKI_GM2(3001:6000,:);
rd_DKI_GM2 = data_DKI_GM2(6001:9000,:);
ad_DKI_GM2 = data_DKI_GM2(9001:12000,:);
mk_DKI_GM2 = data_DKI_GM2(12001:15000,:);
rk_DKI_GM2 = data_DKI_GM2(15001:18000,:);
ak_DKI_GM2 = data_DKI_GM2(18001:21000,:);
% DTI, WLS, NSA=2
data_DTI_GM2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel2_DTI_22.txt']);
fa_DTI_GM2 = data_DTI_GM2(1:3000,:);
md_DTI_GM2 = data_DTI_GM2(3001:6000,:);
rd_DTI_GM2 = data_DTI_GM2(6001:9000,:);
ad_DTI_GM2 = data_DTI_GM2(9001:12000,:);

% DKI and DTI, NSA = 1, 2
fa_pool = [fa_DKI_GM1 fa_DKI_GM2 fa_DTI_GM1 fa_DTI_GM2];
md_pool = [md_DKI_GM1 md_DKI_GM2 md_DTI_GM1 md_DTI_GM2];
rd_pool = [rd_DKI_GM1 rd_DKI_GM2 rd_DTI_GM1 rd_DTI_GM2];
ad_pool = [ad_DKI_GM1 ad_DKI_GM2 ad_DTI_GM1 ad_DTI_GM2];

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
% DTI, WLS, NSA=1
data_DTI_CSF1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/Pool/voxel3_DTI.txt']);
fa_DTI_CSF1 = data_DTI_CSF1(1:3000,:);
md_DTI_CSF1 = data_DTI_CSF1(3001:6000,:);
rd_DTI_CSF1 = data_DTI_CSF1(6001:9000,:);
ad_DTI_CSF1 = data_DTI_CSF1(9001:12000,:);

% DKI, WLS, NSA=2
data_DKI_CSF2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel3_DKI_22.txt']);
fa_DKI_CSF2 = data_DKI_CSF2(1:3000,:);
md_DKI_CSF2 = data_DKI_CSF2(3001:6000,:);
rd_DKI_CSF2 = data_DKI_CSF2(6001:9000,:);
ad_DKI_CSF2 = data_DKI_CSF2(9001:12000,:);
mk_DKI_CSF2 = data_DKI_CSF2(12001:15000,:);
rk_DKI_CSF2 = data_DKI_CSF2(15001:18000,:);
ak_DKI_CSF2 = data_DKI_CSF2(18001:21000,:);
% DTI, WLS, NSA=2
data_DTI_CSF2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/voxel3_DTI_22.txt']);
fa_DTI_CSF2 = data_DTI_CSF2(1:3000,:);
md_DTI_CSF2 = data_DTI_CSF2(3001:6000,:);
rd_DTI_CSF2 = data_DTI_CSF2(6001:9000,:);
ad_DTI_CSF2 = data_DTI_CSF2(9001:12000,:);


% DKI and DTI, NSA = 1, 2
fa_pool = [fa_DKI_CSF1 fa_DKI_CSF2 fa_DTI_CSF1 fa_DTI_CSF2];
md_pool = [md_DKI_CSF1 md_DKI_CSF2 md_DTI_CSF1 md_DTI_CSF2];
rd_pool = [rd_DKI_CSF1 rd_DKI_CSF2 rd_DTI_CSF1 rd_DTI_CSF2];
ad_pool = [ad_DKI_CSF1 ad_DKI_CSF2 ad_DTI_CSF1 ad_DTI_CSF2];

mk_pool = [mk_DKI_CSF1 mk_DKI_CSF2];
rk_pool = [rk_DKI_CSF1 rk_DKI_CSF2];
ak_pool = [ak_DKI_CSF1 ak_DKI_CSF2];


%% Plot DKI and DTI together for NSA = 2

fa_pool2 = [fa_pool(:,11:20) fa_pool(:,31:40)]; % DKI2, DTI2
md_pool2 = [md_pool(:,11:20) md_pool(:,31:40)];
rd_pool2 = [rd_pool(:,11:20) rd_pool(:,31:40)];
ad_pool2 = [ad_pool(:,11:20) ad_pool(:,31:40)];

params_pool = [{fa_pool2} {md_pool2} {rd_pool2} {ad_pool2}];
params_string = ['FA'; 'MD'; 'RD'; 'AD'];
group = [1 3 5 7 9 11 13 15 17 19 2 4 6 8 10 12 14 16 18 20];
positions = [1 1.25 2 2.25 3 3.25 4 4.25 5 5.25 6 6.25 7 7.25 8 8.25 9 9.25 10 10.25];
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/Plots_together'];
fontsize = 20;

color3 = [0.9290, 0.6940, 0.1250]; % DTI color (yellow) for NSA = 2
color4 = [0, 0, 1]; % DKI color (dark blue) for NSA = 2

color = [{color3},{color4},{color3},{color4},{color3},{color4},{color3},{color4},{color3},{color4},{color3},{color4},{color3},{color4},{color3},{color4},{color3},{color4},{color3},{color4}];


for i = 1:length(params_pool)
    x = params_pool{i};
    figure(i), boxplot(x, group, 'positions', positions);
    
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off

    % Group name
    set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8)) ...
        mean(positions(9:10)) mean(positions(11:12)) mean(positions(13:14)) mean(positions(15:16)) ...
        mean(positions(17:18)) mean(positions(19:20))])
    set(gca,'xticklabel',{'10','15','20','25','30','40','50','60','80','100'})


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
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:2), 'DKI, NSA = 2', 'DTI, NSA = 2');


    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 7;
    h = 5;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp2_' tissue '_' params_string(i,:) 'DKIDTI2']), 'pdf');

end

%% Plot only the DKI parameters for NSA = 2

params_pool = [{mk_pool(:,11:20)} {rk_pool(:,11:20)} {ak_pool(:,11:20)}];
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/Plots_together'];
params_string = ['MK'; 'RK'; 'AK'];

for i = 1:length(params_pool)
    figure(i), clf, boxplot(params_pool{i}, SNR)
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off
    
    % Change box colors
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color4,'FaceAlpha',.5);
    end
    
    % Change color of outliers
    outliers = findobj(gcf,'tag', 'Outliers');
    for j = 1:length(outliers)
        outliers(j).MarkerEdgeColor = color4;
    end
    
    % Legend
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:1), 'DKI, NSA = 2');
    
    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 7;
    h = 5;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp2_' tissue '_' params_string(i,:) 'DKI2']), 'pdf');
    
end

%% Plot only the DTI parameters for NSA = 1,2

fa_pool2 = [fa_pool(:,21:30) fa_pool(:,31:40)]; % DTI1 DTI2
md_pool2 = [md_pool(:,21:30) md_pool(:,31:40)];
rd_pool2 = [rd_pool(:,21:30) rd_pool(:,31:40)];
ad_pool2 = [ad_pool(:,21:30) ad_pool(:,31:40)];

params_pool = [{fa_pool2} {md_pool2} {rd_pool2} {ad_pool2}];
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/Plots_together'];
params_string = ['FA'; 'MD'; 'RD'; 'AD'];
group = [1 3 5 7 9 11 13 15 17 19 2 4 6 8 10 12 14 16 18 20];
positions = [1 1.25 2 2.25 3 3.25 4 4.25 5 5.25 6 6.25 7 7.25 8 8.25 9 9.25 10 10.25];

color1 = [0.8500, 0.3250, 0.0980]; % DTI color (orange) for NSA = 1
color3 = [0.9290, 0.6940, 0.1250]; % DTI color (yellow) for NSA = 2
color = [{color3},{color1},{color3},{color1},{color3},{color1},{color3},{color1},{color3},{color1},{color3},{color1},{color3},{color1},{color3},{color1},{color3},{color1},{color3},{color1}];
    
for i = 1:length(params_pool)
    x = params_pool{i};
    figure(i), boxplot(x, group, 'positions', positions);
    
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off

    % Group name
    set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8)) ...
        mean(positions(9:10)) mean(positions(11:12)) mean(positions(13:14)) mean(positions(15:16)) ...
        mean(positions(17:18)) mean(positions(19:20))])
    set(gca,'xticklabel',{'10','15','20','25','30','40','50','60','80','100'})


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
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:2), 'DTI, NSA = 1', 'DTI, NSA = 2');


    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 7;
    h = 5;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp2_' tissue '_' params_string(i,:) 'DTI12']), 'pdf');

end



%% Plot only the DKI parameters for NSA = 1,2

params_pool = [{mk_pool} {rk_pool} {ak_pool}];
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/Plots_together'];
params_string = ['MK'; 'RK'; 'AK'];
group = [1 3 5 7 9 11 13 15 17 19 2 4 6 8 10 12 14 16 18 20];
positions = [1 1.25 2 2.25 3 3.25 4 4.25 5 5.25 6 6.25 7 7.25 8 8.25 9 9.25 10 10.25];

color2 = [0, 0.4470, 0.7410]; % DKI color (blue) for NSA = 1
color4 = [0, 0, 1];  % DKI color (dark blue) for NSA = 2
color = [{color4},{color2},{color4},{color2},{color4},{color2},{color4},{color2},{color4},{color2},{color4},{color2},{color4},{color2},{color4},{color2},{color4},{color2},{color4},{color2}];

for i = 1:length(params_pool)
    x = params_pool{i};
    figure(i), boxplot(x, group, 'positions', positions);
    
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off

    % Group name
    set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8)) ...
        mean(positions(9:10)) mean(positions(11:12)) mean(positions(13:14)) mean(positions(15:16)) ...
        mean(positions(17:18)) mean(positions(19:20))])
    set(gca,'xticklabel',{'10','15','20','25','30','40','50','60','80','100'})


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
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:2), 'DKI, NSA = 1', 'DKI, NSA = 2');


    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 7;
    h = 5;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp2_' tissue '_' params_string(i,:) 'DKI12']), 'pdf');

end

%% Plot DTI and DKI for NSA = 1, 2 (Groups of 4)

%fa_NSA1 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=1/subject002/subject002_voxel1_DKI_rel.txt']);
%fa_NSA2 = importdata(['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/subject002/subject002_voxel1_DKI_rel_sqrt.txt']);

%fa_pool = [fa_NSA1 fa_NSA2 fa_NSA1 fa_NSA2];


params_pool = [{fa_pool} {md_pool} {rd_pool} {ad_pool}]; % DKI1 DKI2 DTI1 DTI2
%params_pool = {fa_pool}; % DKI nsa1, DKIsqrt nsa 2
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/Plots_together'];
%params_string = ['FA'];
params_string = ['FA'; 'MD'; 'RD'; 'AD'];
group = [1 5 9 13 17 21 25 29 33 37 2 6 10 14 18 22 26 30 34 38 3 7 11 15 19 23 27 31 35 39 4 8 12 16 20 24 28 32 36 40];
positions = [1 1.25 1.5 1.75 2.5 2.75 3 3.25 4 4.25 4.5 4.75 5.5 5.75 6 6.25 7 7.25 7.5 7.75 8.5 8.75 9 9.25 ...
    10 10.25 10.5 10.75 11.5 11.75 12 12.25 13 13.25 13.5 13.75 14.5 14.75 15 15.25];
    


color1 = [0.8500, 0.3250, 0.0980]; % DTI color (orange) for NSA = 1
color3 = [0.9290, 0.6940, 0.1250]; % DTI color (yellow) for NSA = 2
color2 = [0, 0.4470, 0.7410];      % DKI color (blue) for NSA = 1
color4 = [0, 0, 1];                % DKI color (dark blue) for NSA = 2

color = [{color3},{color1},{color4},{color2},{color3},{color1},{color4},{color2},{color3},{color1},{color4},{color2},{color3},{color1},{color4},{color2},{color3},{color1},{color4},{color2},...
    {color3},{color1},{color4},{color2},{color3},{color1},{color4},{color2},{color3},{color1},{color4},{color2},{color3},{color1},{color4},{color2},{color3},{color1},{color4},{color2}];
 

for i = 1:length(params_pool)
    x = params_pool{i};
    figure(i), boxplot(x, group, 'positions', positions);
    
    hold on
    plot(linspace(0,105), 0*ones(1,100), '-');
    hold off

    % Group name
    set(gca,'xtick',[mean(positions(1:4)) mean(positions(5:8)) mean(positions(9:12)) mean(positions(13:16)) ...
        mean(positions(17:20)) mean(positions(21:24)) mean(positions(25:28)) mean(positions(28:32)) ...
        mean(positions(33:36)) mean(positions(37:40))])
    set(gca,'xticklabel',{'10','15','20','25','30','40','50','60','80','100'})


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
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:4),'DKI, NSA = 1','DKI, NSA = 2','DTI, NSA = 1','DTI, NSA = 2');
    %legend('location','southeast');

    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 12;
    h = 6;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp2_' tissue '_' params_string(i,:) 'DTIDKI12_22']), 'pdf');

end

%% Plot MK, RK, AK for WM, GM and CSF for NSA = 1, 2 in SAME PLOT (groups of 6)
fontsize = 20;

rk_DKI_WM1(2001:end,:) = NaN;
rk_DKI_WM2(2001:end,:) = NaN;

mk_pool = [mk_DKI_WM1 mk_DKI_WM2 mk_DKI_GM1 mk_DKI_GM2 mk_DKI_CSF1 mk_DKI_CSF2];
rk_pool = [rk_DKI_WM1 rk_DKI_WM2 rk_DKI_GM1 rk_DKI_GM2 rk_DKI_CSF1 rk_DKI_CSF2];
ak_pool = [ak_DKI_WM1 ak_DKI_WM2 ak_DKI_GM1 ak_DKI_GM2 ak_DKI_CSF1 ak_DKI_CSF2];


params_pool = [{mk_pool} {rk_pool} {ak_pool}]; % DKI
folder_pool = ['/Users/live/Documents/Masterthesis/MATLAB/DKI/Results_corrected/simulation_experiment_1/NSA=2/Pool/Plots_together'];
params_string = ['MK'; 'RK'; 'AK'];
group = [1 7 13 19 25 31 37 43 49 55 2 8 14 20 26 32 38 44 50 56 3 9 15 21 27 33 39 45 51 57 4 10 16 22 28 34 40 46 52 58 ...
    5 11 17 23 29 35 41 47 53 59 6 12 18 24 30 36 42 48 54 60];
positions = [1 1.25 1.5 1.75 2 2.25 3 3.25 3.5 3.75 4 4.25 5 5.25 5.5 5.75 6 6.25 7 7.25 7.5 7.75 8 8.25 ...
    9 9.25 9.5 9.75 10 10.25 11 11.25 11.5 11.75 12 12.25 13 13.25 13.5 13.75 14 14.25 15 15.25 15.5 15.75 16 16.25...
    17 17.25 17.5 17.75 18 18.25 19 19.25 19.50 19.75 20 20.25];


colorWM1 =  [0.98, 0.91, 0.71]; % banana mania
colorGM1 = [0.41, 0.41, 0.41]; % dim gray
colorCSF1 = [0.34, 0.63, 0.83]; % blue

colorWM2 =  [1, 0.60, 0.40]; % desert
colorGM2 = [0, 0, 0]; % black
colorCSF2 = [0, 0, 0.9]; % french gray

color = [{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},...
    {colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},...
    {colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1},{colorCSF2},{colorCSF1},{colorGM2},{colorGM1},{colorWM2},{colorWM1}];


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
    set(gca,'xticklabel',{'10','15','20','25','30','40','50','60','80','100'})


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
    set(gca, 'fontsize', fontsize)
    c = get(gca, 'Children');
    hleg1 = legend(c(1:6),'White matter, NSA = 1','White matter, NSA = 2',...
        'Gray matter, NSA = 1','Gray matter, NSA = 2','CSF, NSA = 1','CSF, NSA = 2');
    %legend('location','southeast');

    xlabel('SNR', 'fontsize', fontsize)
    ylabel('Relative error [%]', 'fontsize', fontsize);
    title([params_string(i,:)], 'fontsize', fontsize)
    
    %pause;
    w = 12;
    h = 6;
    set(gcf, 'PaperSize', [w h]);
    set(gcf, 'PaperPosition', [0 0 w h]);
    saveas(figure(i), fullfile(folder_pool, ['simexp2_DKI_' params_string(i,:) '_22']), 'pdf');

end
