% Title: Postprocess H1/H2 Study
% Author: Anthony Su
% Date: 2023-05-25

close all
clear all
clc

%% IMPORT DATA

% identify all files to import
impact_locations = {'2014z','2020z','2024z'};
sensor_locations = {'wingz'};
files = {};
for idx_impact = 1:length(impact_locations)
    for idx_sensor = 1:length(sensor_locations) % both input and output sensors
        files{end+1} = ['impact_',impact_locations{idx_impact},'_sens_',sensor_locations{idx_sensor},'_2023-05-23','.mat'];
    end
end

% check files exist
if(any(~isfile(files)))
    error('file does not exist!')
end

% import data
N_files = length(files);
N_tests = zeros(1,N_files);
N_sensors = zeros(N_files,1);
for idx_file = 1:N_files
    % import data object
    d_ = load(files{idx_file}).H;
    N_tests(idx_file) = length(d_);
    % extract sampling rate
    sampling_rate(idx_file) = d_(1).rate;
    % extract sensor data
    sensors{idx_file} = d_(1).sensors;
    N_sensors(idx_file) = length(sensors{idx_file})-1;
        % ^assuming sensor 1 is the input, this is the # of outputs
    for idx_test = 1:N_tests(idx_file)
        % extract impact index
        idx_impact(idx_file,idx_test) = d_(idx_test).peakIdxs(1);
        % extract time series as array
        data{idx_file,idx_test} = d_(idx_test).data{:,:};
        time{idx_file,idx_test} = seconds(d_(idx_test).data.Properties.RowTimes); % s
    end
end

%% TRUNCATION

for idx_file = 1:N_files
    for idx_test = 1:N_tests(idx_file)
        % truncate any prior to t-0.01s and t+3.70s from impact peak
        [~,idx_max] = max(data{idx_file,idx_test}(:,1));
        idx1 = idx_max-512;
        idx2 = floor(idx_max+51200*3.7);
        data{idx_file,idx_test} = data{idx_file,idx_test}(idx1:idx2,:);

        time{idx_file,idx_test} = time{idx_file,idx_test}(idx1:idx2);
        time{idx_file,idx_test} = time{idx_file,idx_test}-time{idx_file,idx_test}(1);
    end
end

%% REMOVE BAD DATA

% remove impact 2024z test #5
data{3,5} = [];
time{3,5} = [];
N_tests(3) = 4;

%% AVERAGING

% for each file, average each sensor's data across tests NOT by taking
% their average and computing the FRF, but by concatenating them and taking
% the FRF.
for idx_file = 1:N_files
    data_avg{idx_file} = [];
    for idx_test = 1:N_tests(idx_file)
        data_avg{idx_file} = [data_avg{idx_file}; data{idx_file,idx_test}];
    end
end

%% COMPUTATION

% compute frequency response
for idx_file = 1:N_files
    for idx_test = 1:N_tests(idx_file)
        % pull input data
        u = data{idx_file,idx_test}(:,1);
        for idx_sensor = 1:N_sensors(idx_file) % output sensors only
            % pull output data
            y = data{idx_file,idx_test}(:,idx_sensor+1); % skip input at idx=1
            % subtract offset
            u = u-mean(u);
            y = y-mean(y);
            % compute
            w = [1.25,100]; % <-- FREQUENCY BANDWIDTH FOR CZT FRF
            [FRF{idx_test,idx_sensor,idx_file},...
                f{idx_test,idx_sensor,idx_file},...
                coh{idx_test,idx_sensor,idx_file}] = ...
                czt_FRF(u,y,[],w,'off',sampling_rate(idx_file));
        end
    end
    % do also for the averaged data
    u = data_avg{idx_file}(:,1);
    for idx_sensor = 1:N_sensors(idx_file) % output sensors only
        % pull output data
        y = data_avg{idx_file}(:,idx_sensor+1); % skip input at idx=1
        % subtract offset
        u = u-mean(u);
        y = y-mean(y);
        % compute
        w = [1.2,100];
        [~,f_avg{idx_sensor,idx_file},...
            coh_avg{idx_sensor,idx_file},...
            Gxx_hat_avg{idx_sensor,idx_file},...
            Gyy_hat_avg{idx_sensor,idx_file},...
            Gxy_hat_avg{idx_sensor,idx_file},...
            Gyx_hat_avg{idx_sensor,idx_file}] = ...
            czt_FRF(u,y,[],w,'off',sampling_rate(idx_file));
    end
end

%% PLOTTING

% plot frequency response of all impacts, overlaying trials
figure
TL = tiledlayout(N_files,max(N_sensors),'TileSpacing','compact','Padding','tight');

for idx_file = 1:N_files
    for idx_sensor = 1:N_sensors(idx_file)
        % pull data
        f_ = f(:,idx_sensor,idx_file);
        FRF_ = FRF(:,idx_sensor,idx_file);
        coh_ = coh(:,idx_sensor,idx_file);

        f_avg_ = f_avg{idx_sensor,idx_file};
        % FRF_avg_ = FRF_avg{idx_sensor,idx_file};
        FRF_avg_H1_ = Gxy_hat_avg{idx_sensor,idx_file}./Gxx_hat_avg{idx_sensor,idx_file};
        FRF_avg_H2_ = Gyy_hat_avg{idx_sensor,idx_file}./Gyx_hat_avg{idx_sensor,idx_file};
        coh_avg_ = coh_avg{idx_sensor,idx_file};

        ylabels = {'$\left|\frac{y}{u}\right|$ (dB)','$\angle\frac{y}{u}$ (deg)','$\mathrm{coh}\left(u,y\right)$'};

        % sp = stackedplot(T,f_,[mag2db(abs(FRF_)),rad2deg(phase(FRF_)),coh_],'DisplayLabels',ylabels);
        % sp.Layout.Tile = (idx_test-1)*N_sensors(idx_file)+idx_sensor;

        % create a mag/phase/cohere plot for each sensor in each test
        tl = tiledlayout(TL,3,1,'TileSpacing','compact','Padding','tight');
        % place the nested tiled layout in the correct tile of the parent tiled layout
        tl.Layout.Tile = (idx_file-1)*N_sensors(idx_file)+idx_sensor;
        
        % plot averaged data
        nexttile(tl,1)
        plot(f_avg_,mag2db(abs(FRF_avg_H1_)),'k','LineWidth',2)
        ylabel(ylabels{1},'Interpreter','latex','Rotation',0)
        ax = gca;
        ax.XTickLabel = {};
        hold on
        plot(f_avg_,mag2db(abs(FRF_avg_H2_)),'r','LineWidth',1)
        % plot phase
        nexttile(tl,2)
        plot(f_avg_,rad2deg(phase(FRF_avg_H1_)),'k','LineWidth',2)
        ylabel(ylabels{2},'Interpreter','latex','Rotation',0)
        ax = gca;
        ax.XTickLabel = {};
        hold on
        plot(f_avg_,rad2deg(phase(FRF_avg_H2_)),'r','LineWidth',1)
        % plot coherence
        nexttile(tl,3)
        plot(f_avg_,coh_avg_,'k','LineWidth',2)
        ylabel(ylabels{3},'Interpreter','latex','Rotation',0)
        hold on
        plot(f_avg_,coh_avg_,'r','LineWidth',1)
        
        % % plot individual tests
        % for idx_test = 1:N_tests(idx_file)
        %     % plot magnitude
        %     nexttile(tl,1)
        %     plot(f_{idx_test},mag2db(abs(FRF_{idx_test})))
        %     % plot phase
        %     nexttile(tl,2)
        %     plot(f_{idx_test},rad2deg(phase(FRF_{idx_test})))
        %     % plot coherence
        %     nexttile(tl,3)
        %     plot(f_{idx_test},coh_{idx_test})
        % end
        % title nested tiled layout
        sensLoc = sensors{idx_file}(idx_sensor+1).Location{:};
        hammerLoc = impact_locations{idx_file};
        title(tl,['sensor ',sensLoc,' hammer ',hammerLoc])

        legend({'H1','H2'},'Location','southeast')
    end
end
xlabel(TL,'Frequency (Hz)','FontSize',18)

%%

% % plot frequency response of every single impact
% for idx_file = 1:N_files
%     % create a figure for each file
%     figure
%     TL = tiledlayout(N_tests(idx_file),max(N_sensors(idx_file)),'TileSpacing','compact','Padding','tight');
% 
%     for idx_test = 1:N_tests(idx_file)
%         for idx_sensor = 1:N_sensors(idx_file)
%             % pull data
%             f_ = f{idx_test,idx_sensor,idx_file};
%             FRF_ = FRF{idx_test,idx_sensor,idx_file};
%             coh_ = coh{idx_test,idx_sensor,idx_file};
% 
%             ylabels = {'$\left|\frac{y}{u}\right|$','$\angle\frac{y}{u}$','$coh\left(u,y\right)$'};
% 
%             % sp = stackedplot(T,f_,[mag2db(abs(FRF_)),rad2deg(phase(FRF_)),coh_],'DisplayLabels',ylabels);
%             % sp.Layout.Tile = (idx_test-1)*N_sensors(idx_file)+idx_sensor;
% 
%             % create a mag/phase/cohere plot for each sensor in each test
%             tl = tiledlayout(TL,3,1,'TileSpacing','compact','Padding','tight');
%             % place the nested tiled layout in the correct tile of the parent tiled layout
%             tl.Layout.Tile = (idx_test-1)*N_sensors(idx_file)+idx_sensor;
%             % plot magnitude
%             nexttile(tl,1)
%             plot(f_,mag2db(abs(FRF_)))
%             ylabel(ylabels{1},'Interpreter','latex','Rotation',0)
%             ax = gca;
%             ax.XTickLabel = {};
%             hold on
%             % plot phase
%             nexttile(tl,2)
%             plot(f_,rad2deg(phase(FRF_)))
%             ylabel(ylabels{2},'Interpreter','latex','Rotation',0)
%             ax = gca;
%             ax.XTickLabel = {};
%             hold on
%             % plot coherence
%             nexttile(tl,3)
%             plot(f_,coh_)
%             ylabel(ylabels{3},'Interpreter','latex','Rotation',0)
%             hold on
%             % title nested tiled layout
%             sensLoc = sensors{idx_file}(idx_sensor+1).Location{:};
%             title(tl,sensLoc)
%         end
%     end
%     xlabel(TL,'Frequency (Hz)')
%     hammerLoc = impact_locations{idx_file};
%     title(TL,['Impact Hammer at ',hammerLoc])
% end

%%

% % plot raw response of every impact
% for idx_file = 1:N_files
%     % create a figure for each file
%     figure
%     T = tiledlayout(N_tests(idx_file),max(N_sensors(idx_file)),'TileSpacing','compact','Padding','tight');
% 
%     for idx_test = 1:N_tests(idx_file)
%         for idx_sensor = 1:N_sensors(idx_file)
%             % pull data
%             t_ = time{idx_file,idx_test};
%             d_ = data{idx_file,idx_test};
% 
%             tl = tiledlayout(T,2,1,'TileSpacing','compact','Padding','tight');
%             % place the nested tiled layout in the correct tile of the parent tiled layout
%             tl.Layout.Tile = (idx_test-1)*N_sensors(idx_file)+idx_sensor;
%             % plot magnitude
%             nexttile(tl,1)
%             plot(t_,d_(:,1))
%             ylabel('forcing (N)')
%             ax = gca;
%             ax.XTickLabel = {};
%             hold on
%             % plot phase
%             nexttile(tl,2)
%             plot(t_,d_(:,2:4))
%             ylabel('response (m/s2)')
%             ax = gca;
%             hold on
%             % title nested tiled layout
%             sensLoc = sensors{idx_file}(idx_sensor+1).Location{:};
%             title(tl,sensLoc)
% 
%             xlim('tight')
%         end
%     end
%     xlabel(T,'Time (s)')
%     hammerLoc = impact_locations{idx_file};
%     title(T,['Impact Hammer at ',hammerLoc])
% end
