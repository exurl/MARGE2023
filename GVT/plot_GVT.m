% Title: Plot GVT Data
% Author: Anthony Su
% Date: 2023-06-28

close all
clear all
clc

%% IMPORT DATA

% identify files to import
dir_files = dir;
files = [];
for idx = 1:length(dir_files)
    f = dir_files(idx);
    if(contains(f.name,'FRF_impact_') && endsWith(f.name,'.mat') && ~contains(f.name,'GVT1'))
        files = [files,string(f.name)];
    end
end
N_files = length(files);
clear dir_files f

% files = [];
% for in = {'2014z','2020z','2024z'}
%     for out = {'2014z','2020z','2024z'}
%         file = ['FRF_impact_',in{:},'_sens_',out{:},'_PITCHFREEGVT1.mat'];
%         files = [files,string(file)];
%     end
% end
% N_files = length(files)

% check files exist
if(any(~isfile(files)))   
    error('file does not exist!')
end

% initialize struct array
f = load(files(1));
FRFs_wingz1(1) = f;
FRFs_wingz2(1) = f;
FRFs_fusez(1) = f;
FRFs_minor(1) = f;
FRFs_twist(1) = f;
FRFs_long(1) = f;

% import
for idx_file = 1:N_files
    % import data object
    f = load(files(idx_file));
    
    % import objects into groups
    switch f.output
        case {"2214z","2314z","2224z","2324z"}
            FRFs_wingz1(end+1) = f;
        case {"2012z","2020z"}
            FRFs_wingz2(end+1) = f;
        case {"1013z","1028z","3208z","3308z"}
            FRFs_fusez(end+1) = f;
        case {"2214x","2224x","1013y","1028y"}
            FRFs_minor(end+1) = f;
        case {"2214z-2314z","2224z-2324z"}
            FRFs_twist(end+1) = f;
        case {"2014z","2020z","2024z"}
            % FRFs_long(end+1) = f;
        otherwise
            warning(['unexpected output location ',f.output])
    end
end

% delete initialization entry
FRFs_wingz1(1) = [];
FRFs_wingz2(1) = [];
FRFs_fusez(1) = [];
FRFs_minor(1) = [];
FRFs_twist(1) = [];
FRFs_long(1) = [];

% assemble groups
FRFs = {FRFs_wingz1,FRFs_wingz2,FRFs_fusez,FRFs_minor,FRFs_twist,FRFs_long};
FRFs_group_names = ["Wing Bend A","Wing Bend B","Fuselage Bend","In-Plane Bend","Wing Twist","Wing First Bend"];

%% PLOTTING

% identified reference frequencies
f0 = [1.40,2.9,10.25,18.4,19.7,32.6,60.7,36.6,74.9];

%
for idx = 1:length(FRFs)
    f_ = FRFs{idx};
    if(~isempty(f_))
        plot_FRFs(f_,FRFs_group_names(idx),f0)
    end
end
clear f_

%% FUNCTIONS

function plot_FRFs(FRFs,titleStr,f0)
    % INPUTS:
        % FRFs : struct array of structs with fields
            % coh       : 
            % data_from :
            % dt        :
            % freq      :
            % input     :
            % mag_H1    :
            % mag_H2    :
            % output    :
            % w         :
            % window    :
        % titleStr : string
        % f0 : 1-dim array of scalars
    
    N_FRFs = length(FRFs);

    % combine all FRFs into table
    freq = FRFs(1).freq;
    H1_FRF = table(freq);
    H2_FRF = table(freq);
    COH = table(freq);
    for idx_FRF = 1:N_FRFs
        FRF = FRFs(idx_FRF);
        name = [FRF.input,' to ',FRF.output];
        H1_FRF = addvars(H1_FRF,FRF.mag_H1,'NewVariableNames',{name});
        H2_FRF = addvars(H2_FRF,FRF.mag_H2,'NewVariableNames',{name});
        COH = addvars(COH,FRF.coh,'NewVariableNames',{name});
    end
    
    % plot FRF magnitude squared
    fig_mag = figure;
    SP_mag = stackedplot(H1_FRF,H2_FRF,XVariable='freq');
    title(['Magnitude: ',char(titleStr)])
    grid on
    xlabel('Frequency (Hz)')
    set(gcf().Children,'LegendLabels',{'H1 FRF','H2 FRF'})
    
    % % plot mode vertical lines
    % ax = findobj(SP_mag.NodeChildren,'Type','Axes');
    % for idx = 1:length(ax)
    %     xline(ax(idx),f0,':')
    %     set(ax(idx),'YGrid','off')
    %     set(ax(idx),'XMinorTick','on')
    % end

    % adjust figure size
    fig_mag.Position(4) = 34*(N_FRFs+2);
    fig_mag.Children(1).Position(2) = 1/(N_FRFs+2);
    fig_mag.Children(1).Position(4) = (N_FRFs+1)/(N_FRFs+2);

    % print plot
    % print(strcat('mag_',titleStr),'-dpdf','-vector','-fillpage')
    print(strcat('mag_',titleStr),'-dpng','-r300')



    % plot coherence
    fig_coh = figure;
    SP_coh = stackedplot(COH,XVariable='freq');
    title(['Coherence: ',char(titleStr)])
    grid on
    xlabel('Frequency (Hz)')
    
    % % plot mode vertical lines
    % ax = findobj(SP_coh.NodeChildren,'Type','Axes');
    % for idx = 1:length(ax)
    %     xline(ax(idx),f0,':')
    %     set(ax(idx),'YGrid','off')
    %     set(ax(idx),'XMinorTick','on')
    % end

    % adjust figure size
    fig_coh.Position(4) = 34*(N_FRFs+2);
    fig_coh.Children(1).Position(2) = 1/(N_FRFs+2);
    fig_coh.Children(1).Position(4) = N_FRFs/(N_FRFs+2);

    % print plot
    % print(strcat('coh_',titleStr),'-dpdf','-vector','-fillpage')
    print(strcat('coh_',titleStr),'-dpng','-r300')

end