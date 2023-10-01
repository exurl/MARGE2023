% TITLE: Plot Magnitude, All Modes
% AUTHOR: Anthony Su
% DATE: 2023-05-31

close all
clear all
clc

%% import and open figures
f(1) = openfig('wingz_200.fig');
f(2) = openfig('wingt_200.fig');
f(3) = openfig('fusezwingx_200.fig');
f(4) = openfig('fuseytailzReduced_200.fig');

%% extract data from figures
freqs_data = {};
mags_data = {};
titles = {};
labels = {};
% for each figure
for idx_f = 1:length(f)
    f_ = f(idx_f);
    objs = f_.Children.Children;
    plots = objs(isgraphics(objs,'Axes')); % take axes objects
    legs = objs(isgraphics(objs,'Legend'));
    leg = legs(1).String;
    plots = plots(3:3:end); % take only axes w/ magnitude plots
    freqs_data{idx_f} = plots(1).Children(1).XData; % frequency vector
    for idx_p = 1:length(plots)
        p_ = plots(idx_p);
        lines = p_.Children;
        titles{idx_p,idx_f} = string(p_.Title.String); % title (hammer)
        for idx_l = 1:length(lines)
            l_ = lines(idx_l);
            mags_data{idx_l,idx_p,idx_f} = l_.YData; % magnitude vector
            labels{idx_l,idx_f} = string(leg{idx_l}); % line label (sensor)
        end
    end
end

%% re-plot data

fig = figure;
tl = tiledlayout('vertical','TileSpacing','compact','Padding','tight');
colors = {[0,0,1],[0,0.8,0],[1,0,0]};
lineStyles = {'-','-.',':'};
for idx_f = 1:length(freqs_data)
    x = freqs_data{idx_f};
    for idx_p = 1:length(titles(:,idx_f))
        if(~isempty(titles{idx_p,idx_f}))
            ax = nexttile(tl);
            for idx_l = 1:length(labels(:,idx_f))
                y = mags_data{idx_l,idx_p,idx_f};
                if(~isempty(y))
                    plot(x,y,'LineWidth',1.5,'Color',colors{idx_l},'LineStyle',lineStyles{idx_l})
                    hold on
                end
            end
            set(ax,'XTickLabels',[])
            legend(ax,[labels{:,idx_f}],'Location','eastoutside')
            ylabel(ax,titles{idx_p,idx_f},'Rotation',0)
            grid on
            set(ax,'XMinorGrid',true) %
        end
    end
end
linkaxes(tl.Children,'x')
set(ax,'XTickLabels',0:10:100) % re-enable last plot's XTick
ylabel(tl,'Response Magnitude Squared')
xlabel(tl,'Frequency (Hz)')

fig.Position(1) = 0;
fig.Position(2) = 0;
fig.Position(3) = 1600;
fig.Position(4) = 1200;
% print(fig,'-vector','powerPlots','-dsvg')