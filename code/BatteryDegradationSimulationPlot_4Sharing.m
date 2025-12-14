% ============================================================
% Reproduce figures using battery degradation simulation data
% Path of the data is included in the script. 
% ============================================================

clear; clc;

csvDir = "../data/BatteryDegradationSimulationData_csv";
outfile_dir = '../results';

colUC = [0 0.22 0.8];
colV2 = [0.8 0 0];


%% ===============================
% FIGURE — 71 kWh (2×3)
% Row: NMC, LFP
% COlummn: 'Overall capacity','Capacity (cycling)','Capacity (calendar)'
% ===============================
fig1 = figure('Units','centimeters','Position',[2 2 18 12],'Color','w');
tl = tiledlayout(fig1,2,3,'Padding','compact','TileSpacing','compact');

metrics = {'overall','cycling','calendar'};
ylabs   = {'Overall capacity','Capacity (cycling)','Capacity (calendar)'};
chems   = {'NMC','LFP'};

for r = 1:2
    chem = chems{r};
    for c = 1:3
        metric = metrics{c};

        A = readtable(fullfile(csvDir, sprintf("71kWh_%s_%s_UncontrolledCharging.csv", chem, metric)));
        B = readtable(fullfile(csvDir, sprintf("71kWh_%s_%s_V2G.csv",                 chem, metric)));

        ax = nexttile(tl,(r-1)*3+c); hold(ax,'on');
        hA = plot_band_from_summary(ax, A, true, colUC, 'Uncontrolled charging');
        hB = plot_band_from_summary(ax, B, true, colV2, 'V2G');

        yline(ax,0.8,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8);
        ylim(ax,[0 1]); xlim(ax,[A.t_years(1) A.t_years(end)]);
        xlabel(ax,'Time (years)'); ylabel(ax, ylabs{c});
        title(ax, sprintf("%s — 71 kWh", chem));

        if r==1 && c==1
            legend(ax,[hA hB],{'Uncontrolled charging','V2G'},'Location','southwest','Box','off');
        end
    end
end

%{  %}
if ~isempty(outfile_dir)
    if ~exist(outfile_dir,'dir'), mkdir(outfile_dir); end
    exportgraphics(fig1, fullfile(outfile_dir, "BatteryDegradation_71kWh_NMC_LFP"+...
           "_" + string(datetime("now", "Format", "yyyy-MM-dd")) + ".pdf"), ...
           'ContentType','vector','BackgroundColor','none') ;
end



function hLegend = plot_band_from_summary(ax, T, doEOL, colLine, legendLabel)

    t   = T.t_years(:);
    p10 = T.p10(:);
    p90 = T.p90(:);
    mu  = T.median(:);

    % band
    x_band = [t; flipud(t)];
    y_band = [p10; flipud(p90)];
    fill(ax, x_band, y_band, colLine, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
    hold(ax,'on');

    if doEOL
        eol_level = 0.8;
        cross_idx = find(mu(1:end-1) >= eol_level & mu(2:end) < eol_level, 1, 'first');

        if ~isempty(cross_idx)
            t1 = t(cross_idx);   t2 = t(cross_idx+1);
            y1 = mu(cross_idx);  y2 = mu(cross_idx+1);
            if y2 ~= y1
                x_eol = t1 + (eol_level - y1) * (t2 - t1) / (y2 - y1);
            else
                x_eol = t1;
            end

            mask_before = t <= x_eol;
            plot(ax, t(mask_before),  mu(mask_before),  'Color', colLine, 'LineWidth', 1.6);
            plot(ax, t(~mask_before), mu(~mask_before), ':',     'Color', colLine, 'LineWidth', 1.6);

            plot(ax, x_eol, eol_level, 'o', 'MarkerSize', 5, ...
                'MarkerFaceColor', colLine, 'MarkerEdgeColor', 'w');

            dx = 0.02 * range(t) - 5;
            dy = 0.1;
            if contains(legendLabel, "V2G"); dy = -0.1; end

            text(ax, x_eol + dx, eol_level + dy, sprintf('EOL = %.2f', x_eol), ...
                'FontSize', 7, 'Color', colLine, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        else
            plot(ax, t, mu, 'Color', colLine, 'LineWidth', 1.6);
        end
    else
        plot(ax, t, mu, 'Color', colLine, 'LineWidth', 1.6);
    end

    hLegend = plot(ax, NaN, NaN, '-', 'Color', colLine, 'LineWidth', 1.6, 'DisplayName', legendLabel);
end
