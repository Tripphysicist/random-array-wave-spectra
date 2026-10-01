function reproduce_all_figures()
% REPRODUCE_ALL_FIGURES
% Master script to generate Figures 1 through 13 for the paper:
% "Measuring Directional Wave Spectra with Random Spatial Arrays"
% (Collins et al., 2026)
%
% Generates publication figures and saves them in ./figures/ directory.

close all;
repoDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(repoDir, 'src')));
figDir = fullfile(repoDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

dataBase = fullfile(repoDir, 'data');

fprintf('=========================================================\n');
fprintf(' Reproducing Figures 1 - 13 for Random Wave Arrays Paper\n');
fprintf('=========================================================\n\n');

%% --- FIGURE 1: Directional Distribution vs. Scale (Tests 1-11) ---
fprintf('[1/13] Generating Figure 1 (Directional Distributions vs Scale)...\n');
dataFolder = fullfile(dataBase, 'targetL2', 'N4');
datadir = dir(fullfile(dataFolder, '*.mat'));
numScales = numel(datadir);

fig1 = figure(1); clf;
set(fig1, 'Color', 'w', 'Position', [100 100 1200 800]);

for loopi = 1:numScales
    d = load(fullfile(dataFolder, datadir(loopi).name), ...
        'targetL2', 'lambdap', 'parametricDistribution', ...
        'sampledDirectionalSpectrum', 'sampledDirectionConf');

    scale_val = sqrt(d.targetL2) / d.lambdap;
    theta = d.sampledDirectionalSpectrum.theta(:);
    trueDist = d.parametricDistribution(:);

    meanAcrossWaves = squeeze(mean(d.sampledDirectionConf, 2, 'omitnan'));
    
    subplot(3, 4, loopi);
    dist_mu_raw = mean(meanAcrossWaves, 2, 'omitnan');
    dist_sig_raw = std(meanAcrossWaves, 0, 2, 'omitnan');

    shadedErrorBar(theta, dist_mu_raw, dist_sig_raw, 'lineProps', {'b', 'LineWidth', 1.5});
    hold on;
    plot(theta, trueDist, 'k--', 'LineWidth', 2.5);

    grid on;
    axis([-pi pi 0 0.07]);
    xlabel('Direction [rad]');
    ylabel('S(\theta) [m^2/rad]');

    text(-2.8, 0.063, sprintf('Test %d: l_{rms}/\\lambda_p = 10^{%0.1f}', loopi, log10(scale_val)), ...
        'Interpreter', 'tex', 'FontSize', 11, 'FontWeight', 'bold');
    set(gca, 'FontSize', 12);
end
packfig(3,4);
exportgraphics(fig1, fullfile(figDir, 'Figure1.png'), 'Resolution', 200);
exportgraphics(fig1, fullfile(figDir, 'Figure1.eps'), 'ContentType', 'vector');


%% --- FIGURE 2: Error Metrics vs. Scale (N=4 and N=10) ---
fprintf('[2/13] Generating Figure 2 (Performance vs Scale)...\n');
dataFolder_N4 = fullfile(dataBase, 'targetL2', 'N4');
datadir_N4 = dir(fullfile(dataFolder_N4, '*.mat'));
numFiles_N4 = numel(datadir_N4);

dataFolder_N10 = fullfile(dataBase, 'targetL2', 'N10');
datadir_N10 = dir(fullfile(dataFolder_N10, '*.mat'));
numFiles_N10 = numel(datadir_N10);

target_scale_vec_N10 = zeros(numFiles_N10, 1);
peak_std_mean_N10 = zeros(numFiles_N10, 1);
peak_std_err_N10  = zeros(numFiles_N10, 1);
iae_mean_N10 = zeros(numFiles_N10, 1);
iae_std_N10  = zeros(numFiles_N10, 1);

for idx = 1:numFiles_N10
    filePath = fullfile(dataFolder_N10, datadir_N10(idx).name);
    d10 = load(filePath);
    if isfield(d10, 'lambdap'), Lp10 = d10.lambdap; else, Lp10 = (9.81 * d10.Tp^2) / (2 * pi); end
    target_scale_vec_N10(idx) = sqrt(d10.targetL2) / Lp10;
    theta10 = d10.sampledDirectionalSpectrum.theta(:);
    trueDist10 = d10.parametricDistribution(:);
    Area_true10 = trapz(theta10, trueDist10);
    trueDistNorm10 = trueDist10 ./ Area_true10;

    [~, maxIdx10] = max(d10.sampledDirectionConf, [], 1);
    peakThetas10 = theta10(maxIdx10) * (180/pi);
    std_per_geometry10 = squeeze(std(peakThetas10, 0, 2));
    peak_std_mean_N10(idx) = mean(std_per_geometry10, 'omitnan');
    peak_std_err_N10(idx)  = std(std_per_geometry10, 'omitnan');

    meanAcrossWaves10 = squeeze(mean(d10.sampledDirectionConf, 2, 'omitnan'));
    numConfigs10 = size(meanAcrossWaves10, 2);
    iae_per_geometry10 = zeros(numConfigs10, 1);
    for k = 1:numConfigs10
        dist_raw_k = meanAcrossWaves10(:, k);
        if any(isnan(dist_raw_k)), iae_per_geometry10(k) = NaN; continue; end
        Area_k = trapz(theta10, dist_raw_k);
        dist_norm_k = dist_raw_k ./ Area_k;
        iae_per_geometry10(k) = trapz(theta10, abs(dist_norm_k - trueDistNorm10));
    end
    iae_mean_N10(idx) = mean(iae_per_geometry10, 'omitnan');
    iae_std_N10(idx)  = std(iae_per_geometry10, 'omitnan');
end

target_scale_vec_N4 = zeros(numFiles_N4, 1);
peak_std_mean_N4 = zeros(numFiles_N4, 1);
peak_std_err_N4  = zeros(numFiles_N4, 1);
iae_mean_N4 = zeros(numFiles_N4, 1);
iae_std_N4  = zeros(numFiles_N4, 1);

fig2 = figure(2); clf;
set(fig2, 'Color', 'w', 'Position', [100 100 1000 800]);
ax1 = subplot(2, 1, 1); hold on; grid on;
ax2 = subplot(2, 1, 2); hold on; grid on;

clr = colororder;
for i = 1:numFiles_N4
    filePath = fullfile(dataFolder_N4, datadir_N4(i).name);
    data = load(filePath);
    if isfield(data, 'lambdap'), Lp = data.lambdap; else, Lp = (9.81 * data.Tp^2) / (2 * pi); end
    actual_lrms_norm = sqrt(data.l2) / Lp;
    target_scale_vec_N4(i) = sqrt(data.targetL2) / Lp;

    [~, maxIdx] = max(data.sampledDirectionConf, [], 1);
    peakThetas = data.sampledDirectionalSpectrum.theta(maxIdx) * (180/pi);
    std_per_geometry = squeeze(std(peakThetas, 0, 2));
    peak_std_mean_N4(i) = mean(std_per_geometry, 'omitnan');
    peak_std_err_N4(i)  = std(std_per_geometry, 'omitnan');

    thisColor = clr(mod(i-1, size(clr,1))+1, :);
    scatter(ax1, actual_lrms_norm, std_per_geometry, 85, ...
        'MarkerFaceColor', thisColor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');

    trueDist = data.parametricDistribution(:);
    theta = data.sampledDirectionalSpectrum.theta(:);
    Area_true = trapz(theta, trueDist);
    trueDistNorm = trueDist ./ Area_true;

    meanAcrossWaves = squeeze(mean(data.sampledDirectionConf, 2, 'omitnan'));
    numConfigs = size(meanAcrossWaves, 2);
    iae_per_geometry = zeros(numConfigs, 1);
    for k = 1:numConfigs
        dist_raw_k = meanAcrossWaves(:, k);
        if any(isnan(dist_raw_k)), iae_per_geometry(k) = NaN; continue; end
        Area_k = trapz(theta, dist_raw_k);
        dist_norm_k = dist_raw_k ./ Area_k;
        iae_per_geometry(k) = trapz(theta, abs(dist_norm_k - trueDistNorm));
    end
    iae_mean_N4(i) = mean(iae_per_geometry, 'omitnan');
    iae_std_N4(i)  = std(iae_per_geometry, 'omitnan');

    scatter(ax2, actual_lrms_norm, iae_per_geometry, 85, ...
        'MarkerFaceColor', thisColor, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
end

[target_sorted, sIdx] = sort(target_scale_vec_N4);
[target_sorted_N10, sIdx_N10] = sort(target_scale_vec_N10);

errorbar(ax1, target_sorted, peak_std_mean_N4(sIdx), peak_std_err_N4(sIdx), ...
    'k-o', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'CapSize', 8, 'DisplayName', 'N=4');
errorbar(ax1, target_sorted_N10 * 1.1, peak_std_mean_N10(sIdx_N10), peak_std_err_N10(sIdx_N10), ...
    'Color', [0.7 0.7 0.7], 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.7 0.7 0.7], 'CapSize', 8, 'DisplayName', 'N=10');

errorbar(ax2, target_sorted, iae_mean_N4(sIdx), iae_std_N4(sIdx), ...
    'k-s', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'CapSize', 8, 'DisplayName', 'N=4');
errorbar(ax2, target_sorted_N10 * 1.1, iae_mean_N10(sIdx_N10), iae_std_N10(sIdx_N10), ...
    'Color', [0.7 0.7 0.7], 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.7 0.7 0.7], 'CapSize', 8, 'DisplayName', 'N=10');

set(ax1, 'XScale', 'log', 'FontSize', 14);
ylabel(ax1, '\sigma(\theta_p) [^\circ]');
title(ax1, 'Peak Direction Precision');
xlim(ax1, [10^-4.5 10^1.5]); grid(ax1, 'on'); box(ax1, 'on');
legend(ax1, 'Location', 'best');

set(ax2, 'XScale', 'log', 'FontSize', 14);
ylabel(ax2, 'MIAE [rad^{-1}]');
xlabel(ax2, 'Target Scale ($l_{rms}/\lambda_p$)', 'Interpreter', 'latex');
title(ax2, 'Total Distribution Accuracy');
xlim(ax2, [10^-4.5 10^1.5]); ylim(ax2, [0 1.25]);
grid(ax2, 'on'); box(ax2, 'on');
legend(ax2, 'Location', 'southeast');

for j = 1:length(target_sorted)
    text(ax2, target_sorted(j), 1.05, num2str(j), ...
        'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'BackgroundColor', [1 1 1 0.75], 'Margin', 1);
end
exportgraphics(fig2, fullfile(figDir, 'Figure2.png'), 'Resolution', 200);
exportgraphics(fig2, fullfile(figDir, 'Figure2.eps'), 'ContentType', 'vector');


%% --- HELPER FOR CONFIGS & CO-ARRAY FIGURES (Figures 3, 4, 5, 8) ---
plot_config_coarray = @(dataFile, figNum, figTitle, outName) helper_config_coarray(dataFile, figNum, figTitle, outName, figDir);

fprintf('[3/13] Generating Figure 3 (Test 4 Configurations & Co-arrays)...\n');
plot_config_coarray(fullfile(dataBase, 'targetL2', 'N4', 'sensors04_scale_04_ensembles30_configs100.mat'), 3, 'Test 4 (l^* = 10^{-2.5})', 'Figure3');

fprintf('[4/13] Generating Figure 4 (Test 8 Configurations & Co-arrays)...\n');
plot_config_coarray(fullfile(dataBase, 'targetL2', 'N4', 'sensors04_scale_08_ensembles30_configs100.mat'), 4, 'Test 8 (l^* = 10^{-0.5})', 'Figure4');

fprintf('[5/13] Generating Figure 5 (Test 7 Configurations & Co-arrays)...\n');
plot_config_coarray(fullfile(dataBase, 'targetL2', 'N4', 'sensors04_scale_07_ensembles30_configs100.mat'), 5, 'Test 7 (l^* = 10^{-1.0})', 'Figure5');


%% --- FIGURE 6: Error Metrics vs. N (N=3..20) ---
fprintf('[6/13] Generating Figure 6 (Metrics vs Sensor Count N)...\n');
dataFolder_N = fullfile(dataBase, 'N');
datadir_N = dir(fullfile(dataFolder_N, '*.mat'));
numFiles_N = numel(datadir_N);

n_vec = zeros(numFiles_N, 1);
iae_mean_N = zeros(numFiles_N, 1);
iae_std_N  = zeros(numFiles_N, 1);
peak_std_mean_N = zeros(numFiles_N, 1);
peak_std_err_N  = zeros(numFiles_N, 1);

for i = 1:numFiles_N
    d = load(fullfile(dataFolder_N, datadir_N(i).name));
    tokens = regexp(datadir_N(i).name, 'sensors(\d+)', 'tokens');
    if ~isempty(tokens)
        n_vec(i) = str2double(tokens{1}{1});
    elseif isfield(d, 'x') && ~isempty(d.x)
        n_vec(i) = size(d.x, 1);
    else
        n_vec(i) = i;
    end

    theta = d.sampledDirectionalSpectrum.theta(:);
    trueDist = d.parametricDistribution(:);

    [~, maxIdx] = max(d.sampledDirectionConf, [], 1);
    peakThetas = theta(maxIdx) * (180/pi);
    std_per_wave = squeeze(std(peakThetas, 0, 3));
    peak_std_mean_N(i) = mean(std_per_wave, 'omitnan');
    peak_std_err_N(i)  = std(std_per_wave, 'omitnan');

    absErr = abs(d.sampledDirectionConf - trueDist);
    dTh = abs(theta(2) - theta(1));
    mean_err_per_wave = squeeze(mean(sum(absErr, 1) * dTh, 3));
    iae_mean_N(i) = mean(mean_err_per_wave, 'omitnan');
    iae_std_N(i)  = std(mean_err_per_wave, 'omitnan');
end

[n_sorted, sIdx] = sort(n_vec);
fig6 = figure(6); clf;
set(fig6, 'Color', 'w', 'Position', [150 150 800 600]);

subplot(2, 1, 1);
errorbar(n_sorted, peak_std_mean_N(sIdx), peak_std_err_N(sIdx), 'o-k', ...
    'LineWidth', 2, 'MarkerFaceColor', 'k', 'CapSize', 8);
grid on; xlim([2 21]); ylabel('\sigma(\theta_p) [^\circ]');
set(gca, 'XTick', n_sorted, 'FontSize', 12);
title('Peak Direction Precision vs. Sensor Count N');

subplot(2, 1, 2);
errorbar(n_sorted, iae_mean_N(sIdx), iae_std_N(sIdx), 's-k', ...
    'LineWidth', 2, 'MarkerFaceColor', 'k', 'CapSize', 8);
grid on; xlim([2 21]); ylabel('MIAE [rad^{-1}]'); xlabel('Number of Sensors (N)');
set(gca, 'XTick', n_sorted, 'FontSize', 12);
title('Total Distribution Accuracy vs. Sensor Count N');
packfig(2,1);

exportgraphics(fig6, fullfile(figDir, 'Figure6.png'), 'Resolution', 200);
exportgraphics(fig6, fullfile(figDir, 'Figure6.eps'), 'ContentType', 'vector');


%% --- FIGURE 7: Directional Distributions vs. N ---
fprintf('[7/13] Generating Figure 7 (Directional Distributions vs N)...\n');
fig7 = figure(7); clf;
set(fig7, 'Color', 'w', 'Position', [100 100 1200 800]);

for i = 1:numFiles_N
    d = load(fullfile(dataFolder_N, datadir_N(sIdx(i)).name));
    theta = d.sampledDirectionalSpectrum.theta(:);
    trueDist = d.parametricDistribution(:);

    subplot(3, 4, i);
    meanAcrossWaves = squeeze(mean(d.sampledDirectionConf, 2, 'omitnan'));
    mu_dist = mean(meanAcrossWaves, 2, 'omitnan');
    sig_dist = std(meanAcrossWaves, 0, 2, 'omitnan');

    shadedErrorBar(theta, mu_dist, 2*sig_dist, 'lineProps', {'b', 'LineWidth', 1.5});
    hold on;
    plot(theta, trueDist, 'k--', 'LineWidth', 2.0);
    grid on;
    axis([-pi pi 0 0.08]);
    xlabel('\theta [rad]'); ylabel('E(\theta) [m^2/rad]');
    text(-2.8, 0.07, sprintf('N = %d', n_sorted(i)), 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'FontSize', 11);
end
packfig(3,4);
exportgraphics(fig7, fullfile(figDir, 'Figure7.png'), 'Resolution', 200);
exportgraphics(fig7, fullfile(figDir, 'Figure7.eps'), 'ContentType', 'vector');


%% --- FIGURE 8: 3-Element Arrays Configurations & Co-arrays ---
fprintf('[8/13] Generating Figure 8 (N=3 Configurations & Co-arrays)...\n');
plot_config_coarray(fullfile(dataBase, 'N', 'sensors03__targetL2_0.3157_ensembles30_configs30.mat'), 8, 'N = 3 Array Geometries', 'Figure8');


%% --- FIGURE 9: Sensitivity to Noise ---
fprintf('[9/13] Generating Figure 9 (Noise Sensitivity)...\n');
plot_sensitivity(fullfile(dataBase, 'sensitivity', 'noise'), 9, 'Noise Level', ...
    {'0.1%', '0.5%', '5.0%', '10.0%', '50.0%', '100.0%'}, 'Figure9', figDir);


%% --- FIGURE 10: Sensitivity to Time Lag ---
fprintf('[10/13] Generating Figure 10 (Time Lag Sensitivity)...\n');
plot_sensitivity(fullfile(dataBase, 'sensitivity', 'lag'), 10, 'Time Lag Level', ...
    {'0.1% T_p', '0.5% T_p', '1.0% T_p', '5.0% T_p', '10.0% T_p'}, 'Figure10', figDir);


%% --- FIGURE 11: Sensitivity to Spatial Offset ---
fprintf('[11/13] Generating Figure 11 (Spatial Offset Sensitivity)...\n');
plot_sensitivity(fullfile(dataBase, 'sensitivity', 'offset'), 11, 'Spatial Offset Level', ...
    {'0.1% \lambda_p', '0.5% \lambda_p', '1.0% \lambda_p', '5.0% \lambda_p', '10.0% \lambda_p'}, 'Figure11', figDir);


%% --- FIGURE 12: Estimator Comparison (N=7, 5 Methods) ---
fprintf('[12/13] Generating Figure 12 (Estimator Distributions Baseline N=7)...\n');
mKeys = {'BDM', 'EMEM', 'EWDM', 'IMLM', 'bdm_nnls'};
mNames = {'BDM-NR', 'EMEM', 'EWDM', 'IMLM', 'BDM-NNLS'};
mColors = [0.0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.494 0.184 0.556; 0.301 0.745 0.933];

fig12 = figure(12); clf;
set(fig12, 'Color', 'w', 'Position', [100 100 1100 700]);
dirN7 = fullfile(dataBase, 'methods', 'N7');
filesN7 = dir(fullfile(dirN7, '*.mat'));

for m = 1:5
    for f = 1:numel(filesN7)
        if startsWith(filesN7(f).name, mKeys{m}, 'IgnoreCase', true)
            d = load(fullfile(dirN7, filesN7(f).name));
            break;
        end
    end
    
    theta = d.sampledDirectionalSpectrum.theta(:);
    if isfield(d, 'parametricDirectionalSpectrum') && isfield(d.parametricDirectionalSpectrum, 'S')
        trueDist = trapz(d.parametricDirectionalSpectrum.w, d.parametricDirectionalSpectrum.S, 2);
    else
        trueDist = d.parametricDistribution(:);
    end
    
    meanAcrossWaves = squeeze(mean(d.sampledDirectionConf, 2, 'omitnan'));
    mu_dist = mean(meanAcrossWaves, 2, 'omitnan');
    sig_dist = std(meanAcrossWaves, 0, 2, 'omitnan');

    if m <= 3
        subplot(2, 3, m);
    else
        if m == 4, subplot(2, 3, 4); else, subplot(2, 3, 5); end
    end

    shadedErrorBar(theta, mu_dist, 2*sig_dist, 'lineProps', {'Color', mColors(m,:), 'LineWidth', 1.8});
    hold on;
    plot(theta, trueDist, 'k--', 'LineWidth', 2.0);
    grid on;
    axis([-pi pi 0 0.085]);
    title(mNames{m}, 'FontWeight', 'bold', 'Color', mColors(m,:), 'FontSize', 13);
    xlabel('\theta [rad]'); ylabel('E(\theta) [m^2/rad]');
    set(gca, 'FontSize', 11);
end
exportgraphics(fig12, fullfile(figDir, 'Figure12.png'), 'Resolution', 200);
exportgraphics(fig12, fullfile(figDir, 'Figure12.eps'), 'ContentType', 'vector');


%% --- FIGURE 13: Estimator Metrics Across N (N=3, 7, 15) ---
fprintf('[13/13] Generating Figure 13 (Estimator Metrics across N=3, 7, 15)...\n');
nCases = [3, 7, 15];
symbols = {'^', 'o', 's'};
peak_prec = zeros(5, 3);
miae_val = zeros(5, 3);

for ni = 1:3
    Nval = nCases(ni);
    dirN = fullfile(dataBase, 'methods', sprintf('N%d', Nval));
    filesN = dir(fullfile(dirN, '*.mat'));
    for m = 1:5
        for f = 1:numel(filesN)
            if startsWith(filesN(f).name, mKeys{m}, 'IgnoreCase', true)
                d = load(fullfile(dirN, filesN(f).name));
                break;
            end
        end
        theta = d.sampledDirectionalSpectrum.theta(:);
        trueDist = d.parametricDistribution(:);
        trueDistNorm = trueDist ./ trapz(theta, trueDist);

        [~, maxIdx] = max(d.sampledDirectionConf, [], 1);
        peakThetas = theta(maxIdx) * (180/pi);
        std_per_geo = squeeze(std(peakThetas, 0, 2));
        peak_prec(m, ni) = mean(std_per_geo, 'omitnan');

        meanAcrossWaves = squeeze(mean(d.sampledDirectionConf, 2, 'omitnan'));
        numCfg = size(meanAcrossWaves, 2);
        iae_cfg = zeros(numCfg, 1);
        for k = 1:numCfg
            dist_k = meanAcrossWaves(:, k);
            if any(isnan(dist_k)), iae_cfg(k) = NaN; continue; end
            dist_k_norm = dist_k ./ trapz(theta, dist_k);
            iae_cfg(k) = trapz(theta, abs(dist_k_norm - trueDistNorm));
        end
        miae_val(m, ni) = mean(iae_cfg, 'omitnan');
    end
end

fig13 = figure(13); clf;
set(fig13, 'Color', 'w', 'Position', [100 100 900 650]);

subplot(2, 1, 1); hold on; grid on; box on;
for m = 1:5
    for ni = 1:3
        plot(m + (ni-2)*0.15, peak_prec(m, ni), symbols{ni}, ...
            'Color', mColors(m,:), 'MarkerFaceColor', mColors(m,:), 'MarkerSize', 8, 'LineWidth', 1.5);
    end
end
set(gca, 'XTick', 1:5, 'XTickLabel', mNames, 'FontSize', 12);
ylabel('\sigma(\theta_p) [^\circ]'); title('Peak Direction Precision Across Sensor Count (N)');
ylim([0 45]);

subplot(2, 1, 2); hold on; grid on; box on;
for m = 1:5
    for ni = 1:3
        plot(m + (ni-2)*0.15, miae_val(m, ni), symbols{ni}, ...
            'Color', mColors(m,:), 'MarkerFaceColor', mColors(m,:), 'MarkerSize', 8, 'LineWidth', 1.5);
    end
end
set(gca, 'XTick', 1:5, 'XTickLabel', mNames, 'FontSize', 12);
ylabel('MIAE [rad^{-1}]'); title('Total Distribution Accuracy Across Sensor Count (N)');
ylim([0 0.55]);

exportgraphics(fig13, fullfile(figDir, 'Figure13.png'), 'Resolution', 200);
exportgraphics(fig13, fullfile(figDir, 'Figure13.eps'), 'ContentType', 'vector');

fprintf('\n=========================================================\n');
fprintf(' All 13 figures generated successfully in: %s\n', figDir);
fprintf('=========================================================\n');

end


%% ===== SUBFUNCTIONS =====

function helper_config_coarray(dataFile, figNum, figTitle, outName, figDir)
    d = load(dataFile);
    if isfield(d, 'lambdap'), Lp = d.lambdap; else, Lp = (9.81 * d.Tp^2) / (2 * pi); end
    sumerror = sum(abs(d.sampledEnsembleDistribution - d.parametricDistribution'));
    [~, I] = sort(sumerror);
    plotRange = [1, 2, numel(I)-1, numel(I)];

    custom_clr = [0.0 0.2 0.6; 0.3 0.5 0.9; 1.0 0.3 0.3; 0.7 0.0 0.0];
    fig = figure(figNum); clf;
    set(fig, 'Color', 'w', 'Position', [100 100 900 450]);

    % Normalized Array Geometry
    subplot(1, 2, 1); hold on; grid on; box on;
    for k = 1:length(plotRange)
        idx = I(plotRange(k));
        xp = d.x(:, idx) / Lp;
        yp = d.y(:, idx) / Lp;
        cx = mean(xp); cy = mean(yp);
        ang = atan2(yp - cy, xp - cx);
        [~, sIdx] = sort(ang);
        pgon = polyshape(xp(sIdx), yp(sIdx));
        plot(pgon, 'FaceColor', custom_clr(k,:), 'FaceAlpha', 0.12, 'EdgeColor', custom_clr(k,:), 'LineWidth', 1.2);
        plot(xp, yp, '.', 'MarkerSize', 18, 'Color', custom_clr(k,:));
    end
    axis equal; axis square;
    xlabel('x / \lambda_p'); ylabel('y / \lambda_p');
    title([figTitle, ' - Geometry']);
    set(gca, 'FontSize', 12);

    % Co-array Coverage
    subplot(1, 2, 2); hold on; grid on; box on;
    for k = 1:length(plotRange)
        idx = I(plotRange(k));
        coarraySubplot(d.x(:, idx), d.y(:, idx), d.Tp, custom_clr(k,:));
    end
    axis square;
    title([figTitle, ' - Co-array Coverage']);
    set(gca, 'FontSize', 12);

    exportgraphics(fig, fullfile(figDir, [outName, '.png']), 'Resolution', 200);
    exportgraphics(fig, fullfile(figDir, [outName, '.eps']), 'ContentType', 'vector');
end

function plot_sensitivity(sensFolder, figNum, xLabelText, levelLabels, outName, figDir)
    files = dir(fullfile(sensFolder, '*.mat'));
    numL = numel(files);
    clr = turbo(numL);

    fig = figure(figNum); clf;
    set(fig, 'Color', 'w', 'Position', [100 100 1000 450]);

    ax1 = subplot(1, 2, 1); hold on; grid on; box on;
    ax2 = subplot(1, 2, 2); hold on; grid on; box on;

    for i = 1:numL
        d = load(fullfile(sensFolder, files(i).name));
        theta = d.parametricDirectionalSpectrum.theta(:);
        trueDist = trapz(d.parametricDirectionalSpectrum.w, d.parametricDirectionalSpectrum.S, 2);
        trueDistNorm = trueDist ./ trapz(theta, trueDist);

        meanAcrossWaves = squeeze(mean(d.sampledDirectionConf, 2, 'omitnan'));
        dist_mu = mean(meanAcrossWaves, 2, 'omitnan');
        dist_mu_norm = dist_mu ./ trapz(theta, dist_mu);

        plot(ax1, theta, dist_mu_norm, 'Color', clr(i,:), 'LineWidth', 1.8, 'DisplayName', levelLabels{i});
        if i == numL
            plot(ax1, theta, trueDistNorm, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Target');
        end

        % Error metrics
        [~, maxIdx] = max(d.sampledDirectionConf, [], 1);
        peakThetas = theta(maxIdx) * (180/pi);
        std_per_cfg = squeeze(std(peakThetas, 0, 1));
        p_mean = mean(std_per_cfg, 'omitnan');
        p_std = std(std_per_cfg, 'omitnan');

        numCfg = size(meanAcrossWaves, 2);
        iae_cfg = zeros(numCfg, 1);
        for k = 1:numCfg
            dk = meanAcrossWaves(:, k);
            if any(isnan(dk)), iae_cfg(k) = NaN; continue; end
            iae_cfg(k) = trapz(theta, abs((dk ./ trapz(theta, dk)) - trueDistNorm));
        end
        i_mean = mean(iae_cfg, 'omitnan');
        i_std = std(iae_cfg, 'omitnan');

        errorbar(ax2, i, p_mean, p_std, 'o', 'Color', clr(i,:), 'MarkerFaceColor', clr(i,:), 'LineWidth', 1.8, 'CapSize', 6);
    end

    xlabel(ax1, '\theta [rad]'); ylabel(ax1, 'E(\theta) [rad^{-1}]');
    title(ax1, '(a) Directional Distribution');
    legend(ax1, 'Location', 'northeast', 'FontSize', 9);
    set(ax1, 'FontSize', 11);

    ylabel(ax2, '\sigma(\theta_p) [^\circ]');
    title(ax2, ['(b) Precision vs ', xLabelText]);
    set(ax2, 'XTick', 1:numL, 'XTickLabel', levelLabels, 'FontSize', 10);

    exportgraphics(fig, fullfile(figDir, [outName, '.png']), 'Resolution', 200);
    exportgraphics(fig, fullfile(figDir, [outName, '.eps']), 'ContentType', 'vector');
end
