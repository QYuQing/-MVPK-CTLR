clear

dataname = 'nr';
load([' dataname '.mat']);
disp(dataname)

%% para
fold = 1;
maxIter_list = [1, 5, 10, 15, 20, 25, 30, 35, 40];
lambda = 2^0;
tol = 1e-6;
gamma = 0.5;
num_bins = 30;
% bin_limits = [-0.5, 1.2];
bin_limits = [-0.6, 1.2];

aupr_list = zeros(length(maxIter_list), 1);
auc_list = zeros(length(maxIter_list), 1);
time_list = zeros(length(maxIter_list), 1);
score_all_list = cell(length(maxIter_list), 1);
R_all_list = cell(length(maxIter_list), 1);
obj_list = cell(length(maxIter_list), 1);

seed = 2026;
rand('seed', seed);

crossval_idx = crossvalind('Kfold', y(:), 5);
y_vec_true = y(:);
N_total = length(y_vec_true);

test_idx = find(crossval_idx == fold);
train_idx = find(crossval_idx ~= fold);

test_true_idx = test_idx(y_vec_true(test_idx) == 1);
fake_negative_idx = test_true_idx;
real_positive_idx = train_idx(y_vec_true(train_idx) == 1);

for iter_idx = 1:length(maxIter_list)
    current_maxIter = maxIter_list(iter_idx);
    t1 = clock;

    fprintf('Running fold = %d, maxIter = %d\n', fold, current_maxIter);

    %% split
    y_train = y;
    y_train(test_idx) = 0;

    %% kernel
    K1 = kernel_gip(y_train, 1, gamma);
    K2 = kernel_gip(y_train, 2, gamma);
% K1 = K1_fun;
% K2 = K2_fun;
    %% Kronecker kernel
    % K1: nd ¡Á nd
    % K2: nt ¡Á nt
    K = kron(K2, K1);

    %% flatten label
    y_vec_train = y_train(:);
    N = length(y_vec_train);

    %% one-hot label
    Y_train_onehot = zeros(N, 2);
    Y_train_onehot(:, 1) = (y_vec_train == 0);
    Y_train_onehot(:, 2) = (y_vec_train == 1);

    %% train RLS with label relearning
    [Y_pre, A, R, obj] = RLS_label_relearning(Y_train_onehot, K, lambda, current_maxIter, tol);

    %% collect prediction scores and relearned labels
    score_all = Y_pre(:, 2);
    R_all = R;

    %% evaluate on test entries only
    score_test = score_all(test_idx);
    y_test = y_vec_true(test_idx);
    [~, ~, ~, aupr_kronls_ka] = perfcurve(y_test, score_test, 1, ...
        'xCrit', 'reca', 'yCrit', 'prec');
    [~, ~, ~, auc_kronls_ka] = perfcurve(y_test, score_test, 1);

    t2 = clock;
    time_list(iter_idx) = etime(t2, t1);
    aupr_list(iter_idx) = aupr_kronls_ka;
    auc_list(iter_idx) = auc_kronls_ka;
    score_all_list{iter_idx} = score_all;
    R_all_list{iter_idx} = R_all;
    obj_list{iter_idx} = obj;

    fprintf('Finish fold = %d, maxIter = %d: AUPR = %f, AUC = %f\n', ...
        fold, current_maxIter, aupr_kronls_ka, auc_kronls_ka);

    plot_bin_figure(score_all(real_positive_idx), score_all(fake_negative_idx), ...
        R_all(real_positive_idx, 2), R_all(fake_negative_idx, 2), ...
        num_bins, bin_limits, sprintf('Iteration %d', current_maxIter), ...
        sprintf('bin_%s_fold%d_iter%d.pdf', dataname, fold, current_maxIter));
end

save(['relearning_analysis_' dataname '_fold' num2str(fold) '.mat'], ...
    'dataname', 'fold', 'maxIter_list', 'lambda', 'tol', 'gamma', 'num_bins', 'bin_limits', ...
    'crossval_idx', 'test_idx', 'train_idx', 'test_true_idx', 'fake_negative_idx', ...
    'real_positive_idx', 'y', 'y_vec_true', 'aupr_list', 'auc_list', 'time_list', ...
    'score_all_list', 'R_all_list', 'obj_list');

function plot_bin_figure(real_positive_score, fake_positive_score, real_positive_r, fake_positive_r, ...
    num_bins, bin_limits, fig_title, save_name)

    fig = figure('Visible', 'off');
    hold on;

    bin_edges = linspace(bin_limits(1), bin_limits(2), num_bins + 1);

    histogram(real_positive_score, bin_edges, ...
        'DisplayName', 'Predicted Score (Known Positives)', ...
        'FaceColor', [221/256,102/256,112/256], ...
        'FaceAlpha', 0.40, 'EdgeColor', 'none');

    histogram(fake_positive_score, bin_edges, ...
        'DisplayName', 'Predicted Score (Unknown Positives)', ...
        'FaceColor', [12/256,53/256,71/256], ...
        'FaceAlpha', 0.8, 'EdgeColor', 'none');

    histogram(real_positive_r, bin_edges, ...
        'DisplayName', 'Relearned Label (Known Positives)', ...
        'FaceColor', [237/256,174/256,147/256], ...
        'FaceAlpha', 0.40, 'EdgeColor', 'none');

    histogram(fake_positive_r, bin_edges, ...
        'DisplayName', 'Relearned Label (Unknown Positives)', ...
        'FaceColor', [89/256,143/256,145/256], ...
        'FaceAlpha', 0.8, 'EdgeColor', 'none');

    xlim(bin_limits);
    xlabel('Value');
    ylabel('Count');
    title(fig_title);
    legend('show');
    grid on;
    hold off;

    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'PaperPosition', [0 0 8 7]);
    print(fig, save_name, '-dpdf', '-painters');
    close(fig);
end
