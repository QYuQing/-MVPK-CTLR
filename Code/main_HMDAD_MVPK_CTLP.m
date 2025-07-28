clear
seed = 1;
rng(seed, 'twister');
nfolds = 5;
MU = 0.005;
GAMMA = 0.002;
nruns = 20;
% dataname = 'HMDAD';
dataname = 'Disbiome';
load(['F:\KronRLS\github\Dataset\HMDA\all_data_' dataname '.mat']);


crossval_idx = crossvalind('Kfold', length(y(:)), nfolds);

runs_aupr_list = [];
runs_auc_list = [];
for run = 1:nruns
    %% split folds
    crossval_idx = crossvalind('Kfold', length(y(:)), nfolds);
    
    fold_aupr_kronls_ka = [];
    fold_auc_kronls_ka = [];
    fold_times = [];
    for fold = 1:nfolds
        %% mask
        train_idx = find(crossval_idx ~= fold);
        test_idx  = find(crossval_idx == fold);
        
        y_train = y;
        y_train(test_idx) = 0;
        
        %% wknkn
        eta_v=1;k_nn=17;
        [F_new] = wknkn(y_train,Knormalized(K1_fun),Knormalized(K2_fun),k_nn,eta_v);
        %% kernel
        gamma = 16;
        K1(:,:,1) = kernel_cosine(F_new,1,MU,GAMMA);
        K1(:,:,2) = kernel_corr(F_new,1,MU,GAMMA);
        K1(:,:,3) = kernel_gip(F_new,1, gamma);
        K1(:,:,4) = K1_fun;
        
        K2(:,:,1) = kernel_cosine(F_new,2,MU,GAMMA);
        K2(:,:,2) = kernel_corr(F_new,2,MU,GAMMA);
        K2(:,:,3) = kernel_gip(F_new,2, gamma);
        K2(:,:,4) = K2_fun;
        %% model
        t1=clock;
        lambda = 0.03125;
        mu = 0.25;
        rho =0.015625;
        sigma = 0.015625;
        eta = 2^11;
        iteration = 1;
        F_index = y_train;
        
        [beta,omega,theta,F_pre_list] = MVPK_CTLP(F_new,F_index,K1,K2,lambda,mu,eta,rho,sigma,iteration);
        kronrls_A_KA = combine_pre(F_pre_list, beta);
        t2=clock;
        fold_time = etime(t2,t1);
        fold_times = [fold_times,fold_time];
        %% evaluate predictions
        yy = y;
        [~,~,~,aupr_kronrls_A_KA] = perfcurve(yy(test_idx), kronrls_A_KA(test_idx), 1, 'xCrit', 'reca', 'yCrit', 'prec');
        [~,~,~,AUC_kronrls_KA] = perfcurve(yy(test_idx), kronrls_A_KA(test_idx), 1);
        
        fold_aupr_kronls_ka = [fold_aupr_kronls_ka; aupr_kronrls_A_KA];
        fold_auc_kronls_ka = [fold_auc_kronls_ka; AUC_kronrls_KA];
        mean_aupr = mean(fold_aupr_kronls_ka);
        mean_auc = mean(fold_auc_kronls_ka);
        runs_aupr_list = [runs_aupr_list,mean_aupr];
        runs_auc_list = [runs_auc_list,mean_auc];
        fprintf(' mean AUPR=%.4f, mean AUC=%.4f\n', mean_aupr, mean_auc);
        fprintf(' mean time = %.4f\n', mean(fold_times));
    end
    fprintf(' mean AUPR=%.4f std=%.4f, mean AUC=%.4f std=%.4f\n', mean(runs_aupr_list),std(runs_aupr_list), mean(runs_auc_list),std(runs_auc_list));
end