clear
seed = 1;
rand('seed', seed);
nfolds = 10; nruns=10;
MU = 0.005; GAMMA = 0.002;

dataname = 'gpcr';
load(['F:\KronRLS\github\Dataset\DTI\all_data_' dataname '.mat']);
dataname


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
        %% downsample
        num_zeros = sum(y(test_idx) == 0);
        num_ones = sum(y(test_idx) == 1);
        y_test = y(test_idx);
        ones_index = find(y_test==1);
        zero_index = find(y_test==0);
        ones_index = ones_index(1:num_ones);
        zero_index = zero_index(1:1*num_ones);
        test_idx = [test_idx(ones_index);test_idx(zero_index)];
        
        y_train = y;
        y_train(test_idx) = 0;
        
        %% load kernels
        gamma = 0.5;
        K1(:,:,1) = K1_1;
        K1(:,:,2) = kernel_gip(y_train,1,gamma);
        K1(:,:,3) = Knormalized(kernel_cosine(y_train,1,MU,GAMMA));
        K1(:,:,4) = Knormalized(kernel_corr(y_train,1,MU,GAMMA));
        
        K2(:,:,1) = K2_1;
        K2(:,:,2) = K2_2;
        gamma_fp = 4;
        K2(:,:,3) = Knormalized(kernel_gip(Drug_MACCS_fingerprint,1,gamma_fp));
        K2(:,:,4) = Knormalized(kernel_gip(y_train,2,gamma));
        K2(:,:,5) = Knormalized(kernel_cosine(y_train,2,MU,GAMMA));
        K2(:,:,6) = Knormalized(kernel_corr(y_train,2,MU,GAMMA));
        
        %% wknkn
        [~,~,len_k1] = size(K1);
        [~,~,len_k2] = size(K2);
        K1_com = combine_pre(K1, ones(1,len_k1)/len_k1);
        K2_com = combine_pre(K2, ones(1,len_k2)/len_k2);
        eta_v = 0.9; k_nn = 20;
        F_new = wknkn(y_train, Knormalized(K1_com), Knormalized(K2_com), k_nn, eta_v);
        t1=clock;
        
        %% model
        lambda = 2^-5;
        mu = 2^-5;
        rho = 2^-6;
        sigma = 2^0;
        eta = 2^7;
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
    end
    
    mean_aupr = mean(fold_aupr_kronls_ka);
    mean_auc = mean(fold_auc_kronls_ka);
    runs_aupr_list = [runs_aupr_list,mean_aupr];
    runs_auc_list = [runs_auc_list,mean_auc];
    fprintf(' mean AUPR=%.4f, mean AUC=%.4f\n', mean_aupr, mean_auc);
    fprintf(' mean time = %.4f\n', mean(fold_times));
end
fprintf(' mean AUPR=%.4f std=%.4f, mean AUC=%.4f std=%.4f\n', mean(runs_aupr_list),std(runs_aupr_list), mean(runs_auc_list),std(runs_auc_list));