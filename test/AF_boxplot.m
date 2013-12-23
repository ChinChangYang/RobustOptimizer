function AF_boxplot
load('err_con_minmax_minmaxtcjadebin_201311140003.mat')
X = reshape(err(1, :, :), 16, 1);
Y = reshape(err(2, :, :), 16, 1);
load('err_con_minmax_minmaxtcjadebin_201311140048.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311140132.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311140217.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311140301.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311141006.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311141047.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311141127.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311141208.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311141248.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311141402.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311141441.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311141533.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311141652.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
load('err_con_minmax_minmaxtcjadebin_201311141737.mat')
X = [X, reshape(err(1, :, :), 16, 1)];
Y = [Y, reshape(err(2, :, :), 16, 1)];
X = log10(X);
Y = log10(Y);
subplot(211);
boxplot(X);
subplot(212);
boxplot(Y);
end

