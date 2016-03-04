function [K_genes,n_genes] = Best_Hill(T,AG)

    global x_axis y_axis

    x_axis          = T;        % X-axis reference for the best Hill function computation
    y_axis          = AG;       % Y-axis reference for the best Hill function computation

    MYOPT = optimset('fminsearch');
    MYOPT = optimset(MYOPT,'MaxFunEvals',20000,'MaxIter',300, 'TolFun',1.e-8,'TolX',1.e-8, 'Diagnostics','off','Display','off','LargeScale','off','HessUpdate','bfgs');

    mu_0 = [10;4];
    [mu, ~] = fminsearch('Loss_Hill_genes',mu_0,MYOPT);

    K_genes	= mu(1);      % Threshold for the best Hill function
    n_genes = mu(2);      % Hill coeff. for the best Hill function

end