function[cost] = costfun(p, x, y)
    norm_fun = @(mu, sd, x) exp(-(x-mu).^2 / (2*sd^2));

    Ypred = norm_fun(p(1), p(2), x);
    
    cost = sum((y(:)-Ypred(:)).^2);