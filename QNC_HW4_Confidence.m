mu = 10;
sigma = 2; 
 
interval = 0.95;

%simulate the data 

fprintf('Method 1\n');
for N = [5 10 20 40 80 160 1000]
    rng default
    samples = normrnd(mu,sigma,N,1); 

    %Method 1 

    %Calculate sem
    stderror = std(samples)/sqrt(N); 

    %For a 95% confidence interval:
    z_score = -norminv(0.5*(1-interval));
    
    fprintf('N = %d samples :  The %d percent confidence interval = [%2.2f - %2.2f , %2.2f + %2.2f]\n',N, interval*100,mean(samples),z_score*stderror,mean(samples),z_score*stderror);
    
    
end

%% Method 2
 fprintf('Method 2\n');
for N = [5 10 20 40 80 160 1000]
    rng default
    samples = normrnd(mu,sigma,N,1);
    %define degrees of freedom
    nu = N-1; 

    %Find the upper and lower confidence bounds for the 95% confidence interval
    alpha = 1-interval;
    pLo = alpha/2;
    pUp = 1 - alpha/2; 

    %Compute the critical values for the confidence bounds

    crit = tinv([pLo pUp],nu);
   
    fprintf('N = %d samples : The %d percent confidence interval = [%2.2f  %2.2f , %2.2f + %2.2f]\n',N, interval*100,mean(samples),crit(1)*stderror,mean(samples),crit(2)*stderror);
end

%% Method 3
 fprintf('Method 3\n');
 for N = [5 10 20 40 80 160 1000]
    rng default
    samples = normrnd(mu,sigma,N,1);
    %Take 1000 bootstrapped samples from the data and compute the mean
    resamples = 1000;
    m = bootstrp(resamples,@mean,samples); 
    
    delta = mean(samples) - m; 
    delta = sort(delta); 

    %For a 95 %confidence interval, take de 2.5 and the 95 percentile of
    %the differences

    lower = delta(0.025 *  resamples);
    higher = delta(0.95 * resamples); 

   
    fprintf('N = %d samples: The %d percent confidence interval = [%2.2f  %2.2f , %2.2f + %2.2f]\n',N, interval*100,mean(samples),lower,mean(samples),higher);
 end
 %% Method 4
 fprintf('Method 4\n');
for N = [5 10 20 40 80 160 1000]
    rng default
    samples = normrnd(mu,sigma,N,1);
    %As discussed in the noteboook, the bayesian credible intervals is
    %calculated as in Method 1 

    %Calculate sem
    stderror = std(samples)/sqrt(N); 

    %For a 95% confidence interval:
    z_score = -norminv(0.5*(1-interval));

   
    fprintf('N = %d samples: The %d percent confidence interval = [%2.2f - %2.2f , %2.2f + %2.2f]\n',N, interval*100,mean(samples),z_score*stderror,mean(samples),z_score*stderror);

end
