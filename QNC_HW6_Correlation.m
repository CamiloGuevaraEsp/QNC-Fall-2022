clear all
%% Exercise 1
%Define the variables
wing_length = [10.4 10.8 11.1 10.2 10.3 10.2 10.7 10.5 10.8 11.2 10.6 11.4];
tail_length = [7.4 7.6 7.9 7.2 7.4 7.1 7.4 7.2 7.8 7.7 7.8 8.3];

%Rename variables as x and y
x = wing_length;
y = tail_length;

plot(wing_length,tail_length, 'o'); 
xlabel('Wing Length (cm)')
ylabel('Tail Length (cm)')

% Both variables seem positively correlated

%% Exercise 2
rxy = sum((x - mean(x)) .* (y-mean(y))) /  (sqrt(sum((x-mean(x)).^2)) * sqrt(sum((y-mean(y)).^2)));
rxy_b = corrcoef(wing_length,tail_length); 

ryx = sum((y - mean(y)) .* (x-mean(x))) /  (sqrt(sum((y-mean(y)).^2)) * sqrt(sum((x-mean(x)).^2)));
ryx_b = corrcoef(y,x); 

%Both ways to calculate R yield the same result 

%% Exercise 3
%Standard error

sr = sqrt((1-rxy^2)/(length(x)-2)); 

%95 % confidence intervals 

%Calculate zr

zr = 0.5 * log((1 + rxy) / (1- rxy)); 
sz = sqrt(1/(length(x)-3));
zcriterion = norminv(0.025);
zlo = zr + zcriterion * sz;
zup = zr - zcriterion * sz;

rlo = (exp(2*zlo) - 1)/ (exp(2*zlo) + 1); 
rup = (exp(2*zup) - 1)/ (exp(2*zup) + 1);
%% Exercise 4
%compute the statistic
t = rxy/sr; 

%Compute the p-value
p = 2.*(1-tcdf(t, length(x) -2));

fprintf(' p-value = %2.6f\n', p); 
% The null hypothesis is rejected, hence, the r value is different to zero
%% Exercise 5
%Compare with Yales's r
yales_r = 0.75; 
zrs = 0.5 * log((1 + yales_r) / (1- yales_r));

%Calculate the statistic 

lambda = (zr - zrs)/sqrt(1/(length(x) -3 ));

%Compute the p-value
p = 2.*(1-tcdf(lambda, length(x) -2));

fprintf(' p-value = %2.6f\n', p); 

% The null hypothesis is not rejected, there is no difference between my r
% and Yale's R
%% Exercise 6

%Calculate the sample size to detect a difference when r >= 0.5
test_r = 0.5;
%test a range of sample sizes
sample_size = 1:50;
%Calculate t statistic for that value
for i = 1:length(sample_size)
    t = test_r/sqrt((1-test_r^2)/(sample_size(i)-2));
    p(i) = 2.*(1-tcdf(t, sample_size(i) -2));
end

fprintf('Minimun sample size to  reject H0: r = 0 when r >= 0.5 is %d\n',min(find(p <= 0.05)))


