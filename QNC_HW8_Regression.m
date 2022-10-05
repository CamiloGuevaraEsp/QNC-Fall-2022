
close all
alpha = 0.05;
age = [3 4 5 6 7 8 9 11 12 14 15 16 17];
age = age';
wing_length = [1.4 1.5 2.2 2.4 3.1 3.2 3.2 3.9 4.1 4.7 4.5 5.2 5.0];
wing_length  =wing_length';
%% Exercise 1
scatter(age,wing_length)
xlabel('Age')
ylabel('Wing Length (arbitrary units)')
hold on
%% Exercise 2
%Calculate slope
b1 = age\wing_length;
%Calculate y-intercept
a = mean(wing_length) - (b1 * mean(age)); 
yCalc = (b1*age) + a; 
plot(age,yCalc);
%% Exercise 3
mdl = fitlm(age, wing_length);
[p,F] = coefTest(mdl);
fprintf('Exercise 3 \n')
if p < alpha
    fprintf('p-value (%e) < %2.3f , Ho is rejected\n',p,alpha);
else
    fprintf('p-value (%ef) < %2.3f , Ho is not rejected\n',p,alpha);
end
%% Exercise 4
figure
plot(mdl); 
xlabel('Age')
ylabel('Wing Length (arbitrary units')
title('Plot using Matlab custom function');

x = tinv(alpha/2,length(age)-2);
up = b1 + (x * sqrt(mdl.MSE/sum((age - mean(age)).^2)));
low = b1 - (x * sqrt(mdl.MSE/sum((age - mean(age)).^2)));
fprintf('Exercise 4 \n')
fprintf('The 95 percent confidence intervals are %2.3f ; %2.3f\n', low, up);

figure
scatter(age,wing_length)
xlabel('Age')
ylabel('Wing Length (arbitrary units)')
hold on
plot(age, (up*age) + a , '--')
plot(age,(low*age) + a , '--')
plot(age,yCalc);
title('Plot using calculated confidence intervals')

%% Exercise 5
r_squared = mdl.Rsquared.Ordinary;
fprintf('Exercise 5 \n')
fprintf('R^2 = %2.3f\n',r_squared);
%% Exercise 6
x = age;
y = wing_length;
rxy = sum((x - mean(x)) .* (y-mean(y))) /  (sqrt(sum((x-mean(x)).^2)) * sqrt(sum((y-mean(y)).^2)));
rxy_b = corrcoef(age,wing_length);
fprintf('Exercise 6 \n')
fprintf('r(calculated) = %2.3f\n',rxy);
fprintf('r(Matlab) = %2.3f\n',rxy_b(2,1));
%% Exercise 7 
%%Adding noise
%Add noise from a normal distribution
noise = normrnd(0, 1, size(wing_length));
wing_noise = wing_length + noise;

noisy = fitlm(age, wing_noise);
figure
subplot(1,2,1)
plot(mdl); 
xlabel('Age')
ylabel('Wing Length (arbitrary units')
title('Original data');

subplot(1,2,2)
plot(noisy); 
xlabel('Age')
ylabel('Wing Length (arbitrary units')
title('Data with noise');

r_squared_noisy = noisy.Rsquared.Ordinary;
r = corrcoef(age,wing_noise);
fprintf('Exercise 7 \n')
fprintf('R^2(original) = %2.3f\n',r_squared);
fprintf('R^2(noisy) = %2.3f\n',r_squared_noisy);
fprintf('r(Original) = %2.3f\n',rxy);
fprintf('r(noisy) = %2.3f\n',r(2,1));
