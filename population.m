% Read in data
data = readmatrix("population.csv");
year = data(:,1);
pop = data(:,2);

%% Exponential Fit
% We want to fit to pop=ae^(bx) where a and b are unknown.
% To do so, we can make the equation linear and get
% pop_tilda = a_tilda + bx where a_tilda=log(a) and pop_tilda=log(pop).
V = [ones(size(year)) year];% Data Matrix
pop_tilda = log(pop);

% To find a and b, we'll solve Vc=pop_tilda with QR decomp
% which automatically finds the least square solution.
% Note: c = [a_tilda, b]
[Q,R]=qr(V,0);
c = backsolve(R,Q'*pop_tilda);
% Get Parameters
a = exp(c(1));
b = c(2);
              
% Plot
figure(1)
plot(year,pop,'o')
hold on
plot(year, a.*exp(b.*year))
xlabel('Year')
ylabel('Population')
legend('Data','Fit')
title('Exponential Fit with QR Decomposition')
hold off

%% Power-Law Fit
% We need to scale both axes to (0 1] and
% fit to power law pop = a*year^b
% pop_scaled = a*year_adj_scaled^b
% pop/max(pop) = a*( (year+abs(min(year)))/max( (year+abs(min(year))) ) )^b
pop_scaled = pop/max(pop);
year_adj = year + abs(min(year)) + 1; %added +1 to not have any zeros
year_adj_scaled = year_adj/max(year_adj);
                                 
% Make power law linear
% pop_scaled_tilda = a_tilda + b*year_adj_scaled_tilda
% where pop_scaled_tilda = log(pop_scaled)
% a_tilda = log(a)
% year_tilda = log(year_adj_scaled)
pop_scaled_tilda = log(pop_scaled);
year_adj_scaled_tilda = log(year_adj_scaled);
V = [ones(size(year_adj_scaled_tilda)) year_adj_scaled_tilda];
                                 
% lets solve Vc=pop_tilda where c = [a_tilda b] with QR Decomp.
% which finds the least square solution
[Q,R]=qr(V,0);
c = backsolve(R,Q'*pop_scaled_tilda);
              
% Get parameters
a = exp(c(1));
b = c(2);
figure(2)
plot(year_adj_scaled,pop_scaled,'o')
hold on
plot(year_adj_scaled, a.*year_adj_scaled.^b)
xlabel('Scaled Year')
ylabel('Scaled Population')
legend('Data','Fit')
title('Scaled and Shifted Power Law Fit with QR Decomposition')
hold off

%% Piecewise Fit - Exponential and Linear
% I will now do linear fit after 1950 and
% an exponential fit before 1950.
% Fit After 1950
year_after = data(1:end-24,1);
pop_after = data(1:end-24,2);
              
% We'd like to fit the data to pop = m*year+b where m and b are unknown.
% Thus the data matrix is the following:
V = [ones(size(year_after)) year_after];
              
% We will solve V*a=population where a = [b;m] using QR Decomposition
[Q,R]=qr(V,0);
c = backsolve(R,Q'*pop_after);
              
% Get Parameters
b = c(1);
m = c(2);
              
% Fit Before 1950
year_before = data(end-23:end,1);
pop_before = data(end-23:end,2);
              
% We want to fit to pop=ae^(dx) where a and d are unknown.
% To do so, we can make the equation linear and get
% pop_tilda = a_tilda + dx where a_tilda=log(a) and pop_tilda=log(pop).
V = [ones(size(year_before)) year_before];% Data Matrix
pop_tilda = log(pop_before);
              
% To find a and d, we'll solve Vc=pop_tilda with QR decomp
% which automatically finds the least square solution.
% Note: c = [a_tilda, d]
[Q,R]=qr(V,0);
c = backsolve(R,Q'*pop_tilda);
              
% Get Parameters
a = exp(c(1));
d = c(2);
              
% Plot
figure(3)
plot(year_after,pop_after,'o','Color','blue')
hold on
plot(year_before,pop_before,'o','Color','green')
plot(year_after,m*year_after+b,'Color','black')
plot(year_before, a.*exp(d.*year_before),'Color','red')
xlabel('Year')
ylabel('Population')
legend({'Data After 1950','Data Before 1950','Linear Fit After 1950','Exponential Fit Before
1950'},'Location','northwest')
title('Exponential and Linear Fit with QR Decomposition')
hold off
