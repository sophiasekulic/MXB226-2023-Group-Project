%% Ambient Temperature

% Point (0.03, 0.03) corresponds to node 22
x_full(22) % temperature at node 22 with ambient temp of 20 degrees.
%%
syms a
B2 = zeros(33,1);
B2 = symmatrix(B2);
B2(6) = 40;
B2(12) = 40;
B2(18) = 40;
B2(24) = 2.*a.*sqrt(2)/30;
B2(29) = 2.*a.*sqrt(2)/30;
B2(30) = 70;
B2(31) = 70;
B2(32) = 70;
B2(33) = 2.*a.*sqrt(2)/30;
B2 = symmatrix2sym(B2);
%%
sympref('FloatingPointOutput',true);
xsym = A \ B2
%%
node22 = xsym(22)
lowtemp = node22 == 50
lower = solve(lowtemp, a)
hightemp = node22 == 55
higher = solve(hightemp, a)
%% 
% Therefore, an ambient temperature range of 1.7274 to 55 will allow the temperature 
% at point (0.03, 0.03) to remain between 50 and 55.

% Check by calculating the temp at node 22 with ambient temps ranging from
% 0 to 60 degrees

temp = 0:1:60;
length(temp);
comp = zeros(1,length(temp));

for i = 0:length(temp)-1
    B1 = B2;
    B1 = subs(B1, a, i);
    x = A \ B1;
    comp(i+1) = x(22);
end

temp, comp