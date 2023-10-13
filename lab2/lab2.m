%% Zad 1
% Tak samo jak Zad 2, tylko -f, bo chcemy max nie min

%% Zad 2
A = [-2 4 -3 -4 -2
    10 10 3 4 2
    -9 -6 5 -1 -5];

b = [2 1 -2];

Ar = [1 0 1 0 1];
br = [10];

f = [-10 -8 2 9 -3];

my_options = optimoptions('linprog','Algorithm','dual-simplex');
%my_options = optimoptions('linprog','Algorithm','interior-point');
[x,fval,exitflag,output,lambda] = linprog(f,A,b,Ar,br,[-inf 0 0 -inf -inf],[0 inf inf inf inf],my_options)

%% Zad3
%    |s1|s2|s3|
% ---|--|--|--|
% 80 | 2| 1| 0| >= 3*N
% 30 | 1| 4| 6| >= 2*N
% odp|10| 0|20|

N = 7;
A = [2 1 0
    1 4 6];
b = [3*N 2*N];
f = [1 1 1];

[x,fval,exitflag,output,lambda] = linprog(f,-A,-b,[],[],[0 0 0])
% Akurat się zrobiło całkowite. Nie musiało tak być...

%% Zad 4
A = [1 0 0 1 1 1 1 % ci co będą pracować w poniedziałek
    1 1 0 0 1 1 1
    1 1 1 0 0 1 1
    1 1 1 1 0 0 1
    1 1 1 1 1 0 0
    0 1 1 1 1 1 0
    0 0 1 1 1 1 1];
b = [15 10 15 20 15 16 10];
f = [1 1 1 1 1 1 1];
my_options = optimoptions('linprog','Algorithm','dual-simplex');
%my_options = optimoptions('linprog','Algorithm','interior-point');
[x,fval,exitflag,output,lambda] = linprog(f,-A,-b,[],[],[0 0 0 0 0 0 0],[],my_options)

