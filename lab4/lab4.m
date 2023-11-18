% dualne: min (t(b)y_1 + t(g)y_2)
% t(A)y_1 + Iy_2 - Iy_3 = c
% y >= 0

% Bazę początkową robimy tak: Jak c jest dodatni, to weź y_2, a jak c jest
% ujemny, to weź y_3

% Uwaga: jak się zmienia znak w wierszu, to trzeba pamietać, żeby porem
% odwrócić z powrotem. Wspomnieć o tym w raporcie, żeby pokazać Ewci, że
% rozumiemy.

% NOTE: Do raportu można mniejsze n i m.

%% Test
n = 5; % Długość wektora c i g
m = 10; % Liczba wierszy macierzy A i długość wektora b

[c, A, b, g] = drawData(n, m)

[x,fval,exitflag,output,lambda] = linprog(c, A, b, [], [], zeros(1, n), g)


%% Definicje użytych funkcji

% Funkcja drawData generuje dane do testów.
function [c, A, b, g] = drawData(n, m)
    % Generowanie wektora c
    c = randi([-5, 5], 1, n);
    
    % Generowanie macierzy A
    A = randi([-5, 5], m, n);
    
    % Generowanie wektora b
    b = randi([-5, 5], 1, m);
    
    % Generowanie wektora g
    g = randi([1, 30], 1, n);
end
