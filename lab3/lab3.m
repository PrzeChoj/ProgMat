%% Zadanko
% Jak uzyskać postać zadania liniowego?

n = 5;
%x = [1,2,3,4,5];
%y = [100, 101, 99, 98, 102];

x=(1:1:5);
y=[1450, 610, 380, 70, 25];

A = zeros(2*n, 3);
b = zeros(2*n, 1);
for i = 1:n
    A(2*(i-1)+1,:) = [x(i), 1, -1];
    b(2*(i-1)+1,:) = y(i);

    A(2*(i-1)+2,:) = [-x(i), -1, -1];
    b(2*(i-1)+2,:) = -y(i);
end

f = [0, 0, 1];

[x_odp,fval,exitflag,output,lambda] = linprog(f,A,b,[],[],[-inf -inf 0],[inf inf inf])
%fval = 241



a = x_odp(1);
b = x_odp(2);

plot(x, y, 'o');
hold on; % Zatrzymaj aktualny wykres, aby dodać kolejne elementy

% Narysuj linię ax + b
linia_x = min(x):0.1:max(x); % Zakres x dla linii
linia_y = a * linia_x + b; % Wartości y dla linii
plot(linia_x, linia_y, '-', 'LineWidth', 2); % '--' oznacza przerywaną linię

xlabel('Oś X'); % Dodaj etykietę dla osi X
ylabel('Oś Y'); % Dodaj etykietę dla osi Y
title('Wykres x-y'); % Dodaj tytuł wykresu
grid on; % Włącz siatkę

% Opcjonalnie możesz dodać legendę, jeśli masz więcej niż jeden zestaw danych
legend('y(x)');



