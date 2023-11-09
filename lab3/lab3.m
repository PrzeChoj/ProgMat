%% Tests
n = 10;
[x,y] = drawData(n)
[A,b,f,Aeq,beq,lb,ub] = getAandb(x,y)
[x_odp,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub)

a = x_odp(1);
b = x_odp(2);
%abs(max(abs(a*x + b - y)) - fval) < 0.00000001 % powinna byc rownosc
plot_solution(x, y, a, b)

function [x, y] = drawData(n)
    x=(1:1:n)';
    y=randi([0 1500],n,1);
end

function [A,b,f,Aeq,beq,lb,ub] = getAandb(x,y)
    n = length(x);
    A = zeros(2*n, 3);
    b = zeros(2*n, 1);

    for i = 1:n
        A(2*(i-1)+1,:) = [x(i), 1, -1];
        b(2*(i-1)+1,:) = y(i);
    
        A(2*(i-1)+2,:) = [-x(i), -1, -1];
        b(2*(i-1)+2,:) = -y(i);
    end
    
    f = [0,0,1];

    Aeq = []; beq = [];
    lb = [-inf -inf 0]; ub = [inf inf inf];
end

function [out] = singleTest(n)
    [x,y] = drawData(n);
    [A,b,f,Aeq,beq,lb,ub] = getAandb(x,y);

    [x_odp,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub);
    % TODO
end

function [] = plot_solution(x, y, a, b)
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
end
