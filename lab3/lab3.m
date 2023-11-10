%% Tests 2
out = multipleTests(10, 100)

%% Test 3
%x = [1;2;3]
%y = [10;7;5]
%[A_0,b_0] = getAandb(x,y)
n = 5;
[x,y] = drawData(n);
[A,b,f,Aeq,beq,lb,ub] = getAandb(x,y);
x_odp = mySymplex(A,b);
plot_solution(x, y, x_odp)

x_odp = linprog(-f,A,b,Aeq,beq,lb,ub);
plot_solution(x, y, x_odp)

%% Tests
%n = 10;
[x,y] = drawData(n)
%x = [1;2;3]
%y = [10;7;5]
[A,b,f,Aeq,beq,lb,ub] = getAandb(x,y)

[x_odp,fval,exitflag,output,lambda] = linprog(-f,A,b,Aeq,beq,lb,ub)

e = max(y) + fval;
%abs(max(abs(x_odp(1)*x + x_odp(2) - y)) - e) < 0.00000001 % powinna byc rownosc
plot_solution(x, y, x_odp)

%% Functions definitons
function [x, y] = drawData(n)
    x=(1:1:n)';
    y=randi([0 1500],n,1);
end

function [A,b,f,Aeq,beq,lb,ub] = getAandb(x,y)
    n = length(x);
    A = zeros(2*n, 5);
    b = zeros(2*n, 1);

    y_max = max(y);

    for i = 1:n
        A(2*(i-1)+1,:) = [x(i), -x(i), 1, -1, 1];
        b(2*(i-1)+1,:) = y(i) + y_max;
    
        A(2*(i-1)+2,:) = [-x(i), x(i), -1, 1, 1];
        b(2*(i-1)+2,:) = -y(i) + y_max;
    end
    
    f = [0,0,0,0,1];

    Aeq = []; beq = [];
    lb = [0, 0, 0, 0, 0]; ub = [inf, inf, inf, inf, inf];
end

function [x,y] = getDataFromAandb(A,b)
    n = length(A) / 2;
    x = zeros(n,1);
    y = zeros(n,1);

    for i = 1:n
        x(i) = A(2*i-1,1);
        y(i) = b(2*i-1);
    end
end

% Metoda Symplex do minimalizacji funkcji f.
% Działa tak samo jak linprog, ale ma więcej założeń
% Zakłada bowiem postać zadania programowania liniowego taką jak w raporcie
% Czyli m.in.:
% 1. zadanie maksymalizacji piątej zmiennej
% 2. Pięć zmiennych nieujemnych
% 3. Bazą zą pozostałe zmienne od 4 do 5+2*n
% 4. Wektor b jest nieujemny
% 5. Macierz A jest szczególnej postaci i wymiarów 2*n na 5.
function [x_odp] = mySymplex(A,b)
    n = length(A) / 2;

    c = zeros(1, 5+2*n);
    c(5) = 1;

    M = zeros(2*n, 5+2*n); % główna macierz w tabelce
    % Pierwsze 3 kolumny to zmienne a, b, eps'. Kolejne 2n kolumn to
        % dodatkowe zmienne
    M(1:(2*n), 1:5) = A;
    M(1:(2*n), 6:(5+2*n)) = eye(2*n);

    baza = 6:(5+2*n); % początkowe zmienne bazowe.
    kosztyBazy = zeros(2*n, 1); % koszty zmiennych bazowych.

    z = kosztyBazy' * M;
    cMinusZ = c - z;

    % główna pental algorytmu Simplex
    while any(cMinusZ > 0)
        newBaseElement = find(cMinusZ == max(cMinusZ));
        newBaseElement = newBaseElement(1); % Na wypadek remisu, wybierz pierwszy
        
        j = wybierzElementDoWyrzuceniaZBazy(M(1:(2*n), newBaseElement), b);

        baza(j) = newBaseElement;
        kosztyBazy(j) = c(newBaseElement);

        b(j) = b(j) / M(j,newBaseElement); % to najpierw, bo za chwile M(j,newBaseElement) będzie równe 1.
        M(j,1:(5+2*n)) = M(j,1:(5+2*n)) / M(j,newBaseElement);
        for i = 1:(2*n)
            if (i == j) % ten wiersz był zrobiony przed petlą
                continue
            end
            b(i) = b(i) - b(j) * M(i,newBaseElement);
            M(i,1:(5+2*n)) = M(i,1:(5+2*n)) - M(j,1:(5+2*n)) * M(i,newBaseElement);
        end

        z = kosztyBazy' * M;
        cMinusZ = c - z;
    end
    % pentla się zakończyła, czyli mam punkt optymalny. Należy odczytać odpowiedź

    x_odp = zeros(5,1);
    for i = 1:5
        index = find(baza == i);
        if (~isempty(index))
            x_odp(i) = b(index);
        end 
    end
end

% p - pivot column
function [j] = wybierzElementDoWyrzuceniaZBazy(p, b)
    ratio = b./p;

    j = 1;
    currentMin = ratio(1);
    for i = 2:length(p)
        if ((p(i) > 0) && (ratio(i) < currentMin)) % Może być skrajny przypadek, że wszystkie są 0, ale to by oznaczało, że zadanie jest nieograniczone
            j = i;
            currentMin = ratio(i);
        end
    end
end

function [out] = singleTest(n)
    [x,y] = drawData(n);
    [A,b,f,Aeq,beq,lb,ub] = getAandb(x,y);

    options = optimoptions('linprog','Display','none');
    x_odp_MATLAB = linprog(-f,A,b,Aeq,beq,lb,ub,options);
    error_MATLAB = x_odp_MATLAB(5);
    
    x_odp_my = mySymplex(A,b);
    error_my = x_odp_my(5);
    
    out = (abs(error_my - error_MATLAB) < 0.00000001);
end

function [numOfSuccessfulTests] = multipleTests(n, numOfTests)
    numOfSuccessfulTests = 0;
    for i = 1:numOfTests
        disp(i)
        if (singleTest(n))
            numOfSuccessfulTests = numOfSuccessfulTests + 1;
        end
    end
end

function [] = plot_solution(x, y, x_odp)
    a = x_odp(1) - x_odp(2);
    b = x_odp(3) - x_odp(4);

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


