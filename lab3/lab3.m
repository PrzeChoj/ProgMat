%% Tests 2
j = wybierzElementDoWyrzuceniaZBazy(M(1:(2*n), newBaseElement), b)

%% Test 3
%x = [1;2;3]
%y = [10;7;5]
%[A_0,b_0] = getAandb(x,y)
n = 5;
[x,y] = drawData(n);
[A,b,f,Aeq,beq,lb,ub] = getAandb(x,y);
x_odp = mySymplex(A,b);
plot_solution(x, y, x_odp(1), x_odp(2))

%% Tests
%n = 10;
[x,y] = drawData(n)
x = [1;2;3]
y = [10;7;5]
[A,b,f,Aeq,beq,lb,ub] = getAandb(x,y)

[x_odp,fval,exitflag,output,lambda] = linprog(-f,A,b,Aeq,beq,lb,ub)

e = max(y) + fval;
%abs(max(abs(x_odp(1)*x + x_odp(2) - y)) - e) < 0.00000001 % powinna byc rownosc
plot_solution(x, y, x_odp(1), x_odp(2))

function [x, y] = drawData(n)
    x=(1:1:n)';
    y=randi([0 1500],n,1);
end

function [A,b,f,Aeq,beq,lb,ub] = getAandb(x,y)
    n = length(x);
    A = zeros(2*n, 3);
    b = zeros(2*n, 1);

    y_max = max(y);

    for i = 1:n
        A(2*(i-1)+1,:) = [x(i), 1, 1];
        b(2*(i-1)+1,:) = y(i) + y_max;
    
        A(2*(i-1)+2,:) = [-x(i), -1, 1];
        b(2*(i-1)+2,:) = -y(i) + y_max;
    end
    
    f = [0,0,1];

    Aeq = []; beq = [];
    lb = [-inf -inf 0]; ub = [inf inf inf];
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
% 1. zadanie maksymalizacji trzeciej zmiennej
% 2. Trzecia zmienna nieujemna. Pierwsze dwie dowolne rzeczywiste
% 3. Bazą zą pozostałe zmienne od 4 do 3+2*n
% 4. Wektor b jest nieujemny
% 5. Macierz A jest szczególnej postaci i wymiarów 2*n na 3.
function [x_odp] = mySymplex(A,b)
    n = length(A) / 2;

    c = zeros(1, 3+2*n);
    c(3) = 1;

    M = zeros(2*n, 3+2*n); % główna macierz w tabelce
    % Pierwsze 3 kolumny to zmienne a, b, eps'. Kolejne 2n kolumn to
        % dodatkowe zmienne
    M(1:(2*n), 1:3) = A;
    M(1:(2*n), 4:(3+2*n)) = eye(2*n);

    baza = 4:(3+2*n); % początkowe zmienne bazowe.
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
        M(j,1:(3+2*n)) = M(j,1:(3+2*n)) / M(j,newBaseElement);
        for i = 1:(2*n)
            if (i == j) % ten wiersz był zrobiony przed petlą
                continue
            end
            b(i) = b(i) - b(j) * M(i,newBaseElement);
            M(i,1:(3+2*n)) = M(i,1:(3+2*n)) - M(j,1:(3+2*n)) * M(i,newBaseElement);
        end

        z = kosztyBazy' * M;
        cMinusZ = c - z;
    end
    % pentla się zakończyła, czyli mam punkt optymalny. Należy odczytać odpowiedź

    x_odp = zeros(3,1);
    indexA = find(baza == 1);
    if (~isempty(indexA))
        x_odp(1) = b(indexA);
    end
    indexB = find(baza == 2);
    if (~isempty(indexB))
        x_odp(2) = b(indexB);
    end
    indexEps = find(baza == 3);
    if (~isempty(indexEps))
        x_odp(3) = b(indexEps);
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

    [x_odp,fval] = linprog(f,A,b,Aeq,beq,lb,ub);
    
    e = max(y) + fval;
    out = (abs(max(abs(x_odp(1)*x + x_odp(2) - y)) - e) < 0.00000001); % TODO: Porownanie z moim algorytmem
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


