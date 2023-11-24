% dualne: min (t(b)y_1 + t(g)y_2)
% t(A)y_1 + Iy_2 - Iy_3 = c
% y >= 0

% Bazę początkową robimy tak: Jak c jest dodatni, to weź y_2, a jak c jest
% ujemny, to weź y_3

% Uwaga: jak się zmienia znak w wierszu, to trzeba pamietać, żeby potem
% odwrócić z powrotem. Wspomnieć o tym w raporcie, żeby pokazać Ewci, że
% rozumiemy.

% NOTE: Do raportu można mniejsze n i m.

%% MyTest
[c, A, b, g] = drawData(n, m);

options = optimoptions('linprog','Display','none','Algorithm','dual-simplex');
[x,fval,exitflag_linprog,output,lambda] = linprog(c, A, b, [], [], zeros(1, n), g, options)

%% MyTest2
[ZPx, ZDy, exitflag_my] = dualSimplex(c, A, b, g, true)

%% MyTest3
[n_rows, n_cols] = size(A');
A_dual = [A', eye(n_rows), -eye(n_rows)];
c_dual = -[b', g', zeros(1, n_rows)];
b_dual = c';
baseIndexes = zeros(n_rows, 1);
signChanged = zeros(n_rows, 1);

% Zamiana wierszy na ujemny i dopisanie odpowiedniego do bazy:
for i = 1:n_rows
    if(b_dual(i) < 0)
        signChanged(i) = 1;
        baseIndexes(i) = n_cols + n_rows + i;
        A_dual(i, :) = -A_dual(i, :);
        b_dual(i, :) = -b_dual(i, :);
    else
        baseIndexes(i) = n_cols + i;
    end
end

disp("Moj:")
[ZDy, ~, exitflag, A_dual, ~, baseIndexes] = simplex(c_dual, A_dual, b_dual, baseIndexes, true)
    

%% Test
n = 5; % Długość wektora c i g
m = 10; % Liczba wierszy macierzy A i długość wektora b

exitflag_linprog = -2;

while (exitflag_linprog == -2)
    [c, A, b, g] = drawData(n, m);
    
    options = optimoptions('linprog','Display','none','Algorithm','dual-simplex');
    [x,fval,exitflag_linprog,output,lambda] = linprog(c, A, b, [], [], zeros(1, n), g, options)
    
    [ZPx, ZDy, exitflag_my] = dualSimplex(c, A, b, g, true)
end


%% Definicje użytych funkcji

% Metoda simplex do maksymalizacji t(c) * x
% Rozwiązuje metodą dualnego simplexu.
function [ZPx, ZDy, exitflag] = dualSimplex(c, A, b, g, verbose)
    if(verbose)
        disp("Dane wejsciowe: ");
        disp(c);
        disp(A);
        disp(b);
        disp(g);
    end
    exitflag = 1;
    [n_rows, n_cols] = size(A');
    A_dual = [A' eye(n_rows) -eye(n_rows)];
    c_dual = -[b', g', zeros(1, n_rows)];
    b_dual = c';
    baseIndexes = zeros(n_rows, 1);
    signChanged = zeros(n_rows, 1);
    
    % Zamiana wierszy na ujemny i dopisanie odpowiedniego do bazy:
    for i = 1:n_rows
        if(b_dual(i) < 0)
            signChanged(i) = 1;
            baseIndexes(i) = n_cols + n_rows + i;
            A_dual(i, :) = -A_dual(i, :);
            b_dual(i, :) = -b_dual(i, :);
        else
            baseIndexes(i) = n_cols + i;
        end
    end
    if(verbose)
        disp("Dane przekrztalcone do dualnego:")
        disp(signChanged);
        disp(A_dual);
        disp(b_dual);
        disp(c_dual);
    end

    % wywołanie zwykłego solvera simplex
    [ZDy, ~, exitflag, A_dual, ~, baseIndexes] = simplex(c_dual, A_dual, b_dual, baseIndexes, verbose);
    if(verbose)
        disp("Wyniki obliczen:")
        disp(linprog(-c_dual, [], [], A_dual, b_dual, zeros(size(c_dual))')')
        disp(ZDy)
    end
    
    % Dostosowanie wyniku dualnego do wyniku prymalnego:
    A_solution = A_dual(:, baseIndexes);
    for i = 1:n_rows
        if(signChanged(i) == 1)
            A_solution(:, i) = -A_solution(:, i);
        end
    end
    ZPx = c_dual(baseIndexes) * A_solution; % tu rozwiązanie prymalnego
end

% Metoda simplex do maksymalizacji t(c) * x
% Można jej podać z jakich indeksów ma zaczynać
function [x, fval, exitflag, A, bf, baseIndexes] = simplex(c, A, b, baseIndexes, verbose)
    exitflag = 1;
    maxIter = 100;
    bf = b;
    cb = c(baseIndexes); % Pobierz współczynniki z c odpowiadające bazie
    z = cb * A; % Wskaźnik funkcji celu
    zc = z - c; % Wskaźnik funkcji celu dla nowej bazy
    
    for iter = 1:maxIter
        if(verbose)
            disp('Obecna tabela:')
            disp([A bf; z NaN; zc NaN])
        
            disp('Obecne indeksy w bazie:')
            disp(baseIndexes)
        end
        if all(zc >= 0) % Sprawdzenie warunku zakończenia
            if(verbose)
                disp("Znaleziono rozwizanie")
            end
            break;
        end
    
        [~,col] = min(zc); % Znajdź kolumnę o najmniejszym współczynniku
        if not(any(A(:,col) > 0))
            disp("Nie ma rozwiazania problemu dualnego") % Ten simplex jest wywolywany dla problemu dualnego
            x = inf;
            fval = inf;
            exitflag = 0;
            return
        end
    
        mySet = (bf)./A(:,col); % Znajdź ilorazy b/wartość_kolumny dla elementów dodatnich w kolumnie
    
        mySet(A(:,col) < 0) = inf; % Ujemne mnie nie interesują
        [~,row] = min(mySet); % Znajdź indeks wiersza o najmniejszym ilorazie
    
        bf(row) = bf(row) / A(row, col); % Zaktualizuj wartość w bazie
        A(row, :) = A(row, :) / A(row, col); % Znormalizuj wiersz
      
        t = A(:, col) ./ A(row, :); % Oblicz współczynniki dla pozostałych wierszy
        t =  t(:, col);
        t(row) = 0;
    
        Dif = t * A(row, :);
        A = A - Dif;
        
        bf = bf - bf(row) * t;
    
        baseIndexes(row) = col; % Zaktualizuj indeks w bazie
        cb = c(baseIndexes)'; % Zaktualizuj współczynniki c odpowiadające bazie
    
        % Oblicz wskaźnik funkcji celu dla nowej bazy
        z = cb' * A;
        zc = z - c;
    end
    
    % Odczytajmy RO
    x = zeros(1, size(A, 2));
    x(baseIndexes) = bf';
    fval = sum(c.*x); % Oblicz wartość funkcji celu
end

% Funkcja drawData generuje dane do testów.
function [c, A, b, g] = drawData(n, m)
    % Generowanie wektora c
    c = randi([-5, 5], 1, n);
    
    % Generowanie macierzy A
    A = randi([-5, 5], m, n);
    
    % Generowanie wektora b
    b = randi([-5, 5], m, 1);
    
    % Generowanie wektora g
    g = randi([1, 30], n, 1);
end





















