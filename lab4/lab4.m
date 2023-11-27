% dualne: min (t(b)y_1 + t(g)y_2)
% t(A)y_1 + Iy_2 - Iy_3 = c
% y >= 0

% Bazę początkową robimy tak: Jak c jest dodatni, to weź y_2, a jak c jest
% ujemny, to weź y_3

% Uwaga: jak się zmienia znak w wierszu, to trzeba pamietać, żeby potem
% odwrócić z powrotem. Wspomnieć o tym w raporcie, żeby pokazać Ewci, że
% rozumiemy.

% NOTE: Do raportu można mniejsze n i m.

MY_INF = 10^20;

%% MyTest
n = 5;
m = 10;
[c, A, b, g] = drawData(n, m);

options = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');
[x, fval, exitflag_linprog, output, lambda] = linprog(c, A, b, [], [], zeros(1, n), g, options)

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
        A_dual(i, :) = -A_dual(i, :); % zamiana wiersza na ujemny
        b_dual(i) = -b_dual(i);
    else
        baseIndexes(i) = n_cols + i;
    end
end

disp("Moj:")
[x, fval, exitflag, A_after, indices, zBaza] = simplex(c_dual, A_dual, b_dual, baseIndexes, false)
disp("linprog:")
[x, fval, exitflag_linprog, output, lambda] = linprog(-c_dual, [], [], A_dual, b_dual, zeros(1, size(c_dual, 2)))

% odczytanie wyniku
% Chcemy rozwiazać ukałd równań x * B = d,
% gdzie x jest niewiadomą
% B jest macierzą oryginalną A' na elementach starej bazy baseIndexes
A_solution = A_after(:, baseIndexes);
for i = 1:n_rows
    if(signChanged(i) == 1)
        A_solution(:, i) = -A_solution(:, i);
    end
end
my_x = c_dual(indices) * A_solution
% niestety sa tu ujemne :<

% Podejście z odwracaniem
naw_A = [A', eye(n_rows), -eye(n_rows)];
A_B = naw_A(:, indices);
inv(A_B) * c'



%% useful
j = wybierzElementDoWyrzuceniaZBazy(M(:, newBaseElement), b);


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

% Metoda simplex do zadania:
% max (t(c) * y)
% A * y == b
% y >= 0
% 
% Zakładamy, że baseIndexes ma w sobie indeksy nadające się do bazy
function [x, fval, exitflag, A, baza, zBaza] = simplex(c, A, b, baseIndexes, verbose)
    exitflag = 1;
    baza = baseIndexes;
    kosztyBazy = c(baza);
    M = A;

    z = kosztyBazy * M;
    cMinusZ = c - z;
    licznikIteracji = 0;

    while (any(cMinusZ > 0) && licznikIteracji < 100)
        licznikIteracji = licznikIteracji + 1;

        % Krok 1: Wybór nowej zmiennej do wejścia do bazy
        newBaseElement = find(cMinusZ == max(cMinusZ));
        newBaseElement = newBaseElement(1); % Na wypadek remisu, wybierz pierwszy
        
        % Krok 2: Wybór zmiennej do wyjścia z bazy
        j = wybierzElementDoWyrzuceniaZBazy(M(:, newBaseElement), b);

        % Krok 3: Aktualizacja bazy
        baza(j) = newBaseElement;
        kosztyBazy(j) = c(newBaseElement);

        % Krok 4: Przekształcenie wiersza zmiennych bazowych
        b(j) = b(j) / M(j,newBaseElement); % to najpierw, bo za chwile M(j,newBaseElement) będzie równe 1.
        M(j,:) = M(j,:) / M(j,newBaseElement);
        for i = 1:(size(M, 1))
            if (i == j) % ten wiersz był zrobiony przed petlą for
                continue
            end
            b(i) = b(i) - b(j) * M(i,newBaseElement);
            M(i,:) = M(i,:) - M(j,:) * M(i,newBaseElement);
        end

        % Krok 5: Obliczenie nowych kosztów
        z = kosztyBazy * M;
        cMinusZ = c - z;

        if (verbose)
            disp("Kolejna tabelka sympleksowa:")
            disp(M)
            disp("Aktualne zmienne bazowe:")
            disp(baza)
        end
    end
    % pentla się zakończyła, czyli mam punkt optymalny. Należy odczytać odpowiedź

    % Odczytanie rozwiązania
    x_odp = zeros(size(c, 2), 1);
    for i = 1:size(baza, 1)
        x_odp(baza(i)) = b(i);
    end

    fval = c * x_odp;

    if (verbose)
        disp("Rozwiązanie:")
        disp(x_odp)
    end

    x = x_odp;

    zBaza = z(baza);
    A = M;
end

% Funkcja wybierzElementDoWyrzuceniaZBazy wybiera element do usunięcia z
% bazy. Czyli to na którym jako pierwszym dochodzimy do granicy ograniczeń
% p - pivot column
function [j] = wybierzElementDoWyrzuceniaZBazy(p, b)
    ratio = b./p;

    ograniczajace = find(ratio >= 0);

    j = ograniczajace(1);
    currentMin = ratio(j);
    for i = 2:length(ograniczajace)
        o_i = ograniczajace(i);
        if ((p(o_i) > 0) && (ratio(o_i) < currentMin)) % Może być skrajny przypadek, że wszystkie są 0, ale to by oznaczało, że zadanie jest nieograniczone
            j = o_i;
            currentMin = ratio(o_i);
        end
    end
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





















