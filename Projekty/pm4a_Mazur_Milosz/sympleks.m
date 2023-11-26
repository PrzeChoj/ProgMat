function [ZPx, ZDy, ZDexitflag] = sympleks(c, A, b, d, g, draw)
    
    % zadanie dualne z zadania -> min c, A >= b, d <= x <= g
    % postać zadania dualnego: max [b, -g, d] , [A', -I, I] == c, y >= 0
    
    s = size(A);
    AD = [A', -eye(s(2)), eye(s(2))];
    bD = c;
    cD = [b, -g, d];
    
    s = size(A');

    % max cD
    % y >= 0
    % baza dodatnia więc gdzie dodatnie c to z d a gdzie ujemne to z g

    %% Wybór zmiennych bazowych
    xBD = [];
    cBD = [];
    minus = [];

    for i = 1:length(bD)
        if bD(i) >= 0
            xBD = [xBD, s(2) + s(1) + i];
            cBD = [cBD, cD(s(2) + s(1) + i)];
        else
            xBD = [xBD, s(2) + i];
            cBD = [cBD, cD(s(2) + i)];
            minus = [minus, i];
            AD(i,:) = -AD(i,:);
            bD(i) = -bD(i);
        end
    end
    
    %% Początek zadania

    % Zmienna exitflag:
    %  1 - optymalne
    % -3 - nieograniczone
    % -2 - sprzeczne
    % -1 - zbyt duzo iteracji

    % Macierz
    B = bD';
    AA = AD;

    % Funkcja celu
    fC = cD;
    % czeęść bazowa funkcji celu
    cB = cBD';
    % zmienne bazowe
    xB = xBD';
    % optymalne z
    z_opt = cB' * B;
    % wektor z
    z = cB' * AA;
    % wskaznik optymalnosci 
    z_c = z - fC;

    exitflag = -1;
    
    if draw == 1
        fprintf("Startowa tabela zadania pomocniczego:\n");
        print_table(AA, B, fC, cB, xB, z, z_c, z_opt);
    end

    for iter = 1:50
        % z-c
        neg_z_c = find(z_c < 0);

        if isempty(neg_z_c)
            exitflag = 1;
            iterations = iter - 1;
            break; % z-c >=0, koniec
        end

        % Czy AA jest dodatnie dla z_c < 0
        for i=1:length(neg_z_c)
            pos_a_k = find(AA(:,neg_z_c(i)) > 0);

            if isempty(pos_a_k) 
               if draw == 1
                    fprintf("Kolumna %g spelnia warunki nieograniczonosci. \n", neg_z_c(i));
               end
               exitflag = -3; % zadanie nieograniczone, koniec
               iterations = iter - 1;
               break;
            end
        end
        if exitflag == -3 
            if draw == 1
                fprintf("Zadanie nieograniczone. \n");
            end
            ZPx = [];
            ZDy = [];
            ZDexitflag = exitflag;
            return;
        end
        
        %% Algorytm Simplex

        % kolumna do wymiany
        k = find(z_c==min(z_c),1);

        % zmiana bazy
        A_k = AA(:,k);
        A_k(A_k<0) = 0; % bierzemy tylko dodatnie

        bA = B ./ A_k;
        r = find(bA == min(bA)); % minimalna wartosc, wiersz do wymiany
        r = r(1);

        b_r = B(r);
        a_r_k = AA(r,k);

        if draw == 1
            fprintf("Wybieramy wiersz %g i kolumne %g\n",r, k);
            fprintf("Wchodzi x%g za zmienna x%g\n", k, xB(r));
            fprintf("Wartosci a_r_k: %g  b_r: %g\n",a_r_k, b_r);
        end

        % wymiana zmiennych bazy 
        xB(r) = k;
        cB(r) = fC(k);

        % aktualizacja A
        A_old = AA;
        s_AA = size(AA);
        for i=1:s_AA(1)
            for j=1:length(fC)
                if i == r
                    AA(i,j) = A_old(i,j)/a_r_k;
                else
                    AA(i,j) = A_old(i,j) - (A_old(r,j)/a_r_k)*A_old(i,k);
                end
            end
        end
        % aktualizacja B
        for i=1:s_AA(1)
            if i == r
                B(i) = b_r/a_r_k;
            else 
                B(i) = B(i) - ((A_old(i,k)/a_r_k)*b_r);
            end
        end
        
        % aktualizacja
        z_opt = cB' * B;
        z = cB' * AA;
        z_c = z - fC;

        if draw == 1
            print_table(AA, B, fC, cB, xB, z, z_c, z_opt);
        end
    end


    % jesli nie rozwiazano zadania 
    if exitflag ~= 1
        if draw == 1
            fprintf("Zbyt duzo iteracji");
        end
        ZPx = [];
        ZDy = [];
        ZDexitflag = -1;
        return
    end

    ZDexitflag = 1;

    %% Rozwiazanie zadania dualnego

    ZDy = zeros(1, length(fC));
    
    for i=1:length(xB)
        ZDy(xB(i)) = B(i);
    end

    if draw == 1
        fprintf("BRD zadanie dualne\n");
        ZDy
    end

    %% Rozwiazanie zadania pierwotnego
    
    A_b_1 = [];
    for i = 1:length(xBD)
        multi = 1;
        if ismember(i, minus)
            multi = -1;
        end

        A_b_1 = [A_b_1, multi * AA(:,xBD(i))];
    end

    if draw == 1
        fprintf("Macierz A_b^-1\n");
        A_b_1
        fprintf("Wektor c_B\n");
        cB'
    end

    ZPx = cB' * A_b_1;

    if draw == 1
        fprintf("BRD zadanie pierwotne\n");
        ZPx
    end

end