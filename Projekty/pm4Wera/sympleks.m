function [x, fval, exitflag, A, bf, indices] = sympleks(f,A,b,indices)
exitflag = 1;
MAX_ITER = 10000;
c=f;
bf=b;
cb = c(indices);
z = cb * A;
zc = z - c;

% aby uniknąć degeneracji, dodajemy niewielki epsilon przy wyborze 
eps = (0.00001).*(ones(size(b, 1), 1));

for iter = 1:MAX_ITER

    disp('Obecna tabela:')
    disp([A bf; z inf;zc inf])

    disp('Obecne indeksy w bazie:')
    disp(indices)
    if all(zc>=0)
        disp("Znaleziono rozwizanie")
        break;
    end

    % znajdz minimalną wartość w wierszu zc
    [~,col]=min(zc);
    if not(any(A(:,col)>0))
        disp("Nie można znaleźć rozwiązania")
        x = inf;
        fval = inf;
        exitflag = 0;
        return
    end

    %oblicz wartości b/A dla kolumny z wybraną wartością zc
    set=(bf + eps)./A(:,col);

    %znajdz najmniejszą nieujemna wartość w kolumnie
    set(A(:,col)<0)=inf;
    [~,row]=min(set);

    % sprowadzamy A(row, col) do 1
    bf(row)=bf(row)/A(row,col);
    A(row,:)=A(row,:)/A(row,col);
  
    t = A(:,col)./A(row,:);
    t =  t(:, col);
    t(row) = 0;

    % wyzeruj pozostale wspólczynniki w kolumnie
    Dif=t*A(row,:);
    A=A-Dif;
    
    % przemnóż przez to samo dla prawej kolumny
    bf=bf-bf(row)*t;


    % zaaktualizuj lewą kolumnę (przy zmiennych bazowych)
    indices(row)=col;
    cb= c(indices)';
    

    % oblicz z i z-c
    z=cb' * A;
    zc=z-c;
end
% wskaż RO i wartość funkcji
x = zeros(1, size(A, 2));
x(indices) = bf';
fval = sum(c.*x);
end