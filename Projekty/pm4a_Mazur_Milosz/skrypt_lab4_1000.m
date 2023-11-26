clear();

eps = 1e-5;

brak_roz = 0;
brak_roz_i_roz = 0;
roz_i_brak_roz = 0;
roz = 0;
roz_rozne = 0;

for i = 1:1000

    [A, b, c, d, g] = random_task(5, 10);

    opcje=optimset(@linprog);
    opcje=optimset(opcje, 'Display', 'off', 'Algorithm', 'dual-simplex');

    [x, fval, exitflag, output, lambda] = linprog(c, -A, -b, [], [], d, g, opcje);

    [ZPx, ZDy, ZDexitflag] = sympleks(c, A, b, d, g, 0);

    
    if isempty(x)
        if isempty(ZPx)
            brak_roz = brak_roz + 1;
        else
            brak_roz_i_roz = brak_roz_i_roz + 1;
        end
        continue;
    else
        if isempty(ZPx)
            roz_i_brak_roz = roz_i_brak_roz + 1;
            continue;
        end
    end

    if sum(abs(x' - ZPx) < eps) == 5
        roz = roz + 1; 
    else
        roz_rozne = roz_rozne + 1;
    end
end

brak_roz
brak_roz_i_roz
roz_i_brak_roz
roz
roz_rozne