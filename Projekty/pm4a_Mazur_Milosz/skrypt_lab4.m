clear();

[A, b, c, d, g] = random_task(5, 10);

opcje=optimset(@linprog);
opcje=optimset(opcje, 'Display', 'iter', 'Algorithm', 'dual-simplex');

[x, fval, exitflag, output, lambda] = linprog(c, -A, -b, [], [], d, g, opcje);

x

[ZPx, ZDy, ZDexitflag] = sympleks(c, A, b, d, g, 1);


x
ZPx