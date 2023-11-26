function [A, b, c, d, g] = random_task(n, m)

    A = (rand(m,n) - 0.5) * 10;
    b = rand(1, m) * 4 + 1;
    c = (rand(1, n) - 0.5 ) * 10;

    d = rand(1, n) * -29 - 1;
    g = rand(1, n) * 29 + 1;

end