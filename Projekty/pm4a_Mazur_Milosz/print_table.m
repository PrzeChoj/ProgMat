function print_table(AA, B, fC, cB, xB, z, z_c, z_opt)


    table = [xB cB AA B];
    table = [
                NaN NaN fC NaN;
                table;
                NaN NaN z z_opt;
                NaN NaN z_c NaN;
            ];

    T = array2table(table, ...
                    'VariableNames',{'xB','cB','x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11','x12', 'x13','x14','x15','x16','x17','x18','x19','x20' 'B'}, ...
                    'RowName',{'fC','r1', 'r2', 'r3', 'r4', 'r5', 'z', 'z-c'} ...
                    );

    T

end