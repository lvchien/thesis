%% Moore-Penrose pseudoinverse of a number n

function inv = MPInverse(n)
    if n == 0
        inv = 0;
    else
        inv = 1.0/n;
    end;
end
