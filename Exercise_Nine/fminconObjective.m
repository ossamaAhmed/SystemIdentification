function f = fminconObjective(x, nPar)
    f = sqrt(x(nPar+1:end)' * x(nPar+1:end));
end
