function result = periodogram(e)
    N = size(e, 2);
    En = DFT(e);
    result = 1/N * (abs(En) .^ 2);
end