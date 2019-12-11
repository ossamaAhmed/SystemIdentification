k = 0:1:9;
u = cos(k);
Un1 = DFT(u);
Un2 = fft(u);
result = Un1 - Un2;

