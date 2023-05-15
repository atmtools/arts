clear;
A = textread('A_t.txt');
t1 = cputime;
B = A*A;
t2 = cputime;
t = (t2 - t1) * 1000;

