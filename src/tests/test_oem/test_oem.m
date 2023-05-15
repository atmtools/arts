clear;

xa = textread('xa_t.txt');
y = textread('y_t.txt');
J = textread('J_t.txt');
Sa = textread('Sa_t.txt');
Se = textread('Se_t.txt');

f = @( q, r, x, iter ) linear_forward_model( J, q, r, x, iter);

O = oem
O.linear = true;
O.G = true;

% Run OEM
t1 = cputime;
[X,r] = oem(O, struct, [], f, Sa, Se, [], [], xa, y);
t2 = cputime;

x = X.x;
G = X.G;
t = (t2 - t1) * 1000;

