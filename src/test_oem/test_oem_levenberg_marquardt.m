clear;

xa = textread('xa_t.txt');
y = textread('y_t.txt');
Sx = textread('Sx_t.txt');
Se = textread('Se_t.txt');

m = size(y,1);
n = size(xa,1);

H = zeros( n, n, m );
for i = 1:m
    H( :, :, i ) = textread( strcat( strcat( 'H_', int2str( i - 1 ) ), '_t.txt' ) );
end
J = textread('J_t.txt');

f = @( q, r, x, iter ) forward_model( J, H, q, r, x, iter );

O = oem;
O.linear = false;
O.G = true;
O.itermethod = 'ML';

t1 = cputime;
[X,r] = oem(O, struct, [], f, Sx, Se, [], [], xa, y);
t2 = cputime;

x = X.x;
G = X.G;

t = (t2 - t1) * 1000;

