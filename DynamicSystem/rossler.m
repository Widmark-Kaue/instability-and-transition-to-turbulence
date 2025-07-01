function ros = rossler(t,X,a,b,c)

x = X(1); y = X(2); z = X(3);
dxdt = -y - z;
dydt = x+a*y;
dzdt = b + z*(x-c);
ros = [dxdt; dydt; dzdt];

