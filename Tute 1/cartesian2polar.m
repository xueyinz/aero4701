function POLAR = cartesian2polar(cartesian)

x = cartesian(1);
y = cartesian(2);
z = cartesian(3);
r = norm(cartesian);
phi = rad2deg(asin(-z/r));
theta = rad2deg(atan2(y, x));
POLAR = [r; phi; theta];