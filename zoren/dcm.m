function mtr = dcm(euler)

psi = euler(:, 1); theta = euler(:, 2); phi = euler(:, 3); 

mtr(1, 1, :) = cos(psi).*cos(theta);
mtr(1, 2, :) = sin(psi).*cos(theta);
mtr(1, 3, :) = -sin(theta);

mtr(2, 1, :) = cos(psi).*sin(theta).*sin(phi)-sin(psi).*cos(phi);
mtr(2, 2, :) = sin(psi).*sin(theta).*sin(phi)+cos(psi).*cos(phi);
mtr(2, 3, :) = cos(theta).*sin(phi);

mtr(3, 1, :) = cos(psi).*sin(theta).*cos(phi)+sin(psi).*sin(phi);
mtr(3, 2, :) = sin(psi).*sin(theta).*cos(phi)- cos(psi).*sin(phi);
mtr(3, 3, :) = cos(theta).*cos(phi);

end