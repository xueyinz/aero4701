function mtr = dcmQuat(in)
% Direction cosine matrix from quaternions.
% Takes in a [4x1] quaternion vector and outputs [3x3] DCM matrix.
% out = DCM(in)

q0 = in(:, 1);
q1 = in(:, 2);
q2 = in(:, 3);
q3 = in(:, 4);

mtr(1, 1, :) = q0.^2 + q1.^2 - q2.^2 - q3.^2;
mtr(1, 2, :) = 2.*(q1.*q2 + q0.*q3);
mtr(1, 3, :) = 2.*(q1.*q3 - q0.*q2);

mtr(2, 1, :) = 2.*(q1.*q2 - q0.*q3);
mtr(2, 2, :) = q0.^2 - q1.^2 + q2.^2 - q3.^2;
mtr(2, 3, :) = 2.*(q2.*q3 + q0.*q1);

mtr(3, 1, :) = 2.*(q0.*q2 + q1.*q3);
mtr(3, 2, :) = 2.*(q2.*q3 - q0.*q1);
mtr(3, 3, :) = q0.^2 - q1.^2 - q2.^2 + q3.^2;

end