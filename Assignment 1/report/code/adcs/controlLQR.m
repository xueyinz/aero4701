function [m,bdot] = controlLQR(xhat,qint,mode,cont,b,bm1,bdotm1)

% Initialise the bdot vector in case it isn't calculated
bdot = zeros(3,1);

switch mode
    case 'point'

        % Find the control input
        mtot = -cont.K*xhat(1:6);
        mint = cont.Ki*qint;
        
        % Find the total magnetic moment
        m = cross(mtot+mint,b);
        
    case 'bdot'
        
        % Calcualte bdot via filtering
        bdot = (1/cont.tf.den{1}(1))*(cont.tf.num{1}(1)*b + ...
            cont.tf.num{1}(2)*bm1 - cont.tf.den{1}(2)*bdotm1);
        
        % Calculate the magnetic moment
        m = -cont.bdotK*bdot./norm(b)^2;

    case 'none'
        
        % Output an array of zeros
        m = zeros(3,1);
end


