function dangle = rocking_increment_APS(S)
fprintf('rocking_increment_APS');
%% Declaring variables from shape structure
rocking_angle = S.rocking_angle;
if isfield(S, 'Q_lab') 
    Q_lab = S.Q_lab;
else
    s_0 = 2*pi/S.lambda*[0; 0; 1];
    Q_lab = S.R_dqp_12*s_0 - s_0; % O.Q_lab is the q vector for centre of detector.
end

% instrument parameters
lambda = S.lambda;
d = S.d;
D = S.D;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recursion to find dtheta or dphi
if strcmp(rocking_angle, 'dtheta')
    dtheta = 0.0001; % starting guess for dtheta
    dtheta_increment = 0.000001; % dphi loop increment
    LHS = norm(roty(dtheta)*Q_lab-Q_lab); % magnitude of q'_3
    RHS = 2*pi()/lambda*d/D; % magnitude of q'_1 or q'_2
    
    % while loop to recursively solve for dtheta setting the magnitude of q'_3 to equal the magnitude of q'_1/2
    counter = 0;
    counter_limit = 1e5;
    fprintf('\n...calculating dtheta rocking angle...');
    while abs(LHS - RHS)>1e3 && counter < counter_limit        
        counter = counter + 1; % increasing counter for each loop iteration
        dtheta = dtheta + dtheta_increment; % adding an increment to dtheta
        LHS = norm(roty(dtheta)*Q_lab-Q_lab); % magnitude of q'_3
        RHS = 2*pi()/lambda*d/D; % magnitude of q'_1 or q'_2
    end
    if counter == counter_limit
        % if loop exceeds maximum number of iterations, dtheta is set to an arbitrary number to allow the rest of the code to run
        fprintf('\ndtheta calculation has exceeded iteration limit. dtheta set to 0.0025');
        dangle = 0.0025;
    else
        dangle = dtheta;
    end
elseif strcmp(rocking_angle, 'dphi')
    dphi = 0.0001; % starting guess for dphi
    dphi_increment = 0.000001; % dphi loop increment
    LHS = norm(rotx(dphi)*Q_lab-Q_lab); % magnitude of q'_3
    RHS = 2*pi()/lambda*d/D; % magnitude of q'_1 or q'_2
    
    % while loop to recursively solve for dphi setting the magnitude of q'_3 to equal the magnitude of q'_1/2
    counter = 0;
    counter_limit = 1e5;
    fprintf('\n...calculating dphi rocking angle...');
    while abs(LHS - RHS)>1e3 && counter < counter_limit
        counter = counter + 1; % increasing counter for each loop iteration
        dphi = dphi + dphi_increment; % adding an increment to dphi
        LHS = norm(rotx(dphi)*Q_lab-Q_lab); % magnitude of q'_3
        RHS = 2*pi()/lambda*d/D; % magnitude of q'_1 or q'_2
    end
    if counter == counter_limit
        % if loop exceeds maximum number of iterations, dphi is set to an arbitrary number to allow the rest of the code to run
        fprintf('\ndphi calculation has exceeded iteration limit. dphi set to 0.0025');
        dangle = 0.0025;
    else
        dangle = dphi;
    end
end
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end