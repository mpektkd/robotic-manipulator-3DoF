%Trajectory Planning for 3-DOF robotic arm
%First, we plan trajectory for each axis x,y, as z is constant

function coeff = parameters(xA, xB, D, half_tf)
%first phase-5th order polynomial
%%%% Firstly, we are going to design the first part of motion %%%%
x2dot = (xB - xA)/(half_tf - 2*D);
xAdot = x2dot*D + xA;
xBdot = - x2dot*D + xB ;

A1 = [0 0 0 0 0 1;
    (2*D)^5 (2*D)^4 (2*D)^3 (2*D)^2 (2*D) 1;
    0 0 0 0 1 0;
    5*(2*D)^4 4*(2*D)^3 3*(2*D)^2 2*(2*D) 1 0;
    0 0 0 2 0 0 ;
    20*(2*D)^3 12*(2*D)^2 6*(2*D) 2 0 0];
B1 = [xA ; xAdot ; 0 ; x2dot ; 0 ; 0];
x1 = A1\B1;

A2 = [2*D 1;
    1 0];
B2 = [xAdot ; x2dot];

x2 = A2\B2;

A3 = [(half_tf-2*D)^5 (half_tf-2*D)^4 (half_tf-2*D)^3 (half_tf-2*D)^2 (half_tf-2*D) 1;
    half_tf^5 half_tf^4 half_tf^3 half_tf^2 half_tf 1;
    5*(half_tf-2*D)^4 4*(half_tf-2*D)^3 3*(half_tf-2*D)^2 2*(half_tf-2*D) 1 0;
    5*half_tf^4 4*half_tf^3 3*half_tf^2 2*half_tf 1 0;
    20*(half_tf-2*D)^3 12*(half_tf-2*D)^2 6*(half_tf-2*D) 2 0 0;
    20*half_tf^3 12*half_tf^2 6*half_tf 2 0 0];

B3 = [xBdot ; xB ; x2dot ; 0 ; 0; 0];

x3 = A3\B3;
%%%% Secodnly, we are going to design the second part of motion (Palindrome)%%%%
A4 = [(half_tf+2*D)^5 (half_tf+2*D)^4 (half_tf+2*D)^3 (half_tf+2*D)^2 (half_tf+2*D) 1;
    half_tf^5 half_tf^4 half_tf^3 half_tf^2 half_tf 1;
    5*(half_tf+2*D)^4 4*(half_tf+2*D)^3 3*(half_tf+2*D)^2 2*(half_tf+2*D) 1 0;
    5*half_tf^4 4*half_tf^3 3*half_tf^2 2*half_tf 1 0;
    20*(half_tf+2*D)^3 12*(half_tf+2*D)^2 6*(half_tf+2*D) 2 0 0;
    20*half_tf^3 12*half_tf^2 6*half_tf 2 0 0];

B4 = [xBdot ; xB ; -x2dot ; 0 ; 0; 0];

x4 = A4\B4;

A5 = [half_tf+2*D 1;
    1 0];
B5 = [xBdot ; -x2dot];

x5 = A5\B5;

A6 = [(2*half_tf-2*D)^5 (2*half_tf-2*D)^4 (2*half_tf-2*D)^3 (2*half_tf-2*D)^2 (2*half_tf-2*D) 1;
    (2*half_tf)^5 (2*half_tf)^4 (2*half_tf)^3 (2*half_tf)^2 (2*half_tf) 1;
    5*(2*half_tf-2*D)^4 4*(2*half_tf-2*D)^3 3*(2*half_tf-2*D)^2 2*(2*half_tf-2*D) 1 0;
    5*(2*half_tf)^4 4*(2*half_tf)^3 3*(2*half_tf)^2 2*(2*half_tf) 1 0;
    20*(2*half_tf-2*D)^3 12*(2*half_tf-2*D)^2 6*(2*half_tf-2*D) 2 0 0;
    20*(2*half_tf)^3 12*(2*half_tf)^2 6*(2*half_tf) 2 0 0];

B6 = [xAdot ; xA ; -x2dot ; 0 ; 0; 0];

x6 = A6\B6;


coeff = {x1, x2, x3, x4, x5, x6};

end

