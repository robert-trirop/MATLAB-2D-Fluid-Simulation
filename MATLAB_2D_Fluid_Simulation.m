% MIT License
% 
% Copyright (c) 2024 Robert Schütze
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
% If you use this code for your own animations, please give credit to me.

% 2D Fluid Simulation in MATLAB
% Cleaned and commented version
% To increase performance decrease the number of Jacobi iterations, add less inflow
% particles or add them less often (e.g. only every 2nd simulation iteration)

% Simulation parameters
s = 200;                            % Grid size
ar = 2;                             % Aspect ratio
J = [0 1 0; 1 0 1; 0 1 0]/4;        % Stencil for Jacobi method

% Create a grid
[X, Y] = meshgrid(1:s*ar, 1:s);

% Initialize pressure and velocity fields
[p, vx, vy] = deal(zeros(s, s*ar));

% Initial positions of particles
[px, py] = meshgrid(10:15, 95:105);
px = reshape(px, numel(px), 1);
py = reshape(py, numel(py), 1);

% Save these initial positions for the inflow
pxo = px;
pyo = py;

f = figure(1); % Create figure to check if it's closed

% Main simulation loop (stops when closing figure window)
while ishandle(f)
    % Set initial velocity in a specific region
    vx(95:105, 10:15) = 0.1;
    
    % Compute right-hand side for pressure equation
    rhs = -divergence(vx, vy);
    
    % Jacobi iteration to solve for pressure
    % Higher number of iterations yields better solution
    for i = 0:100
        p = conv2(p, J, 'same') + rhs/2;
    end
    
    % Compute velocity gradient and update velocities for non-boundary pixels
    [dx, dy] = gradient(p);
    vx(2:end-1, 2:end-1) = vx(2:end-1, 2:end-1) - dx(2:end-1, 2:end-1);
    vy(2:end-1, 2:end-1) = vy(2:end-1, 2:end-1) - dy(2:end-1, 2:end-1);   
    
    % Advect velocity field using Runge-Kutta 4th order method (-1 = backward)
    [pvx, pvy] = RK4(X, Y, vx, vy, -1);
    vx = interp2(vx, pvx, pvy);
    vy = interp2(vy, pvx, pvy);  
    
    % Advect particles using Runge-Kutta 4th order method (1 = forward)
    [px, py] = RK4(px, py, vx, vy, 1);
    
    % Add the inflow particles
    px = [px; pxo];
    py = [py; pyo];
    
    % Visualization of particle positions
    scatter(px, py, 1, 'filled');
    axis equal; 
    axis([0 s*ar 0 s]);

    drawnow;
end

% Function for Runge-Kutta 4th order method for advection
function [x_new, y_new] = RK4(px, py, vx, vy, h)
   k1x = interp2(vx, px, py);
   k1y = interp2(vy, px, py);
   k2x = interp2(vx, px + h/2 * k1x, py + h/2 * k1y);
   k2y = interp2(vy, px + h/2 * k1x, py + h/2 * k1y);
   k3x = interp2(vx, px + h/2 * k2x, py + h/2 * k2y);
   k3y = interp2(vy, px + h/2 * k2x, py + h/2 * k2y);
   k4x = interp2(vx, px + h * k3x, py + h * k3y);
   k4y = interp2(vy, px + h * k3x, py + h * k3y);
   x_new = px + h/6 * (k1x + 2*k2x + 2*k3x + k4x);
   y_new = py + h/6 * (k1y + 2*k2y + 2*k3y + k4y);
end
