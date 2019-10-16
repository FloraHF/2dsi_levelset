function [ data, g, data0 ] = sdsi3D(accuracy)
% sdsi3D: demonstrate the 3D defence-intrusion game
%
%   [ data, g, data0 ] = sdsi3D(accuracy)
%  
% In this example, the set for the intruder (player b) to reach is circle at the
% origin, the set to avoid is the defender's (player a) capture range.
%
% The relative coordinate dynamics are
%
%     \dot r1     = v_d \cos \phi
%	  \dot r2     = v_i \cos \psi
%	  \dot \theta = v_d / r1 \sin \phi - v_i / r2 \sin \psi
%
% where v_a and v_b are constants, input a is trying to reach the target
% and avoid the avoid set input a is trying to avoid the target and hit the
% avoid set
%
% This file is meant to be compared with my python code of geometric
% solution
%
% Input Parameters:
%
%   accuracy: Controls the order of approximations.
%     'low': Use odeCFL1 and upwindFirstFirst.
%     'medium': Use odeCFL2 and upwindFirstENO2 (default).
%     'high': Use odeCFL3 and upwindFirstENO3.
%     'veryHigh': Use odeCFL3 and upwindFirstWENO5.
%
% Output Parameters:
%
%   data: Implicit surface function at t_max.
%
%   g: Grid structure on which data was computed.
%
%   data0: Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 3/26/04
% Subversion tags for version control purposes.
% $Date: 2012-07-04 14:27:00 -0700 (Wed, 04 Jul 2012) $
% $Id: air3D.m 74 2012-07-04 21:27:00Z mitchell $

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 25;                  % End time.
plotSteps = 10;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 1;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

% What kind of dissipation?
dissType = 'global';


%---------------------------------------------------------------------------
% Problem Parameters, imported from python code
plotter = py.importlib.import_module('plotter');
Config = py.importlib.import_module('Config').Config;
%   targetRadius  Radius of target circle (positive).
%   velocityA	  Speed of the evader (positive constant).
%   velocityB	  Speed of the pursuer (positive constant).
%   inputA	  Maximum turn rate of the evader (positive).
%   inputB	  Maximum turn rate of the pursuer (positive).
level = 2;
filename = ['level', num2str(level), 'veryHigh', '.gif'];
targetRadius = double(Config.TAG_RANGE)+level;
captureRadius = double(Config.CAP_RANGE);
velocityA = double(Config.VD);
velocityB = double(Config.VI);
% inputA = 2*pi;
% inputB = 2*pi;

%---------------------------------------------------------------------------
% What level set should we view?
level = [0, level];

% Visualize the 3D reachable set.
displayType = 'contour';

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 1;

% Visualize the angular dimension a little bigger.
aspectRatio = [ 1 1 0.4 ];

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 0;

%---------------------------------------------------------------------------
% Approximately how many grid cells?
%   (Slightly different grid cell counts will be chosen for each dimension.)
Nx = 29;

% Create the grid.
g.dim = 4;
g.min = [ -1; -1; -1; -1 ]*targetRadius*3;
g.max = [  1; +1; +1; +1 ]*targetRadius*3;
g.bdry = { @addGhostExtrapolate; @addGhostExtrapolate; @addGhostExtrapolate; @addGhostExtrapolate };
% Roughly equal dx in x and y (so different N).
g.N = [ Nx; Nx; Nx; Nx ];
% Need to trim max bound in \psi (since the BC are periodic in this dimension).
% g.max(3) = g.max(3) * (1 - 1 / g.N(3));
g = processGrid(g);

ixd = (g.N(1)+1)/2;
iyd = ceil(0.8*g.N(2));

% For visualization, fix the location of the defender
dLoc = [g.xs{1}(ixd, 1, 1, 1), g.xs{2}(1, iyd, 1, 1)];
tgt = struct(plotter.get_target());
D0 = cell(plotter.get_defender(dLoc(2)));
dctr = struct(plotter.get_constd(dLoc(2)));

g_slice.dim = 2;
g_slice.min = [ g.min(3); g.min(4)];
g_slice.max = [ g.max(3); g.max(4)];
g_slice.bdry = { g.bdry{3}; g.bdry{4} };
g_slice.N = [ g.N(3); g.N(4) ];
g_slice = processGrid(g_slice);

%---------------------------------------------------------------------------
% Create initial conditions (cylinder centered on origin).
% data = reachAvoid(g, targetRadius, captureRadius);
data = max(avoid(g, captureRadius), reach(g, targetRadius));
data0 = data;
data0_slice = squeeze(data(ixd, iyd, :, :));

%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termLaxFriedrichs;
schemeData.hamFunc = @sdsi4DHamFunc;
schemeData.partialFunc = @sdsi4DPartialFunc;
schemeData.grid = g;

% The Hamiltonian and partial functions need problem parameters.
schemeData.velocityA = velocityA;
schemeData.velocityB = velocityB;
% schemeData.inputA = inputA;
% schemeData.inputB = inputB;

%---------------------------------------------------------------------------
% Choose degree of dissipation.

switch(dissType)
 case 'global'
  schemeData.dissFunc = @artificialDissipationGLF;
 case 'local'
  schemeData.dissFunc = @artificialDissipationLLF;
 case 'locallocal'
  schemeData.dissFunc = @artificialDissipationLLLF;
 otherwise
  error('Unknown dissipation function %s', dissFunc);
end

%---------------------------------------------------------------------------
if(nargin < 1)
  accuracy = 'veryHigh';
end

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.75, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  schemeData.derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  schemeData.derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  schemeData.derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  schemeData.derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Restrict the Hamiltonian so that reachable set only grows.
%   The Lax-Friedrichs approximation scheme MUST already be completely set up.
innerFunc = schemeFunc;
innerData = schemeData;
clear schemeFunc schemeData;

% Wrap the true Hamiltonian inside the term approximation restriction routine.
schemeFunc = @termRestrictUpdate;
schemeData.innerFunc = innerFunc;
schemeData.innerData = innerData;
schemeData.positive = 0;

%---------------------------------------------------------------------------
% Initialize Display
f = figure;
axis tight manual;

plot(double(D0{1}), double(D0{2}), 'color','r', 'linewidth', 2.3); hold on;
contour(double(tgt.X),double(tgt.Y),double(tgt.T), [0, 0], 'color','g', 'linewidth', 2.3);
if length(level) == 1
    level = [level(1), level(1)];
end
contour(double(dctr.X),double(dctr.Y),double(dctr.C), level, 'linewidth', 2, 'color', 'b');

% Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

h = visualizeLevelSet(g_slice, data0_slice, displayType, [0], [ 't = ' num2str(t0) ]);

camlight right;  camlight left;
hold on;
axis(g.axis);
daspect(aspectRatio);
drawnow;

frame = getframe(f); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
while(tMax - tNow > small * tMax)

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(tMax, tNow + tPlot) ];
  
  % Take a timestep.
  [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);
  data = max(avoid(g, captureRadius), data);

  if(pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end

  % Get correct figure, and remember its current view.
  figure(f);
  [ view_az, view_el ] = view;

  % Delete last visualization if necessary.
  if(deleteLastPlot)
    delete(h);
  end

  % Move to next subplot if necessary.
  if(useSubplots)
    plotNum = plotNum + 1;
    subplot(rows, cols, plotNum);
  end

  % Create new visualization.
  h = visualizeLevelSet(g_slice, squeeze(data(ixd, iyd, :, :)), displayType, [0], [ 't = ' num2str(tNow) ]);

  % Restore view.
  view(view_az, view_el);
  
  frame = getframe(f); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 
  imwrite(imind,cm,filename, 'DelayTime', 0.1, 'WriteMode','append'); 
  
end

endTime = cputime;
fprintf('Total execution time %g seconds\n', endTime - startTime);


%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function hamValue = sdsi4DHamFunc(t, data, deriv, schemeData)
% sdsi3DHamFunc: analytic Hamiltonian for sdsi reach-avoid game.
%
% hamValue = sdsi3DHamFunc(t, data, deriv, schemeData)
%
% This function implements the hamFunc prototype for the three dimensional
%   sdsi reach-avoid game.
%
% It calculates the analytic Hamiltonian for such a flow field.
%
% Parameters:
%   t            Time at beginning of timestep (ignored).
%   data         Data array.
%   deriv	 Cell vector of the costate (\grad \phi).
%   schemeData	 A structure (see below).
%
%   hamValue	 The analytic hamiltonian.
%
% schemeData is a structure containing data specific to this Hamiltonian
%   For this function it contains the field(s):
%
%   .grid	 Grid structure.
%   .velocityA	 Speed of the evader (positive constant).
%   .velocityB	 Speed of the pursuer (positive constant).
%
% Ian Mitchell 3/26/04

checkStructureFields(schemeData, 'grid', 'velocityA', 'velocityB');

grid = schemeData.grid;

hamValue = -(+ schemeData.velocityA * sqrt(deriv{1}.^2 + deriv{2}.^2) ...
             - schemeData.velocityB * sqrt(deriv{3}.^2 + deriv{4}.^2));



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function alpha = sdsi4DPartialFunc(t, data, derivMin, derivMax, schemeData, dim)
% sdsi3DPartialFunc: Hamiltonian partial fcn for 3D collision avoidance example.
%
% alpha = air3DPartialFunc(t, data, derivMin, derivMax, schemeData, dim)
%
% This function implements the partialFunc prototype for the three dimensional
%   aircraft collision avoidance example (also called the game of
%   two identical vehicles).
%
% It calculates the extrema of the absolute value of the partials of the 
%   analytic Hamiltonian with respect to the costate (gradient).
%
% Parameters:
%   t            Time at beginning of timestep (ignored).
%   data         Data array.
%   derivMin	 Cell vector of minimum values of the costate (\grad \phi).
%   derivMax	 Cell vector of maximum values of the costate (\grad \phi).
%   schemeData	 A structure (see below).
%   dim          Dimension in which the partial derivatives is taken.
%
%   alpha	 Maximum absolute value of the partial of the Hamiltonian
%		   with respect to the costate in dimension dim for the 
%                  specified range of costate values (O&F equation 5.12).
%		   Note that alpha can (and should) be evaluated separately
%		   at each node of the grid.
%
% schemeData is a structure containing data specific to this Hamiltonian
%   For this function it contains the field(s):
%
%   .grid	 Grid structure.
%   .velocityA	 Speed of the evader (positive constant).
%   .velocityB	 Speed of the pursuer (positive constant).
%   .inputA	 Maximum turn rate of the evader (positive).
%   .inputB	 Maximum turn rate of the pursuer (positive).
%
% Ian Mitchell 3/26/04

checkStructureFields(schemeData, 'grid', 'velocityA', 'velocityB');

grid = schemeData.grid;

switch dim
  case 1
    alpha = schemeData.velocityA;

  case 2
    alpha = schemeData.velocityA;

  case 3
    alpha = schemeData.velocityB;
    
  case 4
    alpha = schemeData.velocityB;

  otherwise
    error([ 'Partials for the sdsi game' ...
            ' only exist in dimensions 1-4' ]);
end
