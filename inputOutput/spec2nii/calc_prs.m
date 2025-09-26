function [gp, gr] = calc_prs(gs, phi, debug)
% calc_prs - MATLAB implementation of IDEA-VB17/n4/pkg/MrServers/MrMeasSrv/SeqFW/libGSL/fGSLCalcPRS.cpp
% Calculates the phase encoding (gp) and readout (gr) direction vectors.
%
% Inputs:
%   gs   - Slice normal vector [3x1]
%   phi  - Rotational angle around gs (radians)
%   debug - Boolean flag for printing debug info
%
% Outputs:
%   gp - Phase direction vector [3x1]
%   gr - Read direction vector [3x1]
   
    % PCS axes
    SAGITTAL   = 0;
    CORONAL    = 1;
    TRANSVERSE = 2;

    if nargin<4
        debug=0;
    end
    % Start of function
    orientation = 0;  % will be one of SAGITTAL, CORONAL or TRANSVERSE
    orientation = orientation + class_ori(gs(1), gs(2), gs(3), debug);

    gp = zeros(3,1);

    if orientation == TRANSVERSE
        gp(1) = 0.0;
        gp(2) = gs(3) * sqrt(1.0 / (gs(2)^2 + gs(3)^2));
        gp(3) = -gs(2) * sqrt(1.0 / (gs(2)^2 + gs(3)^2));
    elseif orientation == CORONAL
        gp(1) = gs(2) * sqrt(1.0 / (gs(1)^2 + gs(2)^2));
        gp(2) = -gs(1) * sqrt(1.0 / (gs(1)^2 + gs(2)^2));
        gp(3) = 0.0;
    elseif orientation == SAGITTAL
        gp(1) = -gs(2) * sqrt(1.0 / (gs(1)^2 + gs(2)^2));
        gp(2) = gs(1) * sqrt(1.0 / (gs(1)^2 + gs(2)^2));
        gp(3) = 0.0;
    else
        error('Invalid slice orientation returned from class_ori');
    end

    % Calculate GR = GS x GP
    gr = cross(gs, gp);

    if debug
        fprintf('Before rotation around S:\n');
        fprintf('GP = %10.7f %10.7f %10.7f\n', gp(1), gp(2), gp(3));
        fprintf('GR = %10.7f %10.7f %10.7f\n', gr(1), gr(2), gr(3));
        fprintf('GS = %10.7f %10.7f %10.7f\n', gs(1), gs(2), gs(3));
    end

    % Rotation around GS axis
    if phi ~= 0.0
        if debug
            tmp = phi * 180.0 / pi;
            fprintf('PHI = %10.7f\n', tmp);
        end
    end
    gp(1) = cos(phi) * gp(1) - sin(phi) * gr(1);
    gp(2) = cos(phi) * gp(2) - sin(phi) * gr(2);
    gp(3) = cos(phi) * gp(3) - sin(phi) * gr(3);

    % Recompute GR = GS x GP
    gr = cross(gs, gp);

    if debug
        fprintf('After the Rotation around S:\n');
        fprintf('GP = %10.7f %10.7f %10.7f\n', gp(1), gp(2), gp(3));
        fprintf('GR = %10.7f %10.7f %10.7f\n', gr(1), gr(2), gr(3));
        fprintf('GS = %10.7f %10.7f %10.7f\n', gs(1), gs(2), gs(3));
    end
end