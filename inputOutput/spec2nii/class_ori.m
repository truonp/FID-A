function case_val=class_ori(sag_comp,cor_comp,tra_comp,debug)
% ''' MATLAB implementation of IDEA-VB17/n4/pkg/MrServers/MrMeasSrv/SeqFW/libGSL/fGSLClassOri.cpp
% Function to determine whether a normal vector describes a sagittal, coronal or transverse slice.
% Result:
%     CASE = 0: Sagittal
%     CASE = 1: Coronal
%     CASE = 2: Transverse
%
% :param  sag_comp:   Sagittal component of normal vector
% :param  cor_comp:   Coronal component of normal vector
% :param  tra_comp:   Transverse component of normal vector
%
% :return:    case (0=Sagittal, 1=Coronal or 2=Transverse)

if nargin<4
    debug=0;
end

% Compute some temporary values
abs_sag_comp = abs(sag_comp);
abs_cor_comp = abs(cor_comp);
abs_tra_comp = abs(tra_comp);

eq_sag_cor = abs(abs_sag_comp - abs_cor_comp) < eps;
eq_sag_tra = abs(abs_sag_comp - abs_tra_comp) < eps;
eq_cor_tra = abs(abs_cor_comp - abs_tra_comp) < eps;

if ((eq_sag_cor && eq_sag_tra) || ...
    (eq_sag_cor && (abs_sag_comp < abs_tra_comp)) || ...
    (eq_sag_tra && (abs_sag_comp > abs_cor_comp)) || ...
    (eq_cor_tra && (abs_cor_comp > abs_sag_comp)) || ...
    ((abs_sag_comp > abs_cor_comp) && (abs_sag_comp < abs_tra_comp)) || ...
    ((abs_sag_comp < abs_cor_comp) && (abs_cor_comp < abs_tra_comp)) || ...
    ((abs_sag_comp < abs_tra_comp) && (abs_tra_comp > abs_cor_comp)) || ...
    ((abs_cor_comp < abs_tra_comp) && (abs_tra_comp > abs_sag_comp)))

    if debug
        disp('Mainly transverse.');
    end
    case_val = 2; % Transverse

elseif ((eq_sag_cor && (abs_sag_comp > abs_tra_comp)) || ...
        (eq_sag_tra && (abs_sag_comp < abs_cor_comp)) || ...
        ((abs_sag_comp < abs_cor_comp) && (abs_cor_comp > abs_tra_comp)) || ...
        ((abs_sag_comp > abs_tra_comp) && (abs_sag_comp < abs_cor_comp)) || ...
        ((abs_sag_comp < abs_tra_comp) && (abs_tra_comp < abs_cor_comp)))

    if debug
        disp('Mainly coronal.');
    end
    case_val = 1; % Coronal

elseif ((eq_cor_tra && (abs_cor_comp < abs_sag_comp)) || ...
        ((abs_sag_comp > abs_cor_comp) && (abs_sag_comp > abs_tra_comp)) || ...
        ((abs_cor_comp > abs_tra_comp) && (abs_cor_comp < abs_sag_comp)) || ...
        ((abs_cor_comp < abs_tra_comp) && (abs_tra_comp < abs_sag_comp)))

    if debug
        disp('Mainly sagittal.');
    end
    case_val = 0; % Sagittal

else
    error('Error: Invalid slice orientation');
end

return_val = case_val;
end