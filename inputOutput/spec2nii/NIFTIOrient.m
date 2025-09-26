function out=NIFTIOrient(affine)
out.Q44 = affine;
[qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac] = nifti_mat44_to_quatern(affine);
out.qb=qb;
out.qc=qc;
out.qd=qd;
out.qx=qx;
out.qy=qy;
out.qz=qz;
out.dx=dx;
out.dy=dy;
out.dz=dz;
out.qfac=qfac;
end

function [qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac]=nifti_mat44_to_quatern(R)
% """4x4 affine to quaternion representation."""
% offset outputs are read out of input matrix
qx=R(1,4);
qy=R(2,4);
qz=R(3,4);

% load 3x3 matrix into local variables
r11 = R(1, 1);
r12 = R(1, 2);
r13 = R(1, 3);
r21 = R(2, 1);
r22 = R(2, 2);
r23 = R(2, 3);
r31 = R(3, 1);
r32 = R(3, 2);
r33 = R(3, 3);

% compute lengths of each column; these determine grid spacings
xd = sqrt(r11 * r11 + r21 * r21 + r31 * r31);
yd = sqrt(r12 * r12 + r22 * r22 + r32 * r32);
zd = sqrt(r13 * r13 + r23 * r23 + r33 * r33);

% if a column length is zero, path the trouble
if xd == 0.0
    r11 = 1.0;
    r21 = 0.0;
    r31 = 0.0;
    xd = 1.0;
end
if yd == 0.0
    r22 = 1.0;
    r12 = 0.0;
    r32 = 0.0;
    yd = 1.0;
end
if zd == 0.0
    r33 = 1.0;
    r13 = 0.0;
    r23 = 0.0;
    zd = 1.0;
end

% assign the output lengths
dx = xd;
dy = yd;
dz = zd;

% normalize the columns
r11 = r11 / xd;
r21 = r21 / xd;
r31 = r31 / xd;
r12 = r12 / yd;
r22 = r22 / yd;
r32 = r32 / yd;
r13 = r13 / zd;
r23 = r23 / zd;
r33 = r33 / zd;

zd = r11 * r22 * r33 ...
    - r11 * r32 * r23 ...
    - r21 * r12 * r33 ...
    + r21 * r32 * r13 ...
    + r31 * r12 * r23 ...
    - r31 * r22 * r13;
% zd should be -1 or 1

if zd > 0  % proper
    qfac = 1.0;
else  % improper ==> flip 3rd column
    qfac = -1.0;
    r13 = r13 * -1.0;
    r23 = r23 * -1.0;
    r33 = r33 * -1.0;
end

% now, compute quaternion parameters
a = r11 + r22 + r33 + 1.0;
if a > 0.5  % simplest case
    a = 0.5 * sqrt(a);
    b = 0.25 * (r32 - r23) / a;
    c = 0.25 * (r13 - r31) / a;
    d = 0.25 * (r21 - r12) / a;
else  % trickier case
    xd = 1.0 + r11 - (r22 + r33);  % 4*b*b
    yd = 1.0 + r22 - (r11 + r33);  % 4*c*c
    zd = 1.0 + r33 - (r11 + r22);  % 4*d*d
    if xd > 1.0
        b = 0.5 * sqrt(xd);
        c = 0.25 * (r12 + r21) / b;
        d = 0.25 * (r13 + r31) / b;
        a = 0.25 * (r32 - r23) / b;
    elseif yd > 1.0
        c = 0.5 * sqrt(yd);
        b = 0.25 * (r12 + r21) / c;
        d = 0.25 * (r23 + r32) / c;
        a = 0.25 * (r13 - r31) / c;
    else
        d = 0.5 * sqrt(zd);
        b = 0.25 * (r13 + r31) / d;
        c = 0.25 * (r23 + r32) / d;
        a = 0.25 * (r21 - r12) / d;
    end

    if a < 0.0
        b = -b;
        c = -c;
        d = -d;
    end
end

qb = b;
qc = c;
qd = d;
end