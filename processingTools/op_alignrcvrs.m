% op_alignrcvrs.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,coilcombos]=op_alignrcvrs(in,point,mode,coilcombos);
% 
% DESCRIPTION:
% phase align the receiver channels without combining them.
% 
% INPTUS:
% in            = input spectrum in matlab structure format.
% point         = Index of point in time domain to use for phase reference. 
%                 (optional.  Default = 1);  
% mode          = Method for estimating the coil weights and phases (optional.  Default = 'w').
%                 -'w' performs amplitude weighting of channels based on the
%                 maximum signal of each coil channel.
%                 -'h' performs amplitude weighting of channels based on the
%                 maximum signal of each coil channel divided by the square of
%                 the noise in each coil channel (as described by Hall et al.
%                 Neuroimage 2014). 
% coilcombos	= The predetermined coil phases (in [degrees])and amplitudes 
%                 as generated by the op_getcoilcombos.m function.  If this
%                 argument is provided, the 'point', and 'mode', arguments 
%                 will be ignored.  (Optional.  Default = []).
%
% OUTPUTS:
% out           = Output following alignment of rf channels.  
% coilcombos    = Structure containing two fields:
%                   ph:  Vector of coil phases (in [degrees]) used for alignment.
%                   sig: Vector of coil weights.

function [out,coilcombos]=op_alignrcvrs(in,point,mode,coilcombos);

if in.flags.addedrcvrs
    error('ERROR:  Receivers have already been combined!  Aborting!');
end

%To get best possible SNR, add the averages together (if it hasn't already been done):
if in.dims.averages>0
    av=op_averaging(in);
else
    av=in;
end

%also, for best results, we will combine all subspectra:
if nargin<4
    if in.flags.isFourSteps
        av=op_fourStepCombine(av);
    end
    if in.dims.subSpecs>0
        av=op_combinesubspecs(av,'summ');
    end
    if nargin < 3
        mode='w';
        if nargin < 2
            point=1;
        end
    end
   
end
avfids=av.fids;
avspecs=av.specs;

%initialize phase matrix and the amplitude maxtrix that are the size of nPoints x Coils
ph=ones(in.sz(in.dims.t),in.sz(in.dims.coils));
sig=ones(in.sz(in.dims.t),in.sz(in.dims.coils));

if nargin<4
    coilcombos.ph=zeros(in.sz(in.dims.coils),1);
    coilcombos.sig=zeros(in.sz(in.dims.coils),1);
end

%now start finding the relative phases between the channels and populate
%the ph matrix
for n=1:in.sz(in.dims.coils)
    if nargin<4
        p=phase(avfids(point,n,1,1))*180/pi; %in [degrees]
        coilcombos.ph(n)=p;
        ph(:,n)=p*ph(:,n);
        switch mode
            case 'w'
                S=abs(avfids(point,n,1,1));
                coilcombos.sig(n)=S;
                sig(:,n)=S*sig(:,n);
            case 'h'
                S=abs(avfids(point,n,1,1));
                N=std(avfids(end-100:end,n,1,1));
                sig(:,n)=(S/(N.^2))*sig(:,n);
                coilcombos.sig(n)=(S/(N.^2));
        end
    else
        ph(:,n)=coilcombos.ph(n)*ph(:,n);  %in [degreees]
        sig(:,n)=coilcombos.sig(n)*sig(:,n);
    end
end

%now replicate the phase matrix to equal the size of the original matrix:
replicate=in.sz;
replicate(1)=1;
replicate(2)=1;
ph=repmat(ph,replicate);
%sig=repmat(sig,replicate);
%sig=sig/max(max(max(max(sig))));


%now apply the phases by multiplying the data by exp(-i*ph);
fids=in.fids.*exp(-i*ph*pi/180);
fids_presum=fids;
specs_presum=fftshift(ifft(fids,[],in.dims.t),in.dims.t);


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids_presum;
out.specs=specs_presum;



