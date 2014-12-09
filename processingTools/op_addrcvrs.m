%op_addrcvrs.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [out,fids_presum,specs_presum,ph,sig]=op_addrcvrs(in,point,mode,coilcombos);
% 
% DESCRIPTION:
% Perform weighted coil recombination for MRS data acquired with a reciever
% coil array.
% 
% INPUTS:
% in            = input spectrum in matlab structure format.
% point         = point of fid to use for phase estimation.
% mode          = -'w' performs amplitude weighting of channels based on the
%                 maximum signal of each coil channel.
%                 -'h' performs amplitude weighting of channels based on the
%                 maximum signal of each coil channel divided by the square of
%                 the noise in each coil channel (as described by Hall et al.
%                 Neuroimage 2014). 
% coilcombos    = (optional)  The predetermined coil phases and amplitudes as
%                 generated by the op_getcoilcombos.m function.  If this
%                 argument is provided, the 'point', and 'mode', arguments
%                 will be ignored.
         
function [out,fids_presum,specs_presum,ph,sig]=op_addrcvrs(in,point,mode,coilcombos);

if in.flags.addedrcvrs
    error('ERROR:  Receivers have already been combined!  Aborting!');
end

%To get best possible SNR, add the averages together (if it hasn't already been done):
if ~in.flags.averaged
    av=op_averaging(in);
else
    av=in;
end

%also, for best results, we will combine all subspectra:
if nargin<4
    if in.flags.isISIS
        av=op_ISIScombine(av);
    end
    if in.dims.subSpecs>0
        av=op_combinesubspecs(av,'summ');
    end
end
avfids=av.fids;
avspecs=av.specs;

%initialize phase matrix and the amplitude maxtrix that are the size of nPoints x Coils
ph=ones(in.sz(in.dims.t),in.sz(in.dims.coils));
sig=ones(in.sz(in.dims.t),in.sz(in.dims.coils));

%now start finding the relative phases between the channels and populate
%the ph matrix
for n=1:in.sz(in.dims.coils)
    if nargin<4
        ph(:,n)=phase(avfids(point,n,1,1))*ph(:,n);
        switch mode
            case 'w'
                sig(:,n)=abs(avfids(point,n,1,1))*sig(:,n);
            case 'h'
                S=max(abs(avfids(:,n,1,1)));
                N=std(avfids(end-100:end,n,1,1));
                sig(:,n)=(S/(N.^2))*sig(:,n);
        end
    else
        ph(:,n)=coilcombos.ph(n)*ph(:,n);
        sig(:,n)=coilcombos.sig(n)*sig(:,n);
    end
end

%now replicate the phase matrix to equal the size of the original matrix:
replicate=in.sz;
replicate(1)=1;
replicate(2)=1;
ph=repmat(ph,replicate);
sig=repmat(sig,replicate);
sig=sig/max(max(max(max(sig))));


%now apply the phases by multiplying the data by exp(-i*ph);
fids=in.fids.*exp(-i*ph);
fids_presum=fids;
specs_presum=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%Apply the amplitude factors by multiplying the data by amp;
if mode=='w' || mode=='h'
    fids=fids.*sig;
end


%now sum along coils dimension
fids=sum(fids,in.dims.coils);
fids=squeeze(fids);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%change the dims variables
dims.t=in.dims.t;
dims.coils=0;
dims.averages=in.dims.averages-1+((in.dims.averages-1)<0); %Don't let it go negative
dims.subSpecs=in.dims.subSpecs-1+((in.dims.subSpecs-1)<0); %Don't let it go negative

%re-calculate the sz variable
sz=size(fids);

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.dims=dims;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.addedrcvrs=1;



