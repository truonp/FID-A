% op_addrcvrs.m
% Jamie Near, McGill University 2014.
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
% point         = point of fid to use for phase estimation (optional. Default = 1);
% mode          = Method for estimating the coil weights and phases (optional.  Default = 'w').
%                 -'w' performs amplitude weighting of channels based on the
%                 maximum signal of each coil channel.
%                 -'h' performs amplitude weighting of channels based on the
%                 maximum signal of each coil channel divided by the square of
%                 the noise in each coil channel (as described by Hall et al.
%                 Neuroimage 2014). 
% coilcombos    = The predetermined coil weights and phases and amplitudes as
%                 generated by the op_getcoilcombos.m function.  If this
%                 argument is provided, the 'point', and 'mode', arguments
%                 will be ignored.  (optional.  Default = []).  
%
% OUTPUTS:
% out           = Output dataset with coil channels combined.
% fids_presum   = Input data with coil channels in phase (time domain).
% specs_presum  = Input data with coil channels in phase (frequency domain).
% ph            = Vector of applied coil phases (in degrees).
% sig           = Vector of coil weights.

         
function [out,fids_presum,specs_presum,ph,sig]=op_addrcvrs(in,point,mode,coilcombos);

if in.flags.addedrcvrs || ~in.dims.coils
    disp('WARNING:  Only one receiver channel found!  Returning input without modification!');
    out=in;
    out.flags.addedrcvrs=1;
    fids_presum=in.fids;
    specs_presum=in.specs;
    ph=0;
    sig=1;
else
    
    %To get best possible SNR, add the averages together (if it hasn't already been done):
    if ~in.flags.averaged
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
        if nargin<3
            mode = 'w';
            if nargin<2
                point=1;
            end
        end
    end
    avfids=av.fids;
    avspecs=av.specs;
    
    %initialize phase matrix and the amplitude maxtrix that are the size of nPoints x Coils
    %ph=ones(in.sz(in.dims.t),in.sz(in.dims.coils));
    %sig=ones(in.sz(in.dims.t),in.sz(in.dims.coils));
    
    %now start finding the relative phases between the channels and populate
    %the ph matrix
    for n=1:in.sz(in.dims.coils)
        if nargin<4
            %ph(:,n)=phase(avfids(point,n,1,1))*ph(:,n);
            phs(n)=phase(avfids(point,n,1,1));
            switch mode
                case 'w'
                    %sig(:,n)=abs(avfids(point,n,1,1))*sig(:,n);
                    sigs(n)=abs(avfids(point,n,1,1));
                case 'h'
                    S=max(abs(avfids(:,n,1,1)));
                    N=std(avfids(end-100:end,n,1,1));
                    %sig(:,n)=(S/(N.^2))*sig(:,n);
                    sigs(n)=(S/(N.^2));
            end
        else
            %ph(:,n)=coilcombos.ph(n)*ph(:,n);
            phs(n)=coilcombos.ph(n);
            sigs(n)=coilcombos.sig(n);
        end
    end
    
    %now replicate the phase matrix to equal the size of the original matrix:
    % replicate=in.sz;
    % replicate(1)=1;
    % replicate(in.dims.coils)=1;
    % ph=repmat(ph,replicate);
    % sig=repmat(sig,replicate);
    sigs=sigs/norm(sigs(:));
    
    ph=ones(in.sz);
    sig=ones(in.sz);
    
    if in.dims.coils==1
        for n=1:in.sz(1)
            ph(n,:)=phs(n)*ph(n,:);
            sig(n,:)=sigs(n)*sig(n,:);
        end
    elseif in.dims.coils==2
        for n=1:in.sz(2)
            ph(:,n,:)=phs(n)*ph(:,n,:);
            sig(:,n,:)=sigs(n)*sig(:,n,:);
        end
    elseif in.dims.coils==3
        for n=1:in.sz(3)
            ph(:,:,n,:)=phs(n)*ph(:,:,n,:);
            sig(:,:,n,:)=sigs(n)*sig(:,:,n,:);
        end
    elseif in.dims.coils==4
        for n=1:in.sz(4)
            ph(:,:,:,n,:)=phs(n)*ph(:,:,:,n,:);
            sig(:,:,:,n,:)=sigs(n)*sig(:,:,:,n,:);
        end
    elseif in.dims.coils==5
        for n=1:in.sz(5)
            ph(:,:,:,:,n)=phs(n)*ph(:,:,:,:,n);
            sig(:,:,:,:,n)=sigs(n)*sig(:,:,:,:,n);
        end
    end
    
    
    %now apply the phases by multiplying the data by exp(-i*ph);
    fids=in.fids.*exp(-i*ph);
    fids_presum=fids;
%     specs_presum=fftshift(ifft(fids,[],in.dims.t),in.dims.t);
    specs_presum=FIDAfft(fids,in.dims.t,'t');
    
    %Apply the amplitude factors by multiplying the data by amp;
    if mode=='w' || mode=='h'
        fids=fids.*sig;
    end
    
    
    %now sum along coils dimension
    fids=sum(fids,in.dims.coils);
    fids=squeeze(fids);
    
    %re-calculate Specs using fft
%     specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);
    specs=FIDAfft(fids,in.dims.t,'t');
    
    %change the dims variables
    if in.dims.t>in.dims.coils
        dims.t=in.dims.t-1;
    else
        dims.t=in.dims.t;
    end
    dims.coils=0;
    if in.dims.averages>in.dims.coils
        dims.averages=in.dims.averages-1;
    else
        dims.averages=in.dims.averages;
    end
    if in.dims.subSpecs>in.dims.coils
        dims.subSpecs=in.dims.subSpecs-1;
    else
        dims.subSpecs=in.dims.subSpecs;
    end
    if in.dims.extras>in.dims.coils
        dims.extras=in.dims.extras-1;
    else
        dims.extras=in.dims.extras;
    end
    
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
    
end

end




