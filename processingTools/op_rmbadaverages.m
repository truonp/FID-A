%op_rmbadaverages.m
%Jamie Near, McGill University 2014.
%
%USAGE:
%[out,metric,badAverages]=op_rmbadaverages(in,nsd,domain);
%
%DESCRIPTION:
%Removes motion corrupted averages from a dataset containing multiple
%averages.  Bad averages are identified by calculating a 'likeness' metric
%for each average.  This is done by subtracting each average from the
%average of the averages, and then calculating the root mean squared of
%this difference spectrum.  Averages whose likeness metrics are greater
%than 'nsd' above the mean are discarded.
%
%INPUTS:
%in         = input data in matlab structure format
%nsd        = number of standard deviations to use a rejection threshold
%domain     = domain in which to perform calculations ('t' or 'f')

function [out,metric,badAverages]=op_rmbadaverages(in,nsd,domain);

if in.flags.averaged
    error('ERROR:  Averaging has already been performed!  Aborting!');
end

if ~in.flags.addedrcvrs
    error('ERROR:  Receivers should be combined first!  Aborting!');
end

%first, make a metric by subtracting all averages from the first average, 
%and then taking the sum of all all the spectral points.  
if in.dims.subSpecs>0
    SS=in.sz(in.dims.subSpecs);
else
    SS=1;
end
if domain=='t' || domain=='T'
    infilt=in;
    tmax=0.4;
elseif domain=='f' || domain=='F'
    filt=10;
    infilt=op_filter(in,filt);
end
inavg=op_averaging(infilt);
for n=1:in.sz(in.dims.averages)
    for m=1:SS
        if domain=='t' || domain=='T'
            metric(n,m)=sum((real(infilt.fids(infilt.t>=0 & infilt.t<=tmax,n,m))-(real(inavg.fids(inavg.t>=0 & inavg.t<=tmax,m))/in.sz(in.dims.averages))).^2);
        elseif domain=='f' || domain=='F'
            metric(n,m)=sum((real(infilt.specs(:,n,m))-(real(inavg.specs(:,m))/in.sz(in.dims.averages))).^2);
        end
    end
end

%find the average and standard deviation of the metric
avg=mean(metric);
stdev=std(metric);

for m=1:SS
    
    P(m,:)=polyfit([1:in.sz(in.dims.averages)]',metric(:,m),2);
    figure('position',[0 (m-1)*500 560 420]);
    plot([1:in.sz(in.dims.averages)],metric(:,m),'.',...
        [1:in.sz(in.dims.averages)],polyval(P(m,:),[1:in.sz(in.dims.averages)]),...
        [1:in.sz(in.dims.averages)],(polyval(P(m,:),[1:in.sz(in.dims.averages)])'+(nsd*stdev(m))),':');
    xlabel('Scan Number');
    ylabel('Unlikeness Metric');
    title('Metric for rejection of motion corrupted scans');
end

%Now make a mask that represents the locations of the averages 
%whose metric values are more than nsd standard deviations away from the 
%mean metric value.

for n=1:SS
    %mask(:,n)=metric(:,n)>(avg(n)+(nsd*stdev(n))) | metric(:,n)<(avg(n)-(nsd*stdev(n)));
    %mask(:,n)=metric(:,n)>(polyval(P(n,:),[1:in.sz(in.dims.averages)])'+(nsd*stdev(n))) | metric(:,n)<(polyval(P(n,:),[1:in.sz(in.dims.averages)])'-(nsd*stdev(n)));
    mask(:,n)=metric(:,n)>(polyval(P(n,:),[1:in.sz(in.dims.averages)])'+(nsd*stdev(n)));
end

%Unfortunately, if one average is corrupted, then all of the subspecs
%corresponding to that average have to be thrown away.  Therefore, take the
%minimum intensity projection along the subspecs dimension to find out
%which averages contain at least one corrupted subspec:
if size(mask,2)>1
    mask=sum(mask')'>0;
end

%now the corrupted and uncorrupted average numbers are given by:
badAverages=find(mask);
goodAverages=find(~mask);

%make a new fids array containing only good averages
fids=in.fids(:,goodAverages,:,:);

%%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%re-calculate the sz variable
sz=size(fids);

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.averages=length(goodAverages) * in.rawSubspecs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
