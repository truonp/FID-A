%op_creFit.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% parsFit=op_creFit(in,ph0,ph1);
% 
% DESCRIPTION:
% Perform a Lorentzian lineshape fit to the creatine resonance in a brain 
% proton MRS dataset.
% 
% INPUTS:
% in	= input data in matlab stucture format.
% ph0	= zero order phase to add to the input data.
% ph1	= 1st order phase to add to the input data.

function parsFit=op_creFit(in,ph0,ph1);

if in.flags.isISIS
    error('ERROR:  must have combined subspecs in order to do this!  ABORTING');
end

if ~in.flags.averaged
    error('ERROR:  must have averaged in order to do this!  ABORTING');
end

if ~in.flags.addedrcvrs
    error('ERROR:  must have added receivers in order to do this!  ABORTING');
end
specs=addphase(addphase1(in.specs,ph1),ph0);

ppm=in.ppm((in.ppm>4.6)&(in.ppm<4.85));
spec=specs(((in.ppm>4.6)&(in.ppm<4.85)));

plot(ppm,spec);


parsGuess=zeros(1,4);
parsGuess(1)=max(real(spec)); %Amplitude
parsGuess(2)=3/(42.577*in.Bo); %Linewidth [Hz]
parsGuess(3)= 3; %Frequency [ppm];
parsGuess(4)=(real(spec(1))-real(spec(end)))/(ppm(1)-ppm(end));  %Baseline slope
parsGuess(5)=real(spec(1))-(parsGuess(:,4)*ppm(1));   %Baseline Offset

yGuess=op_lorentz_linbas(parsGuess,ppm);
parsFit=nlinfit(ppm,real(spec'),@op_lorentz_linbas,parsGuess);
yFit=op_lorentz_linbas(parsFit,ppm);

plot(ppm,spec,ppm,yGuess,':',ppm,yFit);
legend('data','guess','fit');
parsFit
parsGuess

area=parsFit(:,1).*parsFit(:,2);
for n=1:size(area,1)
    disp(['Area under the ' num2str(n) 'th fitted curve is: ' num2str(area(n))]);
end
area=sum(area);
disp(['Area under the fitted curve is: ' num2str(area)]);


