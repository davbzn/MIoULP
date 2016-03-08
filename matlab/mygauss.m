function y = mygauss(x,setFWHM,power)
%%mygauss maks a Gaussian with given FWHM and power (2,4,6,8,...)
% y = mygauss(x,setFWHM,power)
%
% now FWHM is correct for all powers
% if you want a given sigma enter FWHM = 2*(2*log(2))^(1/power) * sigma
%
% example:
% clc
% x = makecol(-200:1:200);
% power=2;
% setFWHM=100;
% y=mygauss(x,setFWHM,power);
% out = [ myFWHM(x,y), mymoment(x,y,2), std(y)];
% plot_debug(x,y); grid on
% title(sprintf('setFWHM = %.1f, FWHM = %.3f, sigma = %.3f, std = %.3f', [setFWHM,out] ) )
% 

if nargin<3;power = 2;end
% switch power
%     case 2
%         conversion = 1.665109222315396 ; % sqrt(4*log(2))
%         setsigma = setfwhm./sqrt(8*log(2));
%     case 4
%         conversion = 1.824888611568057 ;
%     case 6
%         conversion = 1.828955772123223 ;
%     otherwise
%         conversion = .5*(2^power*log(4)).^(1/power)*(2).^(1/power) ;%* log(power)*2 ;
%end
%y = exp( -(x./fwhm*conversion ).^power );

% changed: 02.05.2012: now given FWHM is correct for all powers:
setSigma = setFWHM ./ ( 2*(2*log(2))^(1/power) ); 
y = exp( -(abs(x).^power)./(2*setSigma.^power) );

end
