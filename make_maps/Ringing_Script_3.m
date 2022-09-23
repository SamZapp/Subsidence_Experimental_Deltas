function [filtmap] = Ringing_Script_3(gauss_filt_map);
subs=gauss_filt_map;

cnt=1:length(subs(:));
[cntnan]= cnt(isnan(subs(:))==0);

[X,Y]=meshgrid(1:747,1:750);
subscoords=[X(isnan(subs(:))==0) Y(isnan(subs(:))==0) subs(isnan(subs(:))==0) zeros(length(X(isnan(subs(:))==0)),2)];

centerx=840;%Centerx centery is a guess.
centery=811;
%Now we are going to turn the measurements into radial coordinates
%from centerx,centery
vec=[subscoords(:,1)-centerx subscoords(:,2)-centery];
radial=sqrt(vec(:,1).^2+vec(:,2).^2);
theta=atan2(vec(:,2),vec(:,1));

searchradial=10; %search five pixel widths in the radial direction  5
searchtheta=1; %radians.   .15

for i=1:length(subscoords);
    neighborhood=subscoords(and(abs(radial-radial(i))<searchradial,abs(theta-theta(i))<searchtheta),3);
    subscoords(i,4)=median(neighborhood(neighborhood<2 & neighborhood>-2));
end

medfilt=ones(size(cnt))*nan;
medfilt(cntnan)=subscoords(:,4);
medfiltgrid = reshape(medfilt,size(subs));
filtmap=subs-medfiltgrid;