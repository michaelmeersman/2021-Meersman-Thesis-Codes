clear; close; clc

lambda = 65;

alpha_plane = 6;

lam = [sind(lambda) cosd(lambda) 1/sind(lambda) 1/cosd(lambda) tand(lambda)];

alph = [sind(alpha_plane) cosd(alpha_plane) 1/sind(alpha_plane) 1/cosd(alpha_plane) tand(alpha_plane)];

k=1;
for i = 1:5
    for j = 1:5
        alpha_wing(k) = atand(lam(i)*alph(j));
        k = k+1;
    end
end

% scatter(1:25,alpha_wing)


A1 = alpha_wing(10);

% alpha_wing(22)


% this gives the correct result 
A2 = asind(sind(alpha_plane)*cosd(lambda)/sqrt(1-sind(lambda)^2*sind(alpha_plane)^2));


A3 = acosd(cosd(alpha_plane)/cosd(lambda)/sqrt(tand(lambda)^2*cosd(alpha_plane)^2+1));


theta = atand(tand(lambda)*cosd(alpha_plane));

A4 = atand(tand(alpha_plane)/cosd(lambda));