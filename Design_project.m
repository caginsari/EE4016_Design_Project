close all
clear
%% EE4016 Design Project Resit.
%% AF calculation.
freq = 6e9 ;
lambda = 3e8/6e9 ;
N = 19 ;
k = 2*pi / lambda ;
dxx = lambda/4 ;
dyy = sqrt(3)*lambda/4 ;

% Scan Directions
th_0 = deg2rad(45) ;
phi_0 = deg2rad(45) ;

% Theta and Phi Meshgrid
theta = linspace(eps,pi/2,101);
phi = linspace(eps,2*pi,303);
[TH,PHI] = meshgrid(theta,phi) ;

% U-V Coordinate system
u = sin(TH) .* cos(PHI) ;
v = sin(TH) .* sin(PHI) ;

% scan directions in U-V coordinates
u0 = sin(th_0).* cos(phi_0) ;
v0 = sin(th_0).* sin(phi_0) ;

% Meshgrid in cartesian coordinates
% x changes in columns and y changes in rows.
x = -4*dxx:dxx:4*dxx ;
y = -2*dyy:dyy:2*dyy ;
[X,Y] = meshgrid(x,y) ;

fact(:,1) = [ 0 0 1 0 1 0 1 0 0 ] ;
fact(:,2) = [ 0 1 0 1 0 1 0 1 0 ] ;
fact(:,3) = [ 1 0 1 0 1 0 1 0 1 ] ;
fact(:,4) = [ 0 1 0 1 0 1 0 1 0 ] ;
fact(:,5) = [ 0 0 1 0 1 0 1 0 0 ] ;

a = squeeze(fact)';

Beta_n = -k .* ( X.* u0 + Y.*v0 ) ; 

%% AF for generic case.
for jj = 1:length(x)
    for ii = 1:length(y)

        if ii==1 && jj==1
            af_old = 0;
        else
            af_old = AF ;
        end

        af = a(ii,jj).*1./ N .* exp(1i .* Beta_n(ii,jj) ) .* exp(1i.*k.*( X(ii,jj).*u + Y(ii,jj).*v) ) ;

        AF = af + af_old ;
        
    end
end

%% UV surface representation
U = sin(TH) .* cos(PHI);
V = sin(TH) .* sin(PHI);
figure;
surface(U, V, 20*log10(abs(AF)) - max(max(20*log10(abs(AF)))), 'linestyle', 'none') % When the number of point increases, it would be better to hide lines
xlabel('U')
ylabel('V')
grid on
colormap('jet')
colorbar
caxis([-40, 0]) 
view(-37.5,30) % Change the view angle [Az,El]

GratingLobeDiagram(k,2*dxx,2*dyy,th_0,phi_0,'Triangular');