close all;
clear all;
clc;

%<-----------Définition des variables----------->
tmax = 10 ; % Durée 
g=10; % Constante gravitationelle
m=46*10.^-3; %Masse de la balle de golf
r=21.5*10.^-3; %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<mettre ce que c'est


% Constantes liées aux foces extérieures s'applicant sur la balle
env = struct('g',10,'alpha_t',4.5*10.^-5,'alpha_u',1.2*10.^-4,'alpha_s',5*10.^-7,'j',2/5*m*r.^2,'m',m);

%// Position à t=0
x0=0; 
y0=0;

%// Vecteur vitesse (t=0)
v0=52; % Norme de la vitesse en m/s
alpha0=32; % Angle du vecteur vitesse par rapportà l'horizontale 
vx0=v0*cosd(alpha0); % Vitesse selon l'axe x
vy0=v0*sind(alpha0);% Vitesse selon l'axe x

%// Vecteur rotation (t=0)
w0=155; 

%// Nombre de points où l'on définie la trajectoire
nb_time_steps = 1e3+1;


%<--------Calcul des vecteurs position, vitesse et spin en chaque point----->

%// Integration par rapport au temps des vitesses et accélérations obtenues
% grâce à la fonction calcul_vitesse_acceleration afin d'obtenir en tout
% instant t la position, les vitesses et le spin de la balle.  
[T,Z] = ode23( @(t,z) calcul_vitesse_acceleration(t,z,env), linspace (0, tmax,nb_time_steps), [x0; y0; vx0; vy0; w0]);

%Sélection des valeurs pour former la trajectoire décimée (les calculs se
%feront avec ces valeurs sélectionnées)
Z_dec = calcul_dec(T,Z);

%Calcul de la norme du vecteur vitesse
norm_v_dec =norm([Z_dec(:,3),Z_dec(:,4)]);


%<--------Calcul des forces modifiant la trajectoire de la balle-------->

%Couple résistant
couple_resistant_spin=(-env.alpha_s*Z_dec(:,5))';

%Trainée
trainee_dec = -env.alpha_t*norm_v_dec*[Z_dec(:,3),Z_dec(:,4)];

%Effet Magnus
z_cross_v=[-Z_dec(:,4),Z_dec(:,3)]; %produit vectoriel du vecteur unité z par la vitesse
magnus_dec=env.alpha_u*norm(Z_dec(:,5))*z_cross_v;






% <------------Question 1----------->
% On affiche la dernière position de la balle ainsi que sa vitesse et son spin
Aff_1_0=['Question 1) A tmax : '];
Aff_1_1=['Position :  x=', num2str(Z(nb_time_steps,1)),'m y = ',num2str(Z(nb_time_steps,2)),'m'];
Aff_1_2=['Vitesse : ' , num2str(norm([Z(nb_time_steps,3),Z(nb_time_steps,4)])),' m/s'];
Aff_1_3=['Spin : ',num2str(Z(nb_time_steps,5)),' rad/s'];
disp(Aff_1_0);
disp(Aff_1_1);
disp(Aff_1_2);
disp(Aff_1_3);


% <-----------Question 2---------->

%// Graphique de la trajectoire d'équation y en fonction de x
figure(1)
plot(Z_dec(:,1),Z_dec(:,2),'.m','Color','#EDB120','markerSize',15) 
title ('Trajectoire de la balle de golf ');
xlabel('x (m)')
ylabel('y (m)')
saveas(figure(1),'y=f(x).fig')


% <-----------Question 3---------->
[maxY_3,i_3]=max(Z(:,2)); %Recherche du point le plus haut parmis les points de la trajectoire non décimée
hold on 
% Affichage du point
plot(Z(i_3,1),Z(i_3,2),'+r','MarkerSize',10) 
Aff_3_0=['Question 3) L''altitude maximale est de : ', num2str(maxY_3),' m'];
disp(Aff_3_0);

% <-----------Question 4---------->
%Affichage des vecteurs vitesse
quiver(Z_dec(:,1),Z_dec(:,2),Z_dec(:,3),Z_dec(:,4),'r','AutoScaleFactor',0.4,'Color','#D95319'); 

% <-----------Question 5---------->
%Addition des puissances dissipée calculée grâce à la fonction calcul_puissance_dissip puis affichage
puissance_dissip_tot=calcul_puissance_dissip(trainee_dec,[Z_dec(:,3),Z_dec(:,4)])+ calcul_puissance_dissip(couple_resistant_spin,Z_dec(:,5));
Aff_5_0 =['Question 5) L''energie perdue entre T0 et Tmax est de ',num2str(puissance_dissip_tot),' joules'];
disp(Aff_5_0);

% <-----------Question 6---------->



%<------------ Mise en forme de la figure------------->
legend({'Trajectoire décimée','Altitude maximale','Vecteur vitesse'},'Location','northwest')

% <------Fonction donnant la vitesse, l'acceleration et l'acceleration angulaire ---->
function d =calcul_vitesse_acceleration(t,z,env)
v = [ z(3); z(4)] ; % Vecteur colonne Vx Vy
w =z(5) ; % Valeur du spin
vp = - [0 ;env.g]+1/env.m*(-env.alpha_t *norm(v)*v + env.alpha_u*w*[-v(2) ; v(1) ] ) ; % Calcul de l'acceleration grâce au TRD : gamma=Fext/m
wp = - 1/env.j * env.alpha_s * w ; % Calcul de l'accélération basée sur le TMD
d = [v;vp;wp];
end 

%<----------Fonction permettant l'intégration par la méthode desrectangles
function aire = calcul_puissance_dissip(force, vitesse)
    puissance_dissip_trainee = vitesse.*force; 
    aire=0;
    for i=1:length(puissance_dissip_trainee)-1
         moy = (norm(puissance_dissip_trainee(i,:),2)+norm(puissance_dissip_trainee(i+1,:),2))/2;
         aire=aire+moy*0.5;
    end  
end 

%<-------- Fonction calculant la distance totale parcourue par la balle--->
function d = calcul_distance_traj(traj_fournie,traj_calc)
    for i=1:length(traj_fournie)
        d=d+(traj_fournie(i,1)-traj_calc(i,1))^2+(traj_fournie(i,2)-traj_calc(i,2))^2;
    end
end

%<------- Fonction calculant une trajectoire à partir de valeurs initiales--->
function [x,y] = calcul_trajectoire(traj_fournie, valeur_init)
    [A,B] = ode23( @(t,z) calcul_vitesse_acceleration(t,z,env), linspace (0, tmax,nb_time_steps), [traj_fournie(1,1); traj_fournie(1,2); valeur_init(:,1); valeur_init(:,2); valeur_init(:,3)]);
    x=B(:,1);
    y=B(:,2);
end 

%<------- Fonction sélectionnant les éléments de la trajctoire de la balle tous les 0.5s--->
function Traj_dec = calcul_dec(Time, Traj_a_decimer)
    j=1;
    for i=1:length(Time)
        if mod(Time(i),0.5)==0
            Traj_dec(j,1)=Traj_a_decimer(i,1);
            Traj_dec(j,2)=Traj_a_decimer(i,2);
            Traj_dec(j,3)=Traj_a_decimer(i,3);
            Traj_dec(j,4)=Traj_a_decimer(i,4);
            Traj_dec(j,5)=Traj_a_decimer(i,5);
           j=j+1;
        end
    end
end
