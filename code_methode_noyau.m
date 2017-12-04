% ========================================================================
%CALCUL DU TAUX D'ERREUR BINAIRE
% Script de simulation
% .........................METHODE DU NOYAU
% ========================================================================

clear all;
close all;

% ----------------------------------------------------------------------------
% Définitions des Variables
% ----------------------------------------------------------------------------
% _Cp = Ensemble des sorties souples positifs à l'origine
% _Cn = Ensemble des Sorties souples négatifs à l'origine
% _Np = card(_C+) 
% _Nn = card(_C-)
% _N = Nn+Np
% _Sigmap = ecart type de _C+
% _Sigman = ecart type de _C-
% _Lissagep = paramètre de lissage pour les _C+
% _Lissagen = paramètre de lissage pour les _C-
% _Pip = probabilité liée à _C+
% _Pin = probabilité liée à _C-
% _T = le nombre d'itération
%_rop = proba à posteriori _C+
%_ron = proba à posteriori _C-

% ----------------------------------------------------------------------------
% TRANSMISSION
% ----------------------------------------------------------------------------

N=10^3; % Taille du message
d=round(rand(1,N)); % génération pseudo-aléatoire des données sources 
b=2.*d-1; % BPSK modulation 0 -C -1; 1 -> 0

Eb_N0_dB = 0:1; % Plage de simulation Eb/N0 

T=10;     %Nombre d'itération, le choix du nombre d'itération Doit se faire
          %de tel sorte que
          % à terme la valeur entre deux derniers paramètres de lissage
          % optimaux soient très petit ex=|h(t)-h(t-1)| inférieure à
          % 10^(-3)
          
for ii = 1:length(Eb_N0_dB)    
    
    disp(sprintf('Simu pour Eb/No = %.2f dB', Eb_N0_dB(ii)));

        % ----------------------------------------------------------------------------
        % 1 INITIALISATION
        % ----------------------------------------------------------------------------
        EbNo_lin = 10^(Eb_N0_dB(ii)/10);
        sigma2 = 1/(2*EbNo_lin);
        
        X = b + (sqrt(sigma2)*randn(1,N)); % passage dans le canal, application du Bruit Blanc Gaussien 
        jjp=1;
        jjn=1;
        Cp=0;
        Cn=0;
        %1.1 Classification
        for j=1:length(X)
            if X(j)>= 0
                Cp(jjp)=X(j);
                jjp=jjp+1;
            else
                Cn(jjn)=X(j);
                jjn=jjn+1;   
            end
        end
        
        %1.2 Calculs des N+ et N-
        Np=length(Cp);
        Nn=length(Cn);
        
        %1.3 Calcul de sigma, l'ecart type
        sigmap=sqrt(var(Cp));
        sigman=sqrt(var(Cn));
        
        %1.4 calcul du paramètre de lissage optimal
        Lissagep = ((4/(3*Np))^(1/5))*sigmap;
        Lissagen = ((4/(3*Nn))^(1/5))*sigman;
        
        %1.5 calcul des pi, probabilités à priori
        pip = Np/N;
        pin = Nn/N;
        
        %1.6 calcul des fonctions de densité conditionnelles       
        for j=1:length(X)
            for jj=1:length(Cp)
            fp(j,jj) = (1/2)*erfc(((X(j)-Cp(jj))/Lissagep)/sqrt(2));       
            end       
        end
        fp=(1/(Np*Lissagep)).*sum(fp');
               
        for j=1:length(X)
            for jj=1:length(Cn)
            fn(j,jj) = (1/2)*erfc(((X(j)-Cn(jj))/Lissagen)/sqrt(2));        
            end
        end
        fn=(1/(Nn*Lissagen)).*sum(fn');
        
        for t = 2 : 4
            % ----------------------------------------------------------------------------
            % 2 ESTIMATION DES PARAMETRES par itération
            % ----------------------------------------------------------------------------
            
            %2.1 Estimation Step, estimation des probabilités à postériori
            for j=1:length(fp) 
                rop(j)=(pip*fp(j))/(pip*fp(j)+pin*fn(j));
            end
            
            for j=1:length(fn) 
                ron(j)=(pin*fn(j))/(pip*fp(j)+pin*fn(j));
            end
            
            %2.2 étape de maximisation, calcul des pi
            pip=(1/N)*sum(rop);
            pin=(1/N)*sum(ron);
            
            %2.3 Etape de classification
            
            %2.3.1 Classification en utilisant les proba a posteriori
            jjp=1;
            jjn=1;
            Cp=0;
            Cn=0;
            U=rand(1,N); %Distribution uniforme sur 0-1; servant pour la classification

            for j=1:length(X)
                if rop(j)>= U(j)
                    Cp(jjp)=X(j);
                    jjp=jjp+1;
                else
                    Cn(jjn)=X(j);
                    jjn=jjn+1;   
                end
            end
            
            %2.3.2 calcul de la taille des deux classes
            Np=length(Cp);
            Nn=length(Cn);
            
            %2.3.3 calcul de l'ecart type de la distribution pour chaque
            %classe
            sigmap=sqrt(var(Cp));
            sigman=sqrt(var(Cn));
            
            %2.3.4 calcul du paramètre de lissage optimal
            Lissagep = ((4/(3*Np))^(1/5))*sigmap;
            Lissagen = ((4/(3*Nn))^(1/5))*sigman;
            
            %2.3.5 calcul des fonctions de densité conditionnelles
            for j=1:length(X)
                for jj=1:length(Cp)
                fp(j,jj) = (1/2)*erfc(((X(j)-Cp(jj))/Lissagep)/sqrt(2));       
                end       
            end
            fp=(1/(Np*Lissagep)).*sum(fp');
               
            for j=1:length(X)
                for jj=1:length(Cn)
                fn(j,jj) = (1/2)*erfc(((X(j)-Cn(jj))/Lissagen)/sqrt(2));        
                end
            end
            fn=(1/(Nn*Lissagen)).*sum(fn');

            
        end
        
        % ----------------------------------------------------------------------------
        % Calcul du TEB
        % ----------------------------------------------------------------------------
        serfcp=0.5*erfc((Cp/Lissagep)/sqrt(2));
        serfcn=0.5*erfc((-1*(Cn/Lissagen))/sqrt(2));
        TEB(ii)= (pip/Np)*sum(serfcp)+(pin/Nn)*sum(serfcn)
        
end

%Courbe de la probabilité d'erreur "exacte" qui permettra une
%comparaison avec la simulation pas la méthode du noyau afin d'évaluer son
%efficacite.
%pe_mdp2 = 0.5*erfc(sqrt(10.^([0:11]/10)));		


figure(1);
%semilogy(0:11, pe_mdp2, 'r-', Eb_N0_dB, TEB, 'bo-');
semilogy( Eb_N0_dB, TEB);

xlabel('Eb/No (dB)'); ylabel('TEB');
grid on; axis([0 10 1e-6 1]);
%legend('PEB', 'Simulation Méthode du noyau');
