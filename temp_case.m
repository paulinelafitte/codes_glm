%% Réinitialisation
clear all
close all
options = optimset('Display','off');
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
sol=@fzero; % ou fzero ou sol

set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultLineLineWidth',1.5)
set(0,'DefaultAxesColorOrder',[0 0 0])

stau={'. ','-.','-','--'};
seps={'- ','-.','.','--'};
s={'-o','-h','-<','-d','->','-p','-s'};
sp={'-.o','-.h','-.<','-.d','-.>','-.p','-.s'};
Vecsb=["s" "o" "h" "<" "d" ">" "p"]; % décalage de 1 à cause du modulo 7
%sb=Vecsb(mod(j,7)+1);% tikz p -> x

%% Figures souhaitées
Fig=[];
% Tracé de g pour différents kappa
%Fig=[Fig 1];

% Tracé des profils en n et r pour différents kappa
%Fig=[Fig 2];

%Tracé des orbites
%Fig=[Fig 3];

%Tracé de r(n)-tau n pour 1 kappa fixé
%Fig=[Fig 4];

%Tracé de tau_# en fonction de kappa
%Fig=[Fig 5];

%Tracé des profils pour theta>0
%Fig=[Fig 6];

%Tracé des orbites pour theta>0
%Fig=[Fig 7];

%Tracé de Delta et rplus, rminus
Fig=[Fig 8];

%Tracé des profils et des orbites
Fig=[Fig 9 10];

figg=[];
figr=[];
figm=[];
figo=[];
figorb=[];
figtau=[];
figeps=[];
nracine=[];
EQ=[];

for fig=Fig
    figure(fig)
    clf
end

%% Pression p, p'

gamma = 2;

p = @(n) (n.^gamma);
pp = @(n) (gamma*(n.^(gamma-1)));

%% Choix de n_*, calcul de kappa_*, fonctions r_*(tau), u_*(tau,kappa)

n_ast = 1;
k_ast = n_ast*pp(n_ast)/p(n_ast);

r_ast = @(tau) n_ast/tau;
u_ast= @(tau,kappa) sqrt(tau*kappa*p(n_ast)/n_ast);

eps = @(kappa,theta,uast) (kappa*theta/uast.^2);
%k = r_ast*u_ast^2/p(n_ast);

%% Fonctions r(n), r'(n), g(n)

r = @(n,kappa) (kappa./(1+kappa-p(n)));
rp = @(n,kappa) (kappa*pp(n)./(1+kappa-p(n)).^2);
g = @(n,kappa) ((1+kappa-p(n))./kappa-1./n);
gp = @(n,kappa) (-pp(n)/kappa+1./n.^2);

%% Fonctions T, B0, B1, B, C0, C1, C, Delta, rp, rm

%T = @(n,r) (1./r-1./n);
T = @(nr) (1./nr(2)-1./nr(1));
% B0 = @(n,r,kappa) ((1+kappa-p(n))./kappa-1./r);
% B1 = @(n,r,kappa,tau) ((1+tau*(n-1)-r)/kappa);
% B = @(n,r,kappa,tau,epsilon) (B0(n,r,kappa)+epsilon*B1(n,r,kappa,tau));
B0 = @(nr,kappa) ((1+kappa-p(nr(1)))./kappa-1./nr(2));
B1 = @(nr,kappa,tau) ((1+tau*(nr(1)-1)-nr(2))/kappa);
B = @(nr,kappa,tau,epsilon) (B0(nr,kappa)+epsilon*B1(nr,kappa,tau));
C0 = @(n,kappa) (1+kappa-p(n));
C1 = @(n,tau) (1+tau*(n-1));
C = @(n,kappa,tau,epsilon) (C0(n,kappa)+epsilon*C1(n,tau));
Delta = @(n,kappa,tau,epsilon) (C(n,kappa,tau,epsilon).^2-4*kappa*epsilon);
r_plus = @(n,kappa,tau,epsilon) ((C(n,kappa,tau,epsilon)+sqrt(Delta(n,kappa,tau,epsilon)))/2/epsilon);
r_minus = @(n,kappa,tau,epsilon) (2*kappa./(C(n,kappa,tau,epsilon)+sqrt(Delta(n,kappa,tau,epsilon))));

%% Champ de vecteurs

F= @(nr,kappa,tau,epsilon) ([T(nr)+tau*nr(1).*B(nr,kappa,tau,epsilon)./(nr(2)-tau*nr(1)); ...
    (pp(nr(1))+epsilon*tau^2*nr(1)./(nr(2)-tau*nr(1))).*B(nr,kappa,tau,epsilon)/epsilon]);

%% Calcul de nbar

n_bar= @(kappa,n0) (sol(@(n) (p(n)-1-kappa),n0));

%% Calcul du 2e zéro de g -> n_x=n_c (n_cross)

n_c = @(kappa,n0) (sol(@(n) g(n,kappa),n0));

n_c1 = @(kappa,tau,n0) ((1-tau)*(n_c(kappa,n0)-1)/kappa/gp(n_c(kappa,n0),kappa));

%% Calcul de h, n_#=n_ht, tau_#

h = @(n) (p(n)+n.*pp(n));

n_ht = @(kappa,n0) (sol(@(n) (h(n)-1-kappa),n0));
tau_ht = @(kappa,n0) (r(n_ht(kappa,n0),kappa)/n_ht(kappa,n0));

%% Calcul de ntilde qui annule Delta

n_tilde = @(kappa,epsilon,tau,n0) (sol( @(n) (C(n,kappa,tau,epsilon).^2 -4*kappa*epsilon),n0))

nht=[];
tauht=[];
nc=[];

k=[0.5;0.75;1;2;3;4;6];

m=1;

if (sum(Fig==1)==1|sum(Fig==2)==1|sum(Fig==3)==1|sum(Fig==6)==1|sum(Fig==7)==1)
    for kappa=k'

        sb=Vecsb(mod(m,7)+1);

        nbar=n_bar(kappa,1);
        nht = [nht n_ht(kappa,2)];
        tauht = [tauht kappa/(nht(end)^2*pp(nht(end)))];

        if (kappa>k_ast) %
            nc = [nc n_c(kappa,0.7*nbar)];
            ncross = 1.001;
        else
            if (kappa<k_ast)
                nc = [nc n_c(kappa,0.5)];
                ncross = 1.001*nc(end);
            else
                nc = [nc 1];
                ncross = 1;
            end
        end

        %% Figure 1 : plot g

        if (sum(Fig==1)==1)
            figure(1)
            nn=linspace(0.3,nbar,100);
            plot(nn,g(nn,kappa),"-"+sb,'MarkerIndices',[1:5:length(nn) length(nn)])
            hold on
            figg(m)=scatter(nc(end),0,100,'k',sb,'filled','DisplayName',"$\kappa="+num2str(kappa)+"$");
            hold on
        end


        %% Figures 2 et 3 : Tracé des orbites hétéroclines
        if ((sum(Fig==2)==1|sum(Fig==3)==1)|(sum(Fig==6)==1|sum(Fig==7)==1))

            solution=ode23s(@(t,n) g(n,kappa),[0 30],ncross);
            temps = solution.x;
            profil_n = solution.y;



            if (sum(Fig==2)==1)
                figure(2)
                plot(temps,profil_n,'k',temps(1:3:end),profil_n(1:3:end),strcat('k',sb),'linewidth', 1.2);
                hold on
                plot(temps,r(profil_n,kappa),'k:',temps(1:3:end),r(profil_n(1:3:end),kappa),strcat('k',sb),'linewidth', 1.2);
                hold on
                figr(m)=scatter([temps(1); temps(end)],[profil_n(1); profil_n(end)],75,'k',sb,'filled','DisplayName',"$\kappa="+num2str(kappa)+"$");
                hold on
            end
            if (sum(Fig==3)==1)
                figure(3)
                if (kappa<2)
                    subplot(1,2,1)
                    plot(profil_n,r(profil_n,kappa),cell2mat(s(m)),'MarkerIndices',[1:5:length(profil_n) length(profil_n)])
                    hold on
                    scatter(nht(m),r(nht(m),kappa),100,'k',sb,'LineWidth',1.5)
                    hold on
                    figm(m)=scatter([profil_n(1); profil_n(end)],[r(profil_n(1),kappa); r(profil_n(end),kappa)],75,'k',sb,'filled','DisplayName',"$\kappa="+num2str(kappa)+"$");
                    hold on
                    plot(profil_n,tauht(m)*profil_n,"k"+cell2mat(sp(m)),'linewidth',1,'MarkerIndices',[1:5:length(profil_n) length(profil_n)])
                    hold on
                end
                if (kappa>2)
                    subplot(1,2,2)
                    plot(profil_n,r(profil_n,kappa),cell2mat(s(m)),'MarkerIndices',[1:5:length(profil_n) length(profil_n)])
                    hold on
                    scatter(nht(m),r(nht(m),kappa),100,'k',sb,'LineWidth',1.5)
                    hold on
                    figo(m)=scatter([profil_n(1); profil_n(end)],[r(profil_n(1),kappa); r(profil_n(end),kappa)],75,'k',sb,'filled','DisplayName',"$\kappa="+num2str(kappa)+"$");
                    hold on
                    plot(profil_n,tauht(m)*profil_n,"k"+cell2mat(sp(m)),'linewidth',1,'MarkerIndices',[1:5:length(profil_n) length(profil_n)])
                    hold on
                end
            end
        end
        %% Figures 6&7 : avec température, n,r en fonction de y / r en fonction de n
        if (or(sum(Fig==6)==1,sum(Fig==7)==1))
            epsilon=0.00001;
            tau=tauht(end)/(2+epsilon);


            solutionnr=ode15s(@(t,nr) F(nr,kappa,tau,epsilon),[0 10],[ncross;ncross]+epsilon^2);

            tempsnr = solutionnr.x;

            profilnr_n = solutionnr.y(1,:);
            profilnr_r = solutionnr.y(2,:);
            if (sum(Fig==6)==1)
                figure(6)
                plot(tempsnr,profilnr_n,'k',tempsnr(1:3:end),profilnr_n(1:3:end),strcat('k',sb),'linewidth', 1.2);
                hold on
                plot(tempsnr,profilnr_r,'k:',tempsnr(1:3:end),profilnr_r(1:3:end),strcat('k',sb),'linewidth', 1.2);
                hold on
                figr(m)=scatter([tempsnr(1); tempsnr(end)],[profilnr_n(1); profilnr_n(end)],75,'k',sb,'filled','DisplayName',"$\kappa="+num2str(kappa)+"$");
                scatter([tempsnr(1); tempsnr(end)],[profilnr_r(1); profilnr_r(end)],75,'k',sb,'filled','DisplayName',"$\kappa="+num2str(kappa)+"$");
                hold on

                %             plot(temps,profil_n,'r',temps(1:3:end),profil_n(1:3:end),strcat('r',sb),'linewidth', 1.2);
                %             hold on
                %             plot(temps,r(profil_n,kappa),'r:',temps(1:3:end),r(profil_n(1:3:end),kappa),strcat('r',sb),'linewidth', 1.2);
                %             hold on
                %             scatter([temps(1); temps(end)],[profil_n(1); profil_n(end)],75,'r',sb,'filled','DisplayName',"$\kappa="+num2str(kappa)+"$");
                %             hold on
                %
            end
            if (sum(Fig==7)==1)
                figure(7)
                if (kappa<2)
                    subplot(1,2,1)
                    plot(profilnr_n,r(profilnr_n,kappa),cell2mat(s(m)),'MarkerIndices',[1:5:length(profilnr_n) length(profilnr_n)])
                    hold on
                    scatter(nht(m),r(nht(m),kappa),100,'k',sb,'LineWidth',1.5)
                    hold on
                    figm(m)=scatter([profilnr_n(1); profilnr_n(end)],[r(profilnr_n(1),kappa); r(profilnr_n(end),kappa)],75,'k',sb,'filled','DisplayName',"$\kappa="+num2str(kappa)+"$");
                    hold on
                    plot(profilnr_n,tauht(m)*profilnr_n,"k"+cell2mat(sp(m)),'linewidth',1,'MarkerIndices',[1:5:length(profilnr_n) length(profilnr_n)])
                    hold on

                    %                 plot(profil_n,r(profil_n,kappa),"r"+cell2mat(s(m)),'MarkerIndices',[1:5:length(profil_n) length(profil_n)])
                    %                 hold on
                    %                 scatter(nht(m),r(nht(m),kappa),100,'r',sb,'LineWidth',1.5)
                    %                 hold on
                    %                 scatter([profil_n(1); profil_n(end)],[r(profil_n(1),kappa); r(profil_n(end),kappa)],75,'r',sb,'filled','DisplayName',"$\kappa="+num2str(kappa)+"$");
                    %                 hold on
                    %                 plot(profil_n,tauht(m)*profil_n,"r"+cell2mat(sp(m)),'linewidth',1,'MarkerIndices',[1:5:length(profil_n) length(profil_n)])
                    %                 hold on

                end
                if (kappa>2)
                    subplot(1,2,2)
                    plot(profilnr_n,r(profilnr_n,kappa),cell2mat(s(m)),'MarkerIndices',[1:5:length(profilnr_n) length(profilnr_n)])
                    hold on
                    scatter(nht(m),r(nht(m),kappa),100,'k',sb,'LineWidth',1.5)
                    hold on
                    figo(m)=scatter([profilnr_n(1); profilnr_n(end)],[r(profilnr_n(1),kappa); r(profilnr_n(end),kappa)],75,'k',sb,'filled','DisplayName',"$\kappa="+num2str(kappa)+"$");
                    hold on
                    plot(profilnr_n,tauht(m)*profilnr_n,"k"+cell2mat(sp(m)),'linewidth',1,'MarkerIndices',[1:5:length(profilnr_n) length(profilnr_n)])
                    hold on

                    %                 plot(profil_n,r(profil_n,kappa),"r"+cell2mat(s(m)),'MarkerIndices',[1:5:length(profil_n) length(profil_n)])
                    %                 hold on
                    %                 scatter(nht(m),r(nht(m),kappa),100,'r',sb,'LineWidth',1.5)
                    %                 hold on
                    %                 scatter([profil_n(1); profil_n(end)],[r(profil_n(1),kappa); r(profil_n(end),kappa)],75,'r',sb,'filled','DisplayName',"$\kappa="+num2str(kappa)+"$");
                    %                 hold on
                    %                 plot(profil_n,tauht(m)*profil_n,"r"+cell2mat(sp(m)),'linewidth',1,'MarkerIndices',[1:5:length(profil_n) length(profil_n)])
                    %                 hold on
                end
            end
        end
        m=m+1;

    end % fin boucle en kappa
end %fin condition boucle

%% Figure 4 : g, n_ht, tau_ht, positivité de r-tau n
m=1;
if (sum(Fig==4)==1)
    figure(4)

    kappa=2;
    nbar=n_bar(kappa,1)
    n=linspace(0,0.8*nbar,100);

    nht=n_ht(kappa,1);
    tauht = tau_ht(kappa,1.1)

    for t=[0.7*tauht 0.9*tauht tauht 1.2*tauht]
        figtau(m)=plot(n,r(n,kappa)-t*n,cell2mat(stau(m)),'MarkerIndices',[1:5:length(n) length(n)],'linewidth', 1.5,'DisplayName',"$\tau="+num2str(t)+"$")
        hold on
        m=m+1;
    end
    xlim([-0.1 1.5])
    ylim([-0.7 1.6])
    n=linspace(0.3,0.8*nbar,100);
    figtau(m)=plot(n,g(n,kappa),'k-+','MarkerIndices',[1:4:length(n) length(n)],'DisplayName',"g")
    hold on
    m=m+1;
    figtau(m)=plot(n,r(n,kappa)./n,'k-v','MarkerIndices',[1:4:length(n) length(n)],'DisplayName',"n$\,\mapsto\,$r(n)/n")
    m=m+1;
    hold on
    figtau(m)=plot(n(n>0.7),rp(n(n>0.7),kappa),'k-^','MarkerIndices',[1:4:length(n(n>0.7)) length(n(n>0.7))],'linewidth', 1.5,'DisplayName',"n$\,\mapsto\,$r'(n)")
    hold on

end

%% Figure 5 : tau_ht en fonction de kappa

if (sum(Fig==5)==1)
    figure(5)
    ka=linspace(0,12,100);
    % for gamma=[1.4 2 3]
    tau_gam=@(ka) ((ka/gamma).*((1+gamma)./(1+ka)).^(1+1/gamma))
    plot(ka,tau_gam(ka))
    hold on
    %end
    plot([2 2],[0 1],'k:')
    hold on
    area(ka(ka<=2),tau_gam(ka(ka<=2)),'FaceColor',[0.96 0.96 0.96])
    area(ka(ka>=2),tau_gam(ka(ka>=2)),'FaceColor',[0.88 0.88 0.88])
end



%% Figure 8 : avec température, Delta, rplus, rminus en fonction de n a epsilon variable
if (sum(Fig==8)==1|sum(Fig==9)==1|sum(Fig==10)==1)

    kappa=3
    tauht = tau_ht(kappa,1.1)
    tauxtauht=0.9
    tau=tauxtauht*tauht;%0.5%tau_ht(kappa,1.1) ;
    eps=[0.05 0.1 0.1582];% 0.00002 0.00003];

    mineq=1;
    maxeq=1;

    nbar=n_bar(kappa,1);
    nb0=linspace(0,0.9*nbar,100);
    
    m=1;
    %for tau=[0.1*tauht 0.5*tauht tauht 1.05*tauht]

    if (kappa>k_ast) %
        EQ=n_c(kappa,0.7*nbar);
        init=1;
    else
        if (kappa<k_ast)
            EQ = n_c(kappa,0.5);
            init=EQ;
        else
            EQ = 1; % uninteresting case (kappa=2)
        end
    end

    for epsilon=eps
        equil=sol(@(n) r_minus(n,kappa,tau,epsilon)-n,EQ(m))
        EQ=[EQ equil]
    end

    limites=[min(1,0.995*min(EQ)) max(1,1.005*max(EQ))];
    nb1=linspace(limites(1),limites(2),100);    


    if (sum(Fig==8)==1)
        figure(8)
        subplot(2,2,1)
        plot(nb0,Delta(nb0,kappa,tau,0),'r-.',nb0(1:20:end),Delta(nb0(1:20:end),kappa,tau,0),Vecsb(mod(m,7)+1),'linewidth', 1.5);
        hold on
        scatter([nb0(1) nb0(end)],[Delta(nb0(1),kappa,tau,0) Delta(nb0(end),kappa,tau,0)],75,'r',Vecsb(mod(m,7)+1),'filled');
        hold on
        subplot(2,2,2)
        plot(nb0,C(nb0,kappa,tau,0),'r-.',nb0(1:20:end),C(nb0(1:20:end),kappa,tau,0),Vecsb(mod(m,7)+1),'linewidth', 1.5);
        hold on
        scatter([nb0(1) nb0(end)],[C(nb0(1),kappa,tau,0) C(nb0(end),kappa,tau,0)],75,'r',Vecsb(mod(m,7)+1),'filled');
        hold on
        subplot(2,2,3)
        plot(nb0,nb0,'k-','linewidth',0.5)
        hold on
        plot(nb0,r(nb0,kappa),'r-.','linewidth', 1.5);
        hold on
        scatter([nb0(1); nb0(end)],[r(nb0(1),kappa) r(nb0(end),kappa)],75,'r',Vecsb(mod(m,7)+1),'filled');
        hold on
        subplot(2,2,4)
        plot(nb1,r(nb1,kappa),'k',nb1(1:50:end),r(nb1(1:50:end),kappa),Vecsb(mod(m,7)+1),'linewidth', 0.7);
        hold on
        figeps(m)=scatter(EQ(1),r(EQ(1),kappa),75,'k',Vecsb(mod(m,7)+1),'filled','DisplayName',"$\varepsilon=0$");
        hold on
    end

    solution=ode23s(@(t,n) g(n,kappa),[0 55],init+1e-4);
    temps = solution.x;
    profil_n = solution.y;

    if (sum(Fig==9)==1)
        figure(9)

        plot(temps,profil_n,'r',temps(1:3:end),profil_n(1:3:end),'r'+Vecsb(mod(m,7)+1),'linewidth', 1.5);
        hold on
        plot(temps,r(profil_n,kappa),'r:',temps(1:3:end),r(profil_n(1:3:end),kappa),'rx','linewidth', 1.5);
        hold on
        figr(m)=scatter([temps(1); temps(end)],[profil_n(1); profil_n(end)],75,'r',Vecsb(mod(m,7)+1),'filled','DisplayName',"$\varepsilon=0$");
        hold on
    end

    if (sum(Fig==10)==1)
        figure(10)
        plot(profil_n,r(profil_n,kappa),'r',profil_n(1:3:end),r(profil_n(1:3:end),kappa),Vecsb(mod(m,7)+1),'MarkerIndices',[1:5:length(profil_n) length(profil_n)])
        hold on
        figorb(m)=scatter([profil_n(1); profil_n(end)],[r(profil_n(1),kappa); r(profil_n(end),kappa)],75,'r',Vecsb(mod(m,7)+1),'filled','DisplayName',"$\varepsilon=0$");
        hold on
        plot(profil_n,tau*profil_n,'g')
        hold on
        plot(profil_n,tauht*profil_n,'y')
        hold on
    end

    m=m+1;

    mineq=min(mineq, min(solution.y));
    maxeq=max(maxeq, max(solution.y));

    for epsilon=eps
        ntilde=n_tilde(kappa,epsilon,tau,nbar);
        n=linspace(0.01*nbar,ntilde,100);
        nb=linspace(0,1.2*nbar,100);

        if (sum(Fig==8)==1)
            figure(8)
            subplot(2,2,1)
            plot(nb,Delta(nb,kappa,tau,epsilon),'k',nb(7*m:30:end),Delta(nb(7*m:30:end),kappa,tau,epsilon),Vecsb(mod(m,7)+1),'linewidth', 1);
            hold on
            figeps(m)=scatter([nb(1) nb(end)],[Delta(nb(1),kappa,tau,epsilon) Delta(nb(end),kappa,tau,epsilon)],75,'k',Vecsb(mod(m,7)+1),'filled','DisplayName',"$\varepsilon="+num2str(epsilon)+"$");
            hold on
            subplot(2,2,2)
            plot(nb,C(nb,kappa,tau,epsilon),'k',nb(7*m:30:end),C(nb(7*m:30:end),kappa,tau,epsilon),Vecsb(mod(m,7)+1),'linewidth', 1)
            hold on
            scatter([nb(1) nb(end)],[C(nb(1),kappa,tau,epsilon) C(nb(end),kappa,tau,epsilon)],75,'k',Vecsb(mod(m,7)+1),'filled')
            hold on
            subplot(2,2,3)
            plot(n,r_minus(n,kappa,tau,epsilon),'k',n(1:20:end),r_minus(n(1:20:end),kappa,tau,epsilon),Vecsb(mod(m,7)+1),'linewidth', 1);
            scatter([n(1) n(end)],[r_minus(n(1),kappa,tau,epsilon) r_minus(n(end),kappa,tau,epsilon)],75,'k',Vecsb(mod(m,7)+1),'filled');
            hold on
            plot(n,r_plus(n,kappa,tau,epsilon),'k:',n(1:10:end),r_plus(n(1:10:end),kappa,tau,epsilon),Vecsb(mod(m,7)+1),'linewidth',1.5)
            scatter([n(end)],[r_minus(n(end),kappa,tau,epsilon)],75,'k',Vecsb(mod(m,7)+1),'filled')
            hold on
            subplot(2,2,4)
            plot(nb1,r_minus(nb1,kappa,tau,epsilon),'k',nb1(1:50:end),r_minus(nb1(1:50:end),kappa,tau,epsilon),Vecsb(mod(m,7)+1),'linewidth', 0.7);
            figeps(m)=scatter(EQ(m),r_minus(EQ(m),kappa,tau,epsilon),75,'k',Vecsb(mod(m,7)+1),'filled','DisplayName',"$\varepsilon="+num2str(epsilon)+"$");
            hold on
        end

        if (kappa>k_ast) %
            init = (1+1e-4)*[1;1];
        else
            if (kappa<k_ast)
                init=(EQ(m)+1e-4)*[1;1];
            else
                init = 1; % uninteresting case (kappa=2)
            end
        end

        solutionnr=ode23s(@(t,nr) F(nr,kappa,tau,epsilon),[0 55],init);

        tempsnr = solutionnr.x; profilnr_n = solutionnr.y(1,:); profilnr_r = solutionnr.y(2,:);

        if (sum(Fig==9)==1)
            figure(9)
            plot(tempsnr,profilnr_n,'k',tempsnr(1:3:end),profilnr_n(1:3:end),Vecsb(mod(m,7)+1),'linewidth', 1.2);
            hold on
            plot(tempsnr,profilnr_r,'k:',tempsnr(1:3:end),profilnr_r(1:3:end),Vecsb(mod(m,7)+1),'linewidth', 1.2);
            hold on
            figr(m)=scatter([tempsnr(1); tempsnr(end)],[profilnr_n(1); profilnr_n(end)],75,'k',Vecsb(mod(m,7)+1),'filled','DisplayName',"$\varepsilon="+num2str(epsilon)+"$");
            hold on
        end
        if (sum(Fig==10)==1)
            figure(10)
            plot(profilnr_n,profilnr_r,'k',profilnr_n,profilnr_r,Vecsb(mod(m,7)+1),'MarkerIndices',[1:5:length(profilnr_n) length(profilnr_n)])
            hold on
            figorb(m)=scatter([profilnr_n(1); profilnr_n(end)],[profilnr_r(1); profilnr_r(end)],75,'k',Vecsb(mod(m,7)+1),'filled','DisplayName',"$\varepsilon="+num2str(epsilon)+"$");
            hold on
        end

        mineq=min(mineq, min(min(solutionnr.y)));
        maxeq=max(maxeq, max(max(solutionnr.y)));

        
        (EQ(m)-nc)/epsilon-(1-tau)*(nc-1)/kappa/gp(nc,kappa)
        (1-tau)*(nc-1)/kappa/gp(nc,kappa)
        eeps=1e-6;
        %EQ(m)=1
        C1=F([EQ(m),EQ(m)]+eeps*[1,0],kappa,tau,epsilon)/eeps;
        C2=F([EQ(m),EQ(m)]+eeps*[0,1],kappa,tau,epsilon)/eeps;
        A=[C1 C2]

        if (sum(Fig==8)==1)
            figure(8)
            subplot(2,2,1)

            plot(nb,0*nb,"k-",'linewidth',0.5)
            hold on
            subplot(2,2,4)
            plot(nb1,nb1,'r--','linewidth',1)
            hold on
        end
        m=m+1;  
    end




end


%% Légendes, annotations, ...

if (sum(Fig==1)==1)
    figure(1)
    xlim([-0.1 1.1*max(nc)])
    ylim([-0.5 0.7])
    annotation('arrow',[0 0],[0 1]);
    legend(figg)
    xlabel('n')
    ha = annotation('arrow',[0.13 0.91],[0.45 0.45]);  % store the arrow information in ha
    ha.LineWidth  = 1;          % make the arrow bolder for the picture
    ha.HeadWidth  = 10;
    ha.HeadLength = 10;
    hb = annotation('arrow',[0.175 0.175],[0.11 0.92]);  % store the arrow information in ha
    hb.LineWidth  = 1;          % make the arrow bolder for the picture
    hb.HeadWidth  = 10;
    hb.HeadLength = 10;
    grid on
    title("Graph of g")
end

if (sum(Fig==2)==1)
    figure(2)
    legend(figr,'Location','southeast')
    ylim([0.9*min(nc) 1.01*max(nc)])
    grid on
    xlabel("y")
    ylabel("n")
    title("Profiles of n (plain) and of r(n) (dotted)")
end

if (sum(Fig==3)==1)
    figure(3)
    subplot(1,2,1)
    title("$\kappa<2$")
    xlabel("n")
    ylabel("r")
    legend(figm,'Location','northwest')
    ylim([0.2 1])
    grid on
    subplot(1,2,2)
    title("$\kappa>2$")
    xlabel("n")
    ylabel("r")
    ylim([1 2])
    legend(figo(5:7),'Location','northwest')
    grid on
    sgtitle("Heteroclinic orbits")
end

if (sum(Fig==4)==1)
    xlim([-0.099 1.5])
    ylim([-1.15 1.6])
    annotation('arrow',[0 0],[0 1]);
    legend(figtau,'Location','Southeast')
    xlabel('n')

    ha = annotation('arrow',[0.13 0.91],[0.45 0.45]);  % store the arrow information in ha
    ha.LineWidth  = 1;          % make the arrow bolder for the picture
    ha.HeadWidth  = 10;
    ha.HeadLength = 10;
    hb = annotation('arrow',[0.175 0.175],[0.11 0.92]);  % store the arrow information in ha
    hb.LineWidth  = 1;          % make the arrow bolder for the picture
    hb.HeadWidth  = 10;
    hb.HeadLength = 10;
    grid on
    title("Exploration of the case $\kappa=2$")
end

if (sum(Fig==5)==1)
    xlabel("$\kappa$")
    ylabel("$\tau$")
    text(1.9,-0.023,"$\kappa^*$")
    set(gca,'xtick',[])
    yticks([0 1])
    text(0.5,0.3,"$n_\times<n^*$",'FontSize', 12)
    text(5,0.5,"$n_\times>n^*$",'FontSize', 12)
    ha=annotation('arrow',[0.13 0.13],[0.11 0.92]);
    hb=annotation('arrow',[0.13 0.91],[0.11 0.11]);
    %  legend(figtau,'Location','Southeast')

    ha.LineWidth  = 1;          % make the arrow bolder for the picture
    ha.HeadWidth  = 10;
    ha.HeadLength = 10;
    hb.LineWidth  = 1;          % make the arrow bolder for the picture
    hb.HeadWidth  = 10;
    hb.HeadLength = 10;
    grid on
    title("Graph of $\tau_\#$")
end

if (sum(Fig==6)==1)
    figure(6)
    legend(figr,'Location','southeast')
    ylim([0.9*min(nc) 1.01*max(nc)])
    grid on
    xlabel("y")
    ylabel("n")
    title("Profiles of n (plain) and of r(n) (dotted)")
end

if (sum(Fig==7)==1)
    figure(7)
    subplot(1,2,1)
    title("$\kappa<2$")
    xlabel("n")
    ylabel("r")
    legend(figm,'Location','northwest')
    ylim([0.2 1])
    grid on
    subplot(1,2,2)
    title("$\kappa>2$")
    xlabel("n")
    ylabel("r")
    ylim([1 2])
    legend(figo(5:7),'Location','northwest')
    grid on
    sgtitle("Heteroclinic orbits")
end

if (sum(Fig==8)==1)
    figure(8)
    subplot(2,2,1)
    xlabel("n")
    xlim([min(nb) max(nb)])
    ylim([min(Delta(nb,kappa,tau,0.2)) max(Delta(nb,kappa,tau,0.2))])
    grid on
    title("$\Delta$")
    subplot(2,2,2)
    xlabel("n")
    xlim([min(nb) max(nb)])
    ylim([min(C(nb,kappa,tau,0.2)) max(C(nb,kappa,tau,0.2))])
    grid on
    title("$C$")
    xlim([min(nb) max(nb)])
    subplot(2,2,3)
    xlabel("n")
    ylabel("r")
    ylim([0 17])
    title("$r_-(n)$, $r_+(n)$")
    grid on
    subplot(2,2,4)
    legend(figeps,'Location','northwest')
    xlabel("n")
    ylabel("r")
    
    xlim(limites)
    ylim(limites)
    grid on
    title("$n$, $r_-(n)$")
end

if (sum(Fig==9)==1)
    figure(9)
    xlabel("y")
    ylabel("n,r")
    legend(figr,'Location','northwest')
    xlim([0 55])
    ylim([mineq maxeq])
    grid on
    title("Profiles for $\kappa="+num2str(kappa)+"$ and $\tau="+num2str(tauxtauht)+"*\tau_\#$")
end

if (sum(Fig==10)==1)
    figure(10)
    xlabel("n")
    ylabel("r")
    legend(figorb,'Location','northwest')
    ylim([mineq maxeq])
    xlim([mineq maxeq])
    grid on
    title("Heteroclinic orbits for $\kappa="+num2str(kappa)+"$ and $\tau="+num2str(tauxtauht)+"*\tau_\#$")
end
