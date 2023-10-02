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
s={'-o','-h','-<','-d','->','-p','-s'};
sp={'-.o','-.h','-.<','-.d','-.>','-.p','-.s'};
Vecsb=["s" "o" "h" "<" "d" ">" "p"]; % décalage de 1 à cause du modulo 7
%sb=Vecsb(mod(j,7)+1);% tikz p -> x

%% Figures souhaitées
Fig=[];
% Tracé de g pour différents kappa
Fig=[Fig 1];

% Tracé des profils en n et r pour différents kappa
%Fig=[Fig 2];

%Tracé des orbites
%Fig=[Fig 3];

% Tracé de r(n)-tau n pour 1 kappa fixé
%Fig=[Fig 4];

%Tracé de tau_# en fonction de kappa
%Fig=[Fig 5];

figg=[];
figr=[];
figm=[];
figo=[];
figtau=[];

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

%k = r_ast*u_ast^2/p(n_ast);

%% Fonctions r(n), r'(n), g(n)

r = @(n,kappa) (kappa./(1+kappa-p(n)));
rp = @(n,kappa) (kappa*pp(n)./(1+kappa-p(n)).^2);
g = @(n,kappa) ((1+kappa-p(n))./kappa-1./n);

%% Calcul de nbar

n_bar= @(kappa,n0) (sol(@(n) (p(n)-1-kappa),n0));

%% Calcul du 2e zéro de g -> n_x=n_c (n_cross)

n_c = @(kappa,n0) (sol(@(n) g(n,kappa),n0));

%% Calcul de h, n_#=n_ht, tau_#

h = @(n) (p(n)+n.*pp(n));

n_ht = @(kappa,n0) (sol(@(n) (h(n)-1-kappa),n0));
tau_ht = @(kappa,n0) (r(n_ht(kappa,n0))/n_ht(kappa,n0));

nht=[];
tauht=[];
nc=[];

k=[0.5;0.75;1;2;3;4;6];

m=1;

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
    if (or(sum(Fig==2)==1,or(sum(Fig==3)==1,sum(Fig==4)==1)))
        
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
    m=m+1;
end % fin boucle en kappa

%% Figure 4 : g, n_ht, tau_ht, positivité de r-tau n
m=1;
if (sum(Fig==4)==1)
    figure(4)
    
    kappa=2;
    nbar=n_bar(kappa,1)
    n=linspace(0,0.8*nbar,100);
    
    nht=n_ht(kappa,1);
    
    for t=[0.7*nht 0.9*nht nht 1.2*nht]
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
