%{
Article: Identification of Dynamic Systems with Interval Arithmetic

download: https://www.researchgate.net/publication/319056862_Identification_of_Dynamic_Systems_with_Interval_Arithmetic
Uso: perform the identification of a dc-dc motor/generator using interval arithmetic

see also: 
1) More details about how the model were built in the article below.
https://www.researchgate.net/publication/320418710_Identificacao_de_um_motorgerador_CC_por_meio_de_modelos_polinomiais_autorregressivos_e_redes_neurais_artificiais?_sg=tHVRtTZlh0Vyg2RY6xe4kP8xIEIQOJ1DI-bSi0H7pi4A4Jgj5vFGjb89-z_sXZOz5gPaFmXXjJZLZUxA4YEbd_zwRfBPnodQAiBMnxm4.5xRbYXJVv0n5i2E6eDnvs33Zyg0fDJLZOOzjcDWEwXTWjfilPT1SLni15oyl3zdQuhRzsFspkV3TyJFgfGTi4w

2) More in depth studies can be found in the Marcia Peixoto's thesis:
https://ufsj.edu.br/portal2-repositorio/File/ppgel/188-2018-12-17-DissertacaoMarciaPeixoto.pdf

3) Routines of the above thesis
https://ufsj.edu.br/gcom/peixoto2018.php

Other m-files required: intlab toolbox
Subfunctions: none

Authors: Márcia Lucian da Costa Peixoto, Marco Túlio R. Matos, Wilson Rocha {Lacerda Junior} and Samir Angelo Milani Martins
and Erivelton Geraldo Nepomuceno
Website: http://www.ufsj.edu.br/gcom

Please send suggestions for improvement of the above code
to Wilson Rocha at this email address: wilsonrljr@outlook.com

References:
-----------

@InProceedings{PMJ+2017,
  author        = {Peixoto, M{\'{a}}rcia L. C. and Matos, Marco T. R. and J{\'{u}}nior, Wilson R. Lacerda and Martins, Samir A. M. and Nepomuceno, Erivelton G.},
  title         = {{Identification of Dynamic Systems with Interval Arithmetic}},
  booktitle     = {XIII Simp{\'{o}}sio Brasileiro de Automa{\c{c}}{\~{a}}o Inteligente},
  year          = {2017},
  pages         = {1--6},
  address       = {Porto Alegre},
  abstract      = {This paper aims to identify three electrical systems: a series RLC circuit, a motor/generator coupled system, and the Duffing-Ueda oscillator. In order to obtain the system's models was used the error reduction ratio and the Akaike information criterion. Our approach to handle the numerical errors was the interval arithmetic by means of the resolution of the least squares estimation. The routines was implemented in Intlab, a Matlab toolbox devoted to arithmetic interval. Finally, the interval RMSE was calculated to verify the quality of the obtained models. The applied methodology was satisfactory, since the obtained intervals encompass the system's data and allow to demonstrate how the numerical errors affect the answers.},
  annote        = {Lacerda Junior, W. R., Martins, S. A. M. and Nepomuceno, E. G. (2016), "Influence of Sample Rate and Discretization Methods in the Parameter Identification of Systems with Hysteresis", Journal of Applied Nonlinear Dynamics. Rodrigues Junior H. M., Peixoto, M. L. C, Nepomuceno, E. G. and Martins, S. A. M. (2016), "Using Different Interval Extensions to Increase the Accuracy of the Exact Solution on Recursive Functions", Discontinuity, Nonlinearity, and Complexity. Silva, M. R., Nepomuceno, E. G., Amaral, G. F. V. and Martins, S. A. M. (2016), " Exploiting the rounding mode of floating-point in the simulation of Chua's circuit'. Discontinuity, Nonlinearity, and Complexity.},
  archiveprefix = {arXiv},
  arxivid       = {1708.03214},
  eprint        = {1708.03214},
  file          = {:D$backslash$:/User/Google Drive/mendeley/pdfs/Peixoto et al. - 2017 - Identification of Dynamic Systems with Interval Arithmetic.pdf},
}

------------- BEGIN CODE --------------
%}
clear; clc; close all;
intvalinit('displayinfsup');
format long
load sanova1ff.dat
load ennova1.dat
y=sanova1ff(1:60:end);
u=ennova1(1:60:end);
y = y./max(y);
u = u./max(u);
uid = u(1:end/2);
yid = y(1:end/2);
uv = u(end/2+1:end);
yv = y(end/2+1:end);

for k = 3:length(yid)
    
    % Modelo NARX nominal
    psin(k-2,:) = [yid(k-1) yid(k-2) uid(k-1) uid(k-1)*yid(k-1) uid(k-2) uid(k-1)*yid(k-2)...
        uid(k-2)*yid(k-1) uid(k-2)*yid(k-2) yid(k-2)^2];
    
    % Modelo NARX intervalar
    psi(k-2,:) = [intval(yid(k-1)) intval(yid(k-2)) intval(uid(k-1))
        intval(uid(k-1)*yid(k-1)) intval(uid(k-2)) intval(uid(k-1)*yid(k-2))...
        intval(uid(k-2)*yid(k-1)) intval(uid(k-2)*yid(k-2)) intval(yid(k-2)^2)];
end

yMQ = yid(3:length(yid));

for i=1
    teta(:,i) =inv(psi'*psi)*psi'*yMQ;     % Regressores Intervalares
    xi = yMQ - psi*teta(:,1);              % Resíduos
    tetan(:,i)=inv(psin'*psin)*psin'*yMQ;  % Regressores tradicionais
end

%% % Validação Predição livre

g(1)=yv(1);
g(2)=yv(2);      % Nominal

h(1)=intval(yv(1)); % Intervalar
m(1)=h(1);
h(2)=intval(yv(2));
m(2)=h(2);
uval=uv;

for k=3:1:length(yv)
    
    m(k) = teta(1)*m(k-1)+ teta(2)*m(k-2)+teta(3)*uval(k-1)+teta(4)*uval(k-1)*m(k-1)+ ...
        teta(5)*uval(k-2)+ teta(6)*uval(k-1)*m(k-2)+ teta(7)*uval(k-2)*m(k-1)+ ...
        teta(8)*uval(k-2)*m(k-2)+teta(9)*m(k-2)^2; % Intervalar
    
    g(k) = tetan(1)*g(k-1)+ tetan(2)*g(k-2)+tetan(3)*uval(k-1)+tetan(4)*uval(k-1)*g(k-1)+...
        tetan(5)*uval(k-2)+ tetan(6)*uval(k-1)*g(k-2)+ tetan(7)*uval(k-2)*g(k-1)+...
        tetan(8)*uval(k-2)*g(k-2)+tetan(9)*g(k-2)^2; % Nominal
    
    xd(k) = gradientinit(h(k-1));
    
    td(k) = teta(1)*xd(k-1)+ teta(2)*xd(k-2)+ teta(3)*uval(k-1)+teta(4)*uval(k-1)*xd(k-1)+...
        teta(5)*uval(k-2) +teta(6)*uval(k-1)*xd(k-2)+ teta(7)*uval(k-2)*xd(k-1)+ ...
        teta(8)*uval(k-2)*xd(k-2)+teta(9)*xd(k-2)^2;
    
    h(k)  = teta(1)*mid(h(k-1))+ teta(2)*mid(h(k-2))+teta(3)*uval(k-1)+...
        teta(4)*uval(k-1)*mid(h(k-1))+teta(5)*uval(k-2)+ teta(6)*uval(k-1)*mid(h(k-2))+...
        teta(7)*uval(k-2)*mid(h(k-1))+teta(8)*uval(k-2)*mid(h(k-2))+...
        teta(9)*mid(h(k-2))^2+td(k).dx*(h(k-1)-mid(h(k-1))); % MVF
    
end

figure
plot(h)
hold on
plot(yv,'k','LineWidth',1)
hold on
plot(g,'--r','LineWidth',1)
axis([2000 3000 0 1.2])
box off
set(gca,'FontSize',14,'Fontname','Times New Roman')
xlabel('$k$', 'Interpreter','LaTex','Fontsize',16)
ylabel('Velocidade angular (rpm)', 'Interpreter','LaTex','Fontsize',16)
xl = get(gca,'YTickLabel');
new_xl = strrep(xl(:),'.',',');
set(gca,'YTickLabel',new_xl)
print(gcf,'-depsc','fig5')

%%
% RMSE Tradicional

g1=g'; % dados da simulação livre
media=mean(yv); % dados coletados
RMSE=sqrt(sum((yv-g1).^2))/sqrt((sum((yv-media).^2)))

% RMSE Intervalar Simulação Livre
r=h';
media=mean(yv);
RMSEi= sqrt(sum((yv-r).^2))/sqrt((sum((yv-media).^2)));
Rmse=[RMSEi.inf RMSEi.sup]
