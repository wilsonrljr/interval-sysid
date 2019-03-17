%{
Article: Identification of Dynamic Systems with Interval Arithmetic

download: https://www.researchgate.net/publication/319056862_Identification_of_Dynamic_Systems_with_Interval_Arithmetic
Uso: perform the identification of a duffing circuit using interval arithmetic

see also: 

1) More in depth studies can be found in the Marcia Peixoto's thesis:
https://ufsj.edu.br/portal2-repositorio/File/ppgel/188-2018-12-17-DissertacaoMarciaPeixoto.pdf

2) Routines of the above thesis
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
load duffing.dat;
y=duffing(1:1:end);
y=decimate(y,5);
yid = y(1:end/2);
yv = y(end/2+1:end);

for k = 7:length(yid)
    
    psin(k-6,:)= [yid(k-1), yid(k-2), yid(k-3), yid(k-4), yid(k-5), yid(k-6), ...
        yid(k-6)*yid(k-1)^2, yid(k-3)^3, yid(k-1)^3,yid(k-5)^3, yid(k-6)^3, ...
        yid(k-4)^3, yid(k-2)^3, yid(k-1)*yid(k-2)^2, yid(k-5)*yid(k-1)^2,...
        yid(k-3)*yid(k-2)*yid(k-1), yid(k-4)*yid(k-2)*yid(k-1),...
        yid(k-6)*yid(k-2)*yid(k-1)];
    
    psi(k-6,:) = [intval(yid(k-1)), intval(yid(k-2)), intval(yid(k-3)), intval(yid(k-4)), ...
        intval(yid(k-5)), intval(yid(k-6)), intval(yid(k-6)*yid(k-1)^2), ...
        intval(yid(k-3)^3),intval(yid(k-1)^3),intval(yid(k-5)^3),...
        intval(yid(k-6)^3),intval(yid(k-4)^3), intval(yid(k-2)^3),...
        intval(yid(k-1)*yid(k-2)^2),intval(yid(k-5)*yid(k-1)^2),...
        intval(yid(k-3)*yid(k-2)*yid(k-1)), intval(yid(k-4)*yid(k-2)*yid(k-1)),...
        intval(yid(k-6)*yid(k-2)*yid(k-1))];
end

t = 7;
yMQ = yid(t:length(yid));

for i=1
    teta(:,i) =inv(psi'*psi)*psi'*yMQ';     % Regressores Intervalares
    tetan(:,i)=inv(psin'*psin)*psin'*yMQ';  % Regressores Nominais
    %xi = yMQ' - psi*teta(:,1);             % Resíduos
end

%% Validação

for k = 7:length(yv)
    
    psiv(k-6,:) = [yv(k-1), yv(k-2), yv(k-3), yv(k-4), yv(k-5), yv(k-6), yv(k-6)*yv(k-1)^2,...
        yv(k-3)^3, yv(k-1)^3, yv(k-5)^3, yv(k-6)^3, yv(k-4)^3, yv(k-2)^3, ...
        yv(k-1)*yv(k-2)^2, yv(k-5)*yv(k-1)^2, yv(k-3)*yv(k-2)*yv(k-1),...
        yv(k-4)*yv(k-2)*yv(k-1), yv(k-6)*yv(k-2)*yv(k-1)];
end
k = 7:length(yv);
yvMQ = yv(7:length(yv));
m=psiv*teta;
f=psiv*tetan;
figure()
plot(k(500:700),m(500:700)','b')
hold on
plot(k(500:700),yvMQ(500:700),'k','LineWidth',1)
hold on
plot(k(500:700),f(500:700),'--r','LineWidth',1)
axis([510,700 , -4, 4]);
box off
set(gca,'FontSize',14,'Fontname','Times New Roman')
xlabel('$k$', 'Interpreter','LaTex','Fontsize',16)
ylabel('y(k)', 'Interpreter','LaTex','Fontsize',16)
axes('position',[.65 .175 .25 .25])
box on
indexOfInterest = (k>560) & (k<570);
plot(k(indexOfInterest),m(indexOfInterest)')
hold on
plot(k(indexOfInterest),yvMQ(indexOfInterest),'k')
hold on
plot(k(indexOfInterest),f(indexOfInterest),'--r')
axis tight
xl = get(gca,'YTickLabel');
new_xl = strrep(xl(:),'.',',');
set(gca,'YTickLabel',new_xl)

%% RMSE

r=psiv*teta; % dados
r1=r(1:length(yvMQ))';
media=mean(yvMQ);
RMSEi= sqrt(sum((yvMQ-r1).^2))/sqrt((sum((yvMQ-media).^2)));

m=psiv*tetan; % dados
m1=m(1:length(yvMQ))';
media=mean(yvMQ);
RMSEn= sqrt(sum((yvMQ-m1).^2))/sqrt((sum((yvMQ-media).^2)))

%% Validação por simulação livre nominal


g(1)=infsup(yv(1),yv(1));
g(2)=infsup(yv(2),yv(2));
g(3)=infsup(yv(3),yv(3));
g(4)=infsup(yv(4),yv(4));
g(5)=infsup(yv(5),yv(5));
g(6)=infsup(yv(6),yv(6));

h(1)=infsup(yv(1),yv(1));
h(2)=infsup(yv(2),yv(2));
h(3)=infsup(yv(3),yv(3));
h(4)=infsup(yv(4),yv(4));
h(5)=infsup(yv(5),yv(5));
h(6)=infsup(yv(6),yv(6));

for k=7:600
    %Intervalar
    h(k) = teta(1)*h(k-1)+teta(2)*h(k-2) + teta(3)*h(k-3)+ teta(4)*h(k-4)+ teta(5)*h(k-5)+...
        teta(6)*h(k-6)+teta(7)*h(k-6)*h(k-1)*h(k-1)+teta(8)*h(k-3)^3+teta(9)*h(k-1)^3+...
        teta(10)*h(k-5)^3+teta(11)*h(k-6)^3+ teta(12)*h(k-4)^3+...
        teta(13)*h(k-2)*h(k-2)*h(k-2)+teta(14)*h(k-2)*h(k-2)*(h(k-1))+ ...
        teta(15)*h(k-5)*h(k-1)^2 +teta(16)*h(k-3)*h(k-2)*h(k-1)+ ...
        teta(17)*h(k-4)*h(k-2)*h(k-1) +teta(18)*h(k-6)*h(k-2)*h(k-1);
end