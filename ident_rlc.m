%{
Article: Identification of Dynamic Systems with Interval Arithmetic

download: https://www.researchgate.net/publication/319056862_Identification_of_Dynamic_Systems_with_Interval_Arithmetic
Uso: perform the identification of a rlc circuit using interval arithmetic

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
clear; close all; clc
intvalinit('displayinfsup');
load('dados1')

y=alldata(1:end);
u=saida(:,1:end);
y = y./max(y);
u = u./max(u);
ui = u(1:end/2); %ui = input training data
yi = y(1:end/2); % yi = output training data
uv = u(end/2+1:end); %uv = input validation data
yv = y(end/2+1:end); % yv = output validation data
e = (randn(5000,1)); % noise sequence

for k = 3:length(yi)
    psi(k-2,:) = [intval(yi(k-1)), intval(yi(k-2)), intval(ui(k-2))];
end

for k = 3:length(yi) % Matrix of regressors without interval. Matriz de Regressores nominais
    psin(k-2,:) = [yi(k-1), yi(k-2), ui(k-2)];
end


yMQ = yi(3:length(yi));

% Parameter estimation | Estimação dos regressores

for i=1
    teta(:,i) = inv(psi'*psi)*psi'*yMQ; % Without interval
    xi = yMQ - psi*teta(:,1);
    
    tetat(:,i) = inv(psin'*psin)*psin'*yMQ; % Interval approach
end

%% Validação Por Simulação Livre

x(1)=yv(1); %Nominal
x(2)=yv(2);
uval=uv;
h(1) = intval(yv(1));
h(2)=  intval(yv(2));
g(1) = intval(yv(1));
g(2)=  intval(yv(2));

for k=3:length(yv)
    
    % MVF
    xd(k) = gradientinit(g(k-1));
    %
    td(k) = teta(1)*xd(k-1)+teta(2)*xd(k-2)+ teta(3)*uval(k-2);
    %
    g(k)  = teta(1)*mid(g(k-1))+teta(2)*mid(g(k-2))+ teta(3)*uval(k-2) +td(k).dx*(g(k-1)-mid(g(k-1)));
    
    % Nominal
    x(k) = tetat(1)*x(k-1)+tetat(2)*x(k-2)+ tetat(3)*uval(k-2);
    
    % Intervalar
    h(k) = teta(1)*h(k-1)+teta(2)*h(k-2)+ teta(3)*uval(k-2) ;
end

% RMSE tradicional simulação livre
h1=x';
media = mean(yv);
RMSE = sqrt(sum((yv-h1).^2))/sqrt((sum((yv-media).^2)))

% RMSE intervalar tracional
m2=h';
media=mean(yv);
RMSEi= sqrt(sum((yv-m2).^2))/sqrt((sum((yv-media).^2)));

%% Figuras

% Comparação entre dados coletados e modelo nominal
figure
plot(x(5:3000),'--r','LineWidth',1)
hold on
plot(yv(5:300),'-k','LineWidth',1)
set(gca,'FontSize',14,'Fontname','Times New Roman')
xlabel('$k$','Interpreter','LaTex','Fontsize',16)
ylabel('Tens\~ao (V)','Interpreter','LaTex','Fontsize',16)
axis([0 300 0.95 1.01])
xl = get(gca,'YTickLabel');
new_xl = strrep(xl(:),'.',',');
set(gca,'YTickLabel',new_xl)
box off
print(gcf,'-depsc','fig1')

%====================================================================
% Comparação entre dados coletados e modelo nominal, modelo intervalar e nominal

figure
plot(h(5:300),'b')
hold on
plot(x(5:300),'--r','LineWidth',1)
hold on
plot(yv(5:300),'-k','LineWidth',1)
set(gca,'FontSize',14,'Fontname','Times New Roman')
xlabel('$k$','Interpreter','LaTex','Fontsize',16)
ylabel('Tens\~ao (V)','Interpreter','LaTex','Fontsize',16)
axis([0 300 0.94 1.02])
xl = get(gca,'YTickLabel');
new_xl = strrep(xl(:),'.',',');
set(gca,'YTickLabel',new_xl)
box off
print(gcf,'-depsc','fig1')
