clear all
clc

% Importar dados dos arquivos .txt
erro = importdata('erro.txt')';

% Configuração do tempo
T = 0.012; % sample time
t = 0:T:(size(erro, 2)-1)*T;

x = erro(1,1:1000);
y = erro(2,1:1000);
z = erro(3,1:1000);


ise_x = sum(power(x,2))
ise_y = sum(power(y,2))
ise_z = sum(power(z,2))
ise = ise_x + ise_y + ise_z

itae_x = 0;
itae_y = 0;
itae_z = 0;
iae_x = 0;
iae_y = 0;
iae_z = 0;
itse_x = 0;
itse_y = 0;
itse_z = 0;
for k=1:1:1000
    itae_x = itae_x + norm(x(k))*t(k);
    itae_y = itae_y + norm(y(k))*t(k);
    itae_z = itae_z + norm(z(k))*t(k);

    iae_x = iae_x + norm(x(k));
    iae_y = iae_y + norm(y(k));
    iae_z = iae_z + norm(z(k));

    itse_x = itse_x + t(k)*x(k)^2;
    itse_y = itse_y + t(k)*y(k)^2;
    itse_z = itse_z + t(k)*z(k)^2;
    
end
itae_x
itae_y
itae_z

iae_x
iae_y
iae_z

itse_x
itse_y
itse_z

itae = itae_x + itae_y + itae_z
iae = iae_x + iae_y + iae_z
itse = itse_x + itse_y + itse_z

%% Extraçao dos dados de referencia

% Importar dados dos arquivos .txt
erro_ref = importdata('erro_ref.txt')';


x_ref = erro_ref(1,1:1000);
y_ref = erro_ref(2,1:1000);
z_ref = erro_ref(3,1:1000);

ise_x_ref = sum(power(x_ref,2))
ise_y_ref = sum(power(y_ref,2))
ise_z_ref = sum(power(z_ref,2))
ise_ref = ise_x_ref + ise_y_ref + ise_z_ref

itae_x_ref = 0;
itae_y_ref = 0;
itae_z_ref = 0;
iae_x_ref = 0;
iae_y_ref = 0;
iae_z_ref = 0;
itse_x_ref = 0;
itse_y_ref = 0;
itse_z_ref = 0;
for k=1:1:1000
    itae_x_ref = itae_x_ref + norm(x_ref(k))*t(k);
    itae_y_ref = itae_y_ref + norm(y_ref(k))*t(k);
    itae_z_ref = itae_z_ref + norm(z_ref(k))*t(k);

    iae_x_ref = iae_x_ref + norm(x_ref(k));
    iae_y_ref = iae_y_ref + norm(y_ref(k));
    iae_z_ref = iae_z_ref + norm(z_ref(k));

    itse_x_ref = itse_x_ref + t(k)*x_ref(k)^2;
    itse_y_ref = itse_y_ref + t(k)*y_ref(k)^2;
    itse_z_ref = itse_z_ref + t(k)*z_ref(k)^2;
    
end
itae_x_ref
itae_y_ref
itae_z_ref

iae_x_ref
iae_y_ref
iae_z_ref

itse_x_ref
itse_y_ref
itse_z_ref

itae_ref = itae_x_ref + itae_y_ref + itae_z_ref
iae_ref = iae_x_ref + iae_y_ref + iae_z_ref
itse_ref = itse_x_ref + itse_y_ref + itse_z_ref


