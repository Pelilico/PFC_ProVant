%% Código Ajustado com Gráficos Adicionais e Referências - sim 1

close all
clear all
clc

% Importar dados dos arquivos .txt
simulation1 = importdata('in.txt')';
simulation1 = simulation1(:, 1:1000);
inputs1 = importdata('out.txt')';
inputs1 = inputs1(:,1:1000);
ref1 = importdata('ref.txt');
ref1 = ref1(1:1000,:);

% Configuração do tempo
T = 0.012; % sample time
t = 0:T:(size(simulation1, 2)-1)*T;

% Ajustar comprimento das referências, se necessário
if size(ref1, 1) > size(simulation1, 2)
    ref1 = ref1(1:size(simulation1, 2), :); % Truncar se for maior
elseif size(ref1, 1) < size(simulation1, 2)
    ref1 = [ref1; repmat(ref1(end, :), size(simulation1, 2) - size(ref1, 1), 1)]; % Replicar última linha se for menor
end

% Separar componentes da simulação
orientacao = simulation1(1:4, :);

yaw=[];
pitch=[];
roll=[];
x=[];
y=[];
z=[];
vx=[];
vy=[];
vz=[];
p=[];
q=[];
r=[];
for k=1:1:length(simulation1(1,:))
    posicao = mult_qd_casadi([simulation1(5:8, k); zeros(4,1)],inv_qd_casadi([orientacao(:,k); zeros(4,1)]))*2;
    x = [x; posicao(2)]; y = [y; posicao(3)]; z = [z; posicao(4)];

    eulZYX=quat2eul(orientacao(:,k)','zyx');
    yaw=[yaw;eulZYX(1)];
    pitch=[pitch;eulZYX(2)];
    roll=[roll;eulZYX(3)];

    wb = Ad_qd_casadi(Conj_qd_casadi([orientacao(:,k);zeros(4,1)]),[0; simulation1(10:12,k); zeros(4,1)]);
    p = [p; wb(2)]; q = [q; wb(3)]; r = [r; wb(4)];

    velocidade = simulation1(14:16, k) - cross(posicao(2:4), [p(k) q(k) r(k)]);
    vx= [vx; velocidade(1)]; vy= [vy; velocidade(2)]; vz= [vz; velocidade(3)];
end
% Conversão dos ângulos para graus
roll = rad2deg(roll);
pitch = rad2deg(pitch);
yaw = rad2deg(yaw);


% Entradas e forças
u = inputs1;
con = round([1, 1, 1, 1; 0, 0.332, 0, -0.332; -0.332, 0, 0.332, 0; 0.0179, -0.0179, 0.0179, -0.0179] * u, 2);
f1 = u(1, :); f2 = u(2, :); f3 = u(3, :); f4 = u(4, :);

%% Gráficos de Position
figure(1);
subplot(3, 1, 1);
plot(t, x, '-', 'LineWidth', 2); hold on;

ylabel('X (m)'); 
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;


subplot(3, 1, 2);
plot(t, y, '-', 'LineWidth', 2); hold on;

ylabel('Y (m)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;

subplot(3, 1, 3);
plot(t, z, '-', 'LineWidth', 2); hold on;

ylabel('Z (m)'); xlabel('Time (s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;

%% Gráficos de Linear Velocities
figure(2);
subplot(3, 1, 1);
plot(t, vx, '-', 'LineWidth', 2); hold on;

ylabel('Vx (m/s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;


subplot(3, 1, 2);
plot(t, vy, '-', 'LineWidth', 2); hold on;

ylabel('Vy (m/s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;

subplot(3, 1, 3);
plot(t, vz, '-', 'LineWidth', 2); hold on;

ylabel('Vz (m/s)'); xlabel('Time (s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;

%% Gráficos de Orientation (Figura 3)
figure(3);
% Roll
subplot(3, 1, 1);
plot(t, roll, '-', 'LineWidth', 2); hold on;

ylabel('Roll (°)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on;
box off;

% Pitch
subplot(3, 1, 2);
plot(t, pitch, '-', 'LineWidth', 2); hold on;

ylabel('Pitch (°)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on;
box off;

% Yaw
subplot(3, 1, 3);
plot(t, yaw, '-', 'LineWidth', 2); hold on;

ylabel('Yaw (°)');
xlabel('Time (s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on;
box off;
%% Gráficos de Velocidades Angulares
figure(4);
subplot(3, 1, 1);
plot(t, p, '-', 'LineWidth', 2); hold on;grid on;

ylabel('p (rad/s)'); %legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');


subplot(3, 1, 2);
plot(t, q, '-', 'LineWidth', 2); hold on;grid on;

ylabel('q (rad/s)'); %legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');

subplot(3, 1, 3);
plot(t, r, '-', 'LineWidth', 2); hold on; grid on;

ylabel('r (rad/s)'); xlabel('Time (s)'); %legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');

%% Gráficos de Inputs (Sem Referências)
figure(5);


% T (N)
subplot(4, 1, 1);
plot(t, con(1, :), '-', 'LineWidth', 2); grid on; hold on;
ylabel('T (N)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
box off;

% \tau_{\phi} (N.m)
subplot(4, 1, 2);
plot(t, con(2, :), '-', 'LineWidth', 2); grid on; hold on;
ylabel('\tau_{\phi} (N.m)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
box off;

% \tau_{\theta} (N.m)
subplot(4, 1, 3);
plot(t, con(3, :), '-', 'LineWidth', 2); grid on; hold on;
ylabel('\tau_{\theta} (N.m)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
box off;

% \tau_{\psi} (N.m)
subplot(4, 1, 4);
plot(t, con(4, :), '-', 'LineWidth', 2); grid on; hold on;
ylabel('\tau_{\psi} (N.m)');
xlabel('Time (s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
box off;
%% Gráficos de Forces
figure(6);
subplot(4, 1, 1);
plot(t, f1, '-', 'LineWidth', 2); hold on;
%legend('Sim3', 'Sim4',  'Location', 'best');
ylabel('f1 (N)'); grid on;

subplot(4, 1, 2);
plot(t, f2, '-', 'LineWidth', 2); hold on;
%legend('Sim3', 'Sim4',  'Location', 'best');
ylabel('f2 (N)'); grid on;

subplot(4, 1, 3);
plot(t, f3, '-', 'LineWidth', 2); hold on;
%legend('Sim3', 'Sim4',  'Location', 'best');
ylabel('f3 (N)'); 
% xlabel('Time (s)'); grid on;

subplot(4, 1, 4);
plot(t, f4, '-', 'LineWidth', 2); hold on;
%legend('Sim3', 'Sim4',  'Location', 'best');
ylabel('f4 (N)'); xlabel('Time (s)'); grid on;


%% Código Ajustado com Gráficos Adicionais e Referências - sim 1

clear all
clc

% Importar dados dos arquivos .txt
simulation1 = importdata('in_ref.txt')';
simulation1 = simulation1(:, 1:1000);
inputs1 = importdata('out_ref.txt')';
inputs1 = inputs1(:,1:1000);
ref1 = importdata('ref_ref.txt');
ref1 = ref1(1:1000,:);

% Configuração do tempo
T = 0.012; % sample time
t = 0:T:(size(simulation1, 2)-1)*T;

% Ajustar comprimento das referências, se necessário
if size(ref1, 1) > size(simulation1, 2)
    ref1 = ref1(1:size(simulation1, 2), :); % Truncar se for maior
elseif size(ref1, 1) < size(simulation1, 2)
    ref1 = [ref1; repmat(ref1(end, :), size(simulation1, 2) - size(ref1, 1), 1)]; % Replicar última linha se for menor
end

% Separar componentes da simulação
orientacao = simulation1(1:4, :);

yaw=[];
pitch=[];
roll=[];
x=[];
y=[];
z=[];
vx=[];
vy=[];
vz=[];
p=[];
q=[];
r=[];
for k=1:1:length(simulation1(1,:))
    posicao = mult_qd_casadi([simulation1(5:8, k); zeros(4,1)],inv_qd_casadi([orientacao(:,k); zeros(4,1)]))*2;
    x = [x; posicao(2)]; y = [y; posicao(3)]; z = [z; posicao(4)];

    eulZYX=quat2eul(orientacao(:,k)','zyx');
    yaw=[yaw;eulZYX(1)];
    pitch=[pitch;eulZYX(2)];
    roll=[roll;eulZYX(3)];

    wb = Ad_qd_casadi(Conj_qd_casadi([orientacao(:,k);zeros(4,1)]),[0; simulation1(10:12,k); zeros(4,1)]);
    p = [p; wb(2)]; q = [q; wb(3)]; r = [r; wb(4)];

    velocidade = simulation1(14:16, k) - cross(posicao(2:4), [p(k) q(k) r(k)]);
    vx= [vx; velocidade(1)]; vy= [vy; velocidade(2)]; vz= [vz; velocidade(3)];
end
% Conversão dos ângulos para graus
roll = rad2deg(roll);
pitch = rad2deg(pitch);
yaw = rad2deg(yaw);



% Separar componentes das referências
orientacao_d = ref1(:, 1:4);

yaw_d=[];
pitch_d=[];
roll_d=[];
xref1=[];
yref1=[];
zref1=[];
vxref1=[];
vyref1=[];
vzref1=[];
ohmx_ref1=[];
ohmy_ref1=[];
ohmz_ref1=[];
for k=1:1:length(ref1(:,1))
    posicao_d = mult_qd_casadi([ref1(k ,5:8), zeros(1,4)]',inv_qd_casadi([orientacao_d(k,:), zeros(1,4)]'))*2;
    xref1 = [xref1; posicao_d(2)]; yref1 = [yref1; posicao_d(3)]; zref1 = [zref1; posicao_d(4)];

    eulZYX=quat2eul(orientacao_d(k,:),'zyx');
    yaw_d=[yaw_d;eulZYX(1)];
    pitch_d=[pitch_d;eulZYX(2)];
    roll_d=[roll_d;eulZYX(3)];

    wb = Ad_qd_casadi(Conj_qd_casadi([orientacao_d(k,:),zeros(1,4)]'),[0, ref1(k,10:12), zeros(1,4)]);
    ohmx_ref1 = [ohmx_ref1; wb(2)]; ohmy_ref1 = [ohmy_ref1; wb(3)]; ohmz_ref1 = [ohmz_ref1; wb(4)];

    velocidade = ref1(k, 14:16)' - cross(posicao_d(2:4), [ohmx_ref1(k) ohmy_ref1(k) ohmz_ref1(k)]);
    vxref1 = [vxref1; velocidade(1)]; vyref1 = [vyref1; velocidade(2)]; vzref1 = [vzref1; velocidade(3)];
end
% Conversão dos ângulos para graus
roll_d = rad2deg(roll_d);
pitch_d = rad2deg(pitch_d);
yaw_d = rad2deg(yaw_d);




% Entradas e forças
u = inputs1;
con = round([1, 1, 1, 1; 0, 0.332, 0, -0.332; -0.332, 0, 0.332, 0; 0.0179, -0.0179, 0.0179, -0.0179] * u, 2);
f1 = u(1, :); f2 = u(2, :); f3 = u(3, :); f4 = u(4, :);

%% Gráficos de Position
figure(1);
subplot(3, 1, 1);
plot(t, x, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
ylabel('X (m)'); 
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;


subplot(3, 1, 2);
plot(t, y, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
ylabel('Y (m)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;
legend('FNMPC','NMPC')
subplot(3, 1, 3);
plot(t, z, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC');
ylabel('Z (m)'); xlabel('Time (s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;

%% Gráficos de Linear Velocities
figure(2);
subplot(3, 1, 1);
plot(t, vx, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
ylabel('Vx (m/s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;


subplot(3, 1, 2);
plot(t, vy, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
ylabel('Vy (m/s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;

subplot(3, 1, 3);
plot(t, vz, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
ylabel('Vz (m/s)'); xlabel('Time (s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on; box off;

%% Gráficos de Orientation (Figura 3)
figure(3);
% Roll
subplot(3, 1, 1);
plot(t, roll, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
ylabel('Roll (°)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on;
box off;

% Pitch
subplot(3, 1, 2);
plot(t, pitch, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
ylabel('Pitch (°)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on;
box off;

% Yaw
subplot(3, 1, 3);
plot(t, yaw, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
ylabel('Yaw (°)');
xlabel('Time (s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
grid on;
box off;
%% Gráficos de Velocidades Angulares
figure(4);
subplot(3, 1, 1);
plot(t, p, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
ylabel('p (rad/s)'); %legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');


subplot(3, 1, 2);
plot(t, q, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
ylabel('q (rad/s)'); %legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');

subplot(3, 1, 3);
plot(t, r, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
ylabel('r (rad/s)'); xlabel('Time (s)'); %legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');

%% Gráficos de Inputs (Sem Referências)
figure(5);
% T (N)
subplot(4, 1, 1);
plot(t, con(1, :), '--', 'LineWidth', 2); grid on;
legend('FNMPC','NMPC')
ylabel('T (N)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');

% \tau_{\phi} (N.m)
subplot(4, 1, 2);
plot(t, con(2, :), '--', 'LineWidth', 2); grid on;
legend('FNMPC','NMPC')
ylabel('\tau_{\phi} (N.m)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');

% \tau_{\theta} (N.m)
subplot(4, 1, 3);
plot(t, con(3, :), '--', 'LineWidth', 2); grid on;
legend('FNMPC','NMPC')
ylabel('\tau_{\theta} (N.m)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');

% \tau_{\psi} (N.m)
subplot(4, 1, 4);
plot(t, con(4, :), '--', 'LineWidth', 2); grid on;
legend('FNMPC','NMPC')
ylabel('\tau_{\psi} (N.m)');
xlabel('Time (s)');
%legend('Sim3', 'Ref', 'Sim4',  'Location', 'best');
%% Gráficos de Forces
figure(6);
subplot(4, 1, 1);
plot(t, f1, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
%legend('Sim3', 'Sim4',  'Location', 'best');
ylabel('f1 (N)'); grid on;

subplot(4, 1, 2);
plot(t, f2, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
%legend('Sim3', 'Sim4',  'Location', 'best');
ylabel('f2 (N)'); grid on;

subplot(4, 1, 3);
plot(t, f3, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
%legend('Sim3', 'Sim4',  'Location', 'best');
ylabel('f3 (N)'); grid on;
% xlabel('Time (s)'); grid on;

subplot(4, 1, 4);
plot(t, f4, '--', 'LineWidth', 2); hold on;
legend('FNMPC','NMPC')
%legend('Sim3', 'Sim4',  'Location', 'best');
ylabel('f4 (N)'); xlabel('Time (s)'); grid on;
