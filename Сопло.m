T_c = 700;      % Температура внутри (К, по условию)
P_c = 500e3;    % Давление внутри (Па, по условию)
P_amb = 101e3;  % Давление снаружи (Па, атмосферное)
T_amb = 300;    % Температура снаружи (К, стандартное)
gamma = 1.66;   % Гамма
h_th = 0.04;    % Высота сопла (м)

r = 0.01;
num = 100;
noz = zeros(num, 1);    
noz(1,1) = 0.02;        % Верхняя по высоте граница сопла
noz(1,2) = -0.02;       % Нижняя по высоте граница сопла
step = 0.05/num;        

for i = 2 : num*3/5     % Построение первой части сопла
    noz(i, 1) = r * (1 - (i * step - 3 * r)/r)^(1/2);
    noz(i, 2) = - r * (1 - (i * step - 3 * r)/r)^(1/2);
end

for i = num*3/5 : num   % Построение второй части сопла
    noz(i, 1) = r * (1 + (i * step - 3 * r)/r);
    noz(i, 2) = - r * (1 + (i * step - 3 * r)/r);
end

X = 1:num;
Y = num * 7/10 : num;
Z = num * 4/5 : num;

figure(1);          
clf;
plot(X.*step, noz(:,1), 'k', 'LineWidth', 3);

figure(1);
hold on;
plot(X.*step, noz(:,2), 'k', 'LineWidth', 3);
xlabel('Длина сопла (м)');
ylabel('Высота сопла (м)');
axis([0 0.05 -0.04 0.04]);
title('Профиль сопла');
grid on;
hold on;

plot(0, 0.02, '.', 'color', 'black', 'MarkerSize', 15);
hold on;

plot(0, -0.02, '.', 'color', 'black', 'MarkerSize', 15);

A_1 = zeros(num, 1);
mp_sw = zeros(num, 1);
pp_sw = zeros(num, 1);

syms M;
G = (gamma + 1)/(2 * gamma - 2);
g1 = (2/(gamma + 1))^G;
g2 = (1 + ((gamma - 1)/2) * (M^2))^G;
m_1 = zeros(num, 1);
m_2 = zeros(num, 1);

for i = 1 : num * 3/5   % Вычисление чисел Маха и приведённых площадей
                        % для первой части сопла
    A = 100 * (r * (1 - (i * step - 3 * r)/r)^(1/2));
    A_1(i, 1) = A;
    eqn = g1 * g2/M == A;
    m_1(i, 1) = vpasolve(eqn, M, [-inf 1]);
    m_2(i, 1) = m_1(i, 1);
end

m_0 = m_1(1, 1);        % Получение числа Маха в начале сопла

for i = num * 3/5 : num % Вычисление чисел Маха и приведённых площадей
                        % для второй части сопла
    A = 100 * (r * (1 + (i * step - 3 * r)/r));
    A_1(i, 1) = A;
    eqn = g1 * g2/M == A;
    m_1(i, 1) = vpasolve(eqn, M, [1 inf]);
    m_2(i, 1) = vpasolve(eqn, M, [-inf 1]);
end 

m_f1 = m_1(num, 1);     % Получение чисел Маха в конце сопла в
m_f2 = m_2(num, 1);     % сверхзвуковом и дозвуковом случаях

m = ((1 + ((m_f1)^2) * (gamma - 1)/2)/((gamma * (m_f1^2)) - (gamma - 1)/2))^(1/2);
p_21 = ((2/(gamma + 1)) * (1 + (m^2) * (gamma - 1)/ 2) )^(-gamma/(gamma - 1));

figure(2);       
clf;
plot(X.*step, A_1(:,1), 'r', 'LineWidth', 3);
xlabel('Длина сопла (м)');
ylabel('A/A*');
title('Зависимость приведённой площади сечения от расстояния');
grid on;
hold on;

plot(0, 2, '.', 'color', 'red', 'MarkerSize', 15);

P_1 = zeros(num, 2);   
P_2 = zeros(num, 2);

for i = 1 : num * 3/5   % Вычисление приведённых давлений в первой части сопла
    P_1(i, 1) = ((2/(gamma + 1)) * (1 + (m_1(i, 1)^2) * (gamma - 1)/ 2) )^(-gamma/(gamma - 1));
    P_2(i, 1) = ((2/(gamma + 1)) * (1 + (m_2(i, 1)^2) * (gamma - 1)/ 2) )^(-gamma/(gamma - 1));
end

p_0 = P_1(1, 1);        % Получение приведённого давления в начале сопла

for i = num * 3/5 : num % Вычисление приведённых давлений во второй части сопла
    P_1(i, 1) = ((2/(gamma + 1)) * (1 + (m_1(i, 1)^2) * (gamma - 1)/ 2) )^(-gamma/(gamma - 1));
    P_2(i, 1) = ((2/(gamma + 1)) * (1 + (m_2(i, 1)^2) * (gamma - 1)/ 2) )^(-gamma/(gamma - 1));
end 

p_f1 = P_1(num, 1);     % Получение приведённых давлений в конце сопла
p_f2 = P_2(num, 1);     % в сверхзвуковом и дозвуковом случаях

figure(3);
clf;
plot(X.*step, P_1(:,1), 'b', 'LineWidth', 3);
xlabel('Длина сопла (м)');
ylabel('P/P*');
title('Зависимость величины приведённого давления от расстояния');
grid on;
hold on;
plot(X.*step, P_2(:,1), 'c', 'LineWidth', 3);
hold on;

plot(0, p_0, '.', 'color', 'cyan', 'MarkerSize', 15);
text(-0.005, p_0, {'P_0/P*' p_0});
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%% (P_a)' %%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(0.05, p_f1, '.', 'color', 'blue', 'MarkerSize', 15);
text(0.0505, p_f1, {'P_1/P*' p_f1});
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%% (P_a)'' %%%%%%%%%%%%%%%%%%%%%%%%%%
plot(0.05, p_f2, '.', 'color', 'cyan', 'MarkerSize', 15);
text(0.0505, p_f2, {'P_2/P*' p_f2});
hold on;

% plot(0.0351, 0.259476, '.', 'color', 'blue', 'MarkerSize', 15);
% text(0.0351, 0.259476, {'    P/P* = 0.259476'});

% Выбираем случайное значение, расположенное между (P_a)'' 
% и приведённым давлением при М = 0 (P0 в коде программы).
% По графику видно диапазон, пусть для простоты значение = 2

C = (((1/2.02)^((gamma - 1)/gamma) * (gamma + 1)/2 - 1) * 2/(gamma - 1))^(1/2);
mp = zeros(num, 1);
pp = zeros(num, 1);

for i = 1 : num * 3/5         % Вычисление чисел Маха при давлении на выходе
                              % выше (P_a)''/P* 
    A = 100 * (r * (1 - (i * step - 3 * r)/r)^(1/2));
    A_1(i, 1) = A;
    eqn1 = (C/M) * ((1 + (M^2)*(gamma - 1)/2)/(1 + (C^2)*(gamma - 1)/2))^(G) == A/3;
    mp(i, 1) = vpasolve(eqn1, M, [0 5]);
    pp(i, 1) = ((2/(gamma + 1)) * (1 + (mp(i, 1)^2) * (gamma - 1)/ 2) )^(-gamma/(gamma - 1));
end 

for i = num * 3/5 : num    
                              
    A = 100 * (r * (1 + (i * step - 3 * r)/r));
    A_1(i, 1) = A;
    eqn1 = (C/M) * ((1 + (M^2)*(gamma - 1)/2)/(1 + (C^2)*(gamma - 1)/2))^(G) == A/3;
    mp(i, 1) = vpasolve(eqn1, M, [0 5]);
    pp(i, 1) = ((2/(gamma + 1)) * (1 + (mp(i, 1)^2) * (gamma - 1)/ 2) )^(-gamma/(gamma - 1));
end 

plot(X.*step, pp(:, 1), '.', 'color', 'red', 'LineWidth', 3);
hold on;

% for i = 1 : 47
%     plot(0.035, 0.266 + i/30, '.', 'color', 'black', 'MarkerSize', 6);
% end
% hold on;

% plot(0.05, pp(num, 1), '.', 'color', 'red', 'MarkerSize', 15);
% hold on;

plot(0.05, p_21, '.', 'color', 'black', 'MarkerSize', 15);
text(0.0505, p_21, {'P_3/P*' p_21});
hold on;

m_sw = ((1 + ((2.43198)^2) * (gamma - 1)/2)/((gamma * (2.43198^2)) - (gamma - 1)/2))^(1/2);

for i = num * 4/5 : num    
                              
    A = 100 * (r * (1 + (i * step - 3 * r)/r));
    A_1(i, 1) = A;
    eqn1 = (m_sw/M) * ((1 + (M^2)*(gamma - 1)/2)/(1 + (m_sw^2)*(gamma - 1)/2))^(G) == A/1.5;
    mp_sw(i, 1) = vpasolve(eqn1, M, [0 5]);
    pp_sw(i, 1) = ((2/(gamma + 1)) * (1 + (mp_sw(i, 1)^2) * (gamma - 1)/ 2) )^(-gamma/(gamma - 1));
    disp (pp_sw(i, 1));
end

plot(Z.*step, pp_sw(num * 4/5: num, 1), '.', 'color', 'black', 'LineWidth', 3);
hold on;

for i = 1 : 50
    plot(0.04, 0.158 + i/30, '.', 'color', 'black', 'MarkerSize', 6);
end
hold on;

plot(0.04, 0.141468, '.', 'color', 'blue', 'MarkerSize', 15);
hold on;

P0 = (2/(gamma + 1))^(-gamma/(gamma - 1));
plot(X.*step, P0, '.', 'color', 'green', 'LineWidth', 3);
text(0.018, P0 + 0.05, {'P/P* = 2.04883 при М = 0 '});
hold on;

P_us = 2.02;
plot(X.*step, P_us, '.', 'color', 'yellow', 'LineWidth', 3);
text(0.019, P_us - 0.05, {'Выбранное P/P* = 2.02'});
legend({'Сверхзвуковое течение', 'Дозвуковое течение', '', '', '', 'Профиль при давлении на выходе выше P_1', '', '', 'Профиль УВ'}, 'Location', 'southwest');

figure(4);
clf;
plot(X.*step, m_1(:,1),X.*step, m_2(:,1), 'c', 'LineWidth', 3);
axis([0 0.05 0 3]);
xlabel('Длина сопла (м)');
ylabel('M');
title('Зависимость числа Маха потока от расстояния');
grid on;
hold on;

plot(0, m_0, '.', 'color', 'cyan', 'MarkerSize', 15);
text(-0.005, m_0, {'M_0' m_0});
hold on;         


%%%%%%%%%%%%%%%%%%%%%%%%%%% (M_a)' %%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(0.05, m_f1, '.', 'color', 'blue', 'MarkerSize', 15);
text(0.0505, m_f1, {'M_1' m_f1});
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%% (M_a)'' %%%%%%%%%%%%%%%%%%%%%%%%%%
plot(0.05, m_f2, '.', 'color', 'cyan', 'MarkerSize', 15);
text(0.0505, m_f2, {'M_2' m_f2});
hold on;

plot(X.*step, mp(:, 1), '.', 'color', 'red', 'LineWidth', 3);
hold on;

% plot(0.0351, 2, '.', 'color', 'blue', 'MarkerSize', 15);
% text(0.0351, 2, {'    М = 2'});
% hold on;

% for i = 1 : 48
%     plot(0.035, 0.365709 + i/30, '.', 'color', 'black', 'MarkerSize', 6);
% end
% hold on;

% plot(0.05, mp(num, 1), '.', 'color', 'red', 'MarkerSize', 15);
% hold on;

M_3 = ((1 + ((m_f1)^2) * (gamma - 1)/2)/((gamma * (m_f1^2)) - (gamma - 1)/2))^(1/2);
plot(0.05, M_3, '.', 'color', 'black', 'MarkerSize', 15);
text(0.0505, M_3, {'M_3' M_3});
hold on;


m_sw = ((1 + ((2.43198)^2) * (gamma - 1)/2)/((gamma * (2.43198^2)) - (gamma - 1)/2))^(1/2);
p_21 = ((2/(gamma + 1)) * (1 + (m_sw^2) * (gamma - 1)/ 2) )^(-gamma/(gamma - 1));

for i = num * 4/5 : num    
                              
    A = 100 * (r * (1 + (i * step - 3 * r)/r));
    A_1(i, 1) = A;
    eqn1 = (m_sw/M) * ((1 + (M^2)*(gamma - 1)/2)/(1 + (m_sw^2)*(gamma - 1)/2))^(G) == A/1.5;
    mp_sw(i, 1) = vpasolve(eqn1, M, [0 5]);
    pp_sw(i, 1) = ((2/(gamma + 1)) * (1 + (mp_sw(i, 1)^2) * (gamma - 1)/ 2) )^(-gamma/(gamma - 1));
end

plot(Z.*step, mp_sw(num * 4/5: num, 1), '.', 'color', 'black', 'LineWidth', 3);
hold on;

for i = 1 : 60
    plot(0.04, 0.374733 + i/30, '.', 'color', 'black', 'MarkerSize', 6);
end
hold on;

plot(0.04, 2.39588, '.', 'color', 'blue', 'MarkerSize', 15);
legend({'Сверхзвуковое течение', 'Дозвуковое течение', '', '', '', 'Профиль при давлении на выходе выше P_1', '', 'Профиль УВ'}, 'Location', 'northwest');
