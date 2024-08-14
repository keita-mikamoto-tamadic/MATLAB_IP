clear;

syms t_symbol lamda tau ex_tau
syms theta(t_symbol) phi(t_symbol) fd v
syms ref z

syms g Mp Jp l u_1 u_3 Mw Jw a u_2 u_4 


% 一般化座標ベクトル定義
q = [theta(t_symbol); phi(t_symbol)];
dq = diff(q, t_symbol);
ddq = diff(dq, t_symbol);

% 外力をθ周りのトルクに換算
% ex_torq = fd * cos(q(1)) * 2 * 1;
torq = [ex_tau; tau];


% 振子とタイヤの摩擦トルクとタイヤの転がり抵抗
tau_fp = u_3 * dq(2);
tau_fw = -u_3 * dq(2) - u_4 * (dq(1) + dq(2));
dissipation = [tau_fp, tau_fw];

% 直交座標変換
% タイヤの並進方向
xw = a*(q(1)+q(2));

% 全微分を用いた文字式の常微分
dxw = diff(xw, q(1))*dq(1) + diff(xw, q(2))*dq(2);

% 振子の直交座標表記
xp = l*sin(q(1)) + xw;
yp = l*cos(q(1));

% 振子直交座標の微分
dxp = diff(xp,q(1))*dq(1) + diff(xp,q(2))*dq(2);
dyp = diff(yp, q(1))*dq(1) + diff(yp,q(2))*dq(2);

% ==========================
% 振子のエネルギー
% ==========================
% 並進エネルギー
Kpt = 1/2 * Mp * ((dxp^2) + (dyp^2));

% 回転エネルギー
Kpr = 1/2 * (Jp * dq(1)^2);

% 位置エネルギー(真下が0になるように)
Up = Mp * g * (l + yp);

% 損失エネルギー
Dp = 1/2 * u_1 * dq(1)^2;

% ==========================
% タイヤのエネルギー
% ==========================
% 並進エネルギー
Kwt = 1/2 * Mw * dxw^2;

% 回転エネルギー
Kwr = 1/2 * Jw * (dq(1) + dq(2))^2;

% 位置エネルギー
% 損失エネルギー
Dw = 1/2 * u_2 *dq(2)^2;

% ラグランジュの運動方程式
% エネルギーの総和を求める
K = Kpt + Kpr + Kwt + Kwr;
U = Up;

% 損失エネルギーDをベクトル形式に
D = [Dp; Dw];

% ラグラジアン
L = K-U;

% ============================
% ラグランジェの運動方程式
% ============================
dLq = sym('dLq', [2 1]);
ddLq = sym('ddLq', [2 1]);
eq = sym('eq', [2 1]);

N = 2;

% ラグランジュの運動方程式を立てる
for i = 1:N
    dLq(i) = diff(L,dq(i));
    temp = 0;
    for j = 1:N
        temp = temp + diff(dLq(i),dq(j))*ddq(j) + diff(dLq(i),q(j))*dq(j);
    end
    ddLq(i) = temp;

    eq(i) = ddLq(i) - diff(L, q(i)) + diff(D(i), dq(i)) == torq(i) + dissipation(i);
end
eq_view =eq;

% 状態ベクトル表示にする。
% X1  = dtheta, X2  = dX1, X3 =  dphi, X4 = dX3
% 状態変数
statenum = N;

% 状態ベクトル数
n_state = N * 2;


% 一旦角加速度を単一のシンボリック変数に置く
ddq_temp = sym("ddq_temp", [statenum, 1]);

for i = 1:statenum
    eq(i) = subs(eq(i), ddq, ddq_temp);
end

% 状態について整理
% 各状態量を単一のシンボリック変数に置き換え
for i = 1:statenum
    num = 2 * (i-1) + 1;
    x(num, 1) = q(i);
    x(num+1,1) = dq(i);
end

% 仮の状態ベクトルで置き換え
x_temp = sym("x_temp", [n_state, 1]);
for i = 1:statenum
    eq(i) = subs(eq(i), x(2:2:end), x_temp(2:2:end));
    eq(i) = subs(eq(i), x(1:2:end), x_temp(1:2:end));

end

% 角加速度に対する連立一次方程式を解いて、ddq=の形にする
state_eq_right = solve(eq, ddq_temp);

if statenum == 1
    temp = state_eq_right;
    clear state_eq_right;
    state_eq_right{1} = temp;
else
    state_eq_right = struct2cell(state_eq_right);
end

for i = 1:statenum
    num = 2 * (i - 1) + 1;
    dx_nonline(num,1) = x_temp(num+1);
    dx_nonline(num+1, 1) = simplify(state_eq_right{i});
end

% 非線形運動方程式
matlabFunction(dx_nonline, "File", "state_eq_nonline", "Vars", {g, Mp, Jp, l, u_1, u_3, Mw, Jw, a, u_2, u_4, x_temp, tau, ex_tau})

% ============================
% 状態方程式
% ============================
% AとBを求めるために線形近似
linear_point = [0; 0; 0; 0];
linear_point_in = [0; 0];

jacobi_A = jacobian(dx_nonline, x_temp);
jacobi_B = jacobian(dx_nonline, torq);

% A
A = subs(jacobi_A, x_temp, linear_point);
% B
B = subs(jacobi_B, x_temp, linear_point);
B = subs(B, torq, linear_point_in);
% 外力の項を取り除く
B = B(:, 2);
% C
C = [1, 0, 0, 0];

% パラメータ読み込み
run("param.m");
A = subs(A);
B = subs(B);

% ============================
% 最適サーボゲイン算出
% ============================
% sym型からdouble型に変換
% Aを拡大系へ拡張
A = double(A);
A = vertcat(A, -C);
A = horzcat(A, zeros(5,1));
% Bを拡大系へ拡張
B = double(B);
B = vertcat(B, 0);

Q = double(Q_svr);
R = double(R_lqr);

[K, S, P] = lqr(A, B, Q, R);
Ig = K(5)
K(5) = []
