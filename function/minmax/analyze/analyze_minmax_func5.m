%% Two-Two dimension
syms x1 x2 y1 y2
f = x1^2 + x2^2 - ...
	((1 + x1 - y1)^2 + 100 * (y2 - x2 - (y1 - x1)^2)^2);
df_x1 = diff(f, x1);
df_x2 = diff(f, x2);
df_y1 = diff(f, y1);
df_y2 = diff(f, y2);
eq1 = df_x1 == 0;
eq2 = df_x2 == 0;
eq3 = df_y1 == 0;
eq4 = df_y2 == 0;
sol = solve(eq1, eq2, eq3, eq4);
fprintf('sol.x1 = %f\n', double(sol.x1));
fprintf('sol.x2 = %f\n', double(sol.x2));
fprintf('sol.y1 = %f\n', double(sol.y1));
fprintf('sol.y2 = %f\n', double(sol.y2));

%% Three-Three dimension
clear
syms x1 x2 x3 y1 y2 y3
f = x1^2 + x2^2 + x3^2 - ...
	(1 + x1 - y1)^2 + 100 * (y2 - x2 - (y1 - x1)^2)^2 - ...
	(1 + x2 - y2)^2 + 100 * (y3 - x3 - (y2 - x2)^2)^2;
df_x1 = diff(f, x1);
df_x2 = diff(f, x2);
df_x3 = diff(f, x3);
df_y1 = diff(f, y1);
df_y2 = diff(f, y2);
df_y3 = diff(f, y3);
eq1 = df_x1 == 0;
eq2 = df_x2 == 0;
eq3 = df_x3 == 0;
eq4 = df_y1 == 0;
eq5 = df_y2 == 0;
eq6 = df_y3 == 0;
sol = solve(eq1, eq2, eq3, eq4, eq5, eq6);
fprintf('sol.x1 = \n');
disp(sol.x1);
fprintf('sol.x2 = \n');
disp(sol.x2);
fprintf('sol.x3 = \n');
disp(sol.x3);
fprintf('sol.y1 = \n');
disp(sol.y1);
fprintf('sol.y2 = \n');
disp(sol.y2);
fprintf('sol.y3 = \n');
disp(sol.y3);
