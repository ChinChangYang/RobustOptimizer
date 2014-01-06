function displayitermessages(X, U, f, countiter, XX, YY, ZZ, varargin)
%DISPLAYITERMESSAGES Display iterative messages
persistent fPrev
stdX = std(X, 0, 2);
[maxstdX, maxstdXidx] = max(stdX);
[minstdX, minstdXidx] = min(stdX);
fprintf(...
	['Iter: %d,\tFmin: %0.4E,\tFmean: %0.4E,\t', ...
	'Fstd: %0.4E,\t', ...
	'Max Std X(%d):%0.4E,\t', ...
	'Min Std X(%d):%0.4E\n'], ...
	countiter, min(f(:)), mean(f(:)), std(f(:)), maxstdXidx, maxstdX, ...
	minstdXidx, minstdX);

if ~isempty(fPrev) && isequal(size(fPrev), size(numel(f)))
	pr = sum(fPrev ~= f) / numel(fPrev);
	fprintf('Progressive rate = %f%%\n', pr * 100);
end

fPrev = f;

if nargin >= 7
	[D, ~] = size(X);
	[~, ~, cD] = size(XX);
	
	if D >= 2 && cD >= 2
		for i = 1 : cD
			subplot(1, cD, i);
			hold off;
			contour(XX(:, :, i), YY(:, :, i), ZZ(:, :, i), 7);
			hold on;
			plot(X(i,:), X(i+1,:), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
			plot(U(i,:), U(i+1,:), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
			pause(0.01);
		end
	elseif D >= 2 && cD == 1
		hold off;
		contour(XX, YY, ZZ, 7);
		hold on;
		plot(X(1,:), X(2,:), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
		plot(U(1,:), U(2,:), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
		pause(0.01);
	else
		hold off;
		plot(XX, YY);
		hold on;
		plot(X, f, 'kx');
		plot(U, zeros(1, numel(U)), 'rx');
		pause(0.01);
	end
end

if ~isempty(varargin)
	nvarargin = numel(varargin);
	for i = 1 : floor(nvarargin / 2)
		fprintf('%s = %.4E\n', varargin{2*i-1}, varargin{2*i});
	end
end

end
