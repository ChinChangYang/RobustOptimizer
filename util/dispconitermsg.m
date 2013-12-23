function dispconitermsg(X, U, f, countiter, XX, YY, ZZ, CC, varargin)
%DISPCONSITERMESSAGES Display constraint iteration messages
persistent fPrev

hold off;
plot(0, 0, 'w');
if ~isempty(CC)
	[D, ~] = size(X);
	if D >= 2
		contour(XX, YY, CC);
	end
end

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
	if D == 2
		hold on;
		contour(XX, YY, ZZ, 7);
		plot(X(1,:), X(2,:), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
		plot(U(1,:), U(2,:), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
		xlabel('x');
		ylabel('y');
		title(sprintf('Generation %d', countiter));
		axis([min(XX(:)), max(XX(:)), min(YY(:)), max(YY(:))]);
	elseif D >= 4
		subplot(121);
		hold off;
		contour(XX, YY, ZZ, 7);
		hold on;
		plot(X(1,:), X(2,:), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
		plot(U(1,:), U(2,:), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
		xlabel('x1');
		ylabel('x2');
		title(sprintf('Generation %d', countiter));
		axis([min(XX(:)), max(XX(:)), min(YY(:)), max(YY(:))]);
		subplot(122);		
		hold off;
		plot(X(3,:), X(4,:), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
		hold on;
		plot(U(3,:), U(4,:), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
		xlabel('y1');
		ylabel('y2');
		title(sprintf('Generation %d', countiter));
	else
		hold on;
		plot(XX, YY);
		plot(X, f, 'kx');
		plot(U, zeros(1, numel(U)), 'rx');
	end
	
	pause(0.01);
end

if ~isempty(varargin)
	nvarargin = numel(varargin);
	for i = 1 : floor(nvarargin / 2)
		fprintf('%s = %.4E\n', varargin{2*i-1}, varargin{2*i});
	end
end

end
