function linespec = getlinespec( index )
% GETLINESPEC Get the index-th line specification string syntax.
switch index
	case 1
		linespec = 'b';
	case 2
		linespec = 'r:';
	case 3
		linespec = 'g--';
	case 4
		linespec = 'k-.';
	case 5
		linespec = 'c';
	case 6
		linespec = 'm';
	case 7
		linespec = 'y';
	case 8
		linespec = 'b-.';
	case 9
		linespec = 'r';
	case 10
		linespec = 'g-.';
	case 11
		linespec = 'k--';
	case 12
		linespec = 'c-.';
	case 13
		linespec = 'm-.';
	case 14
		linespec = 'y-.';
	case 15
		linespec = 'b:';
	case 16
		linespec = 'r--';
	case 17
		linespec = 'g';
	case 18
		linespec = 'k:';
	case 19
		linespec = 'c:';
	case 20
		linespec = 'm:';
	case 21
		linespec = 'y:';
	otherwise
		linespec = 'b--';
end
end

