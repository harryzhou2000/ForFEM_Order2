T = [-270 -260 -240 -220 -200 -180 -160 -140 -120 -100 -80 -60 -40 -20 0 ...
    20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320]+273.15;
E_copper = [139 139 138 138 138 137 136 136 135 134 133 132 131 130 130 129 ...
    128 127 126 125 124 123 122 121 120 119 118 117 116 115 114];

A = [T'.^0,T'.^1,T'.^2,T'.^3,T'.^4];
E_copper_coeffs = A\E_copper';
Ts = linspace(0,600,1000);
Es = polyval(E_copper_coeffs(end:-1:1),Ts);
plot(T,E_copper,Ts,Es);