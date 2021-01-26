% Algoritmo: Symbolic GNSS Atttitude Determination
% Autor: Felipe Oliveira e Silva
% Instituto: DAT/UFLA

clc, clear, close all, format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% triedro i = ECI
% triedro e = ECEF
% triedro n = ENU
% triedro l = NED
% triedo b = BODY
% _X_Y = receiver x, satellite Y
% _Xx_Yy = receiver x w.r.t. receiver X, satellite y w.r.t. satellite Y

%%%%%%%%%%%%%%%%%%%%%%%%% Dual-antenna Baseline %%%%%%%%%%%%%%%%%%%%%%%%%%%

syms rN rE rD real
syms a b c real

r2l_12 = [rN;rE;rD];
r2b_12 = [a;b;c];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gravity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms gp positive
syms ax ay az real

gp2l = [0;0;gp];                                                            % Jiang(1998) (1) / Rogers (11.6)
gp2b = -[ax;ay;az];                                                         % Jiang (10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRIAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = [gp2l r2l_12  cross(gp2l,r2l_12)];                                      % Zhao (1)
B = [gp2b r2b_12  cross(gp2b,r2b_12)];                                      % Zhao (1)

% Cb2l = L*inv(B);                                                            % Zhao (3)
Cb2l = transpose(inv(L))*transpose(B);                                      % Zhao (3)

Cb2l = simplify(Cb2l);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Error analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms dax day daz da db dc dgp real

dCb2l = diff(Cb2l,ax)*dax + diff(Cb2l,ay)*day + diff(Cb2l,az)*daz + diff(Cb2l,a)*da + diff(Cb2l,b)*db + diff(Cb2l,c)*dc + diff(Cb2l,gp)*dgp;  % Savage (3.5.2-1)

E = dCb2l*transpose(Cb2l);                                                  % Savage (3.5.1-2) - adaptado

%%%%%%%%%%%%%%%%%%%%%%%% Particular case Cb2l = I %%%%%%%%%%%%%%%%%%%%%%%%%

% E = subs(E,{ax,ay,az,rN,rE,rD},{0,0,-gp,a,b,c});
% E = subs(E,{ax,ay,az,a,b,c},{0,0,-gp,rN,rE,rD});

%%%%%%%%%%%%%%%%%%%%%%%%% Generic case Cb2l ~= I %%%%%%%%%%%%%%%%%%%%%%%%%%

E_sym = (E + transpose(E))/2;                                               % Savage (3.5.1-6)
E_sksym = (E - transpose(E))/2;                                             % Savage (3.5.1-6)

varphi_skew = -E_sksym;                                                     % Jiang (19)

eta_n = E_sym(1,1);                                                         % Jiang (41) - erros de normalidade
eta_e = E_sym(2,2);                                                         % Jiang (42)
eta_d = E_sym(3,3);                                                         % Jiang (43)
o_n = E_sym(2,3);                                                           % Jiang (40) - erros de ortogonalidade
o_e = E_sym(1,3);                                                           % Jiang (40)
o_d = E_sym(1,2);                                                           % Jiang (40)
varphi_n = E_sksym(2,3);                                                    % Jiang (37) - erros de alinhamento
varphi_e = -E_sksym(1,3);                                                   % Jiang (38)
varphi_d = E_sksym(1,2);                                                    % Jiang (39)

% Analyis as a function of daN daE daD drN drE drD

syms daN daE daD drN drE drD real

da2b = inv(Cb2l)*[daN;daE;daD];
dr2b = inv(Cb2l)*[drN;drE;drD];

eta_n = simplify(subs(eta_n,{dax,day,daz,da,db,dc},{da2b(1),da2b(2),da2b(3),dr2b(1),dr2b(2),dr2b(3)}));
eta_e = simplify(subs(eta_e,{dax,day,daz,da,db,dc},{da2b(1),da2b(2),da2b(3),dr2b(1),dr2b(2),dr2b(3)}));
eta_d = simplify(subs(eta_d,{dax,day,daz,da,db,dc},{da2b(1),da2b(2),da2b(3),dr2b(1),dr2b(2),dr2b(3)}));
o_n = simplify(subs(o_n,{dax,day,daz,da,db,dc},{da2b(1),da2b(2),da2b(3),dr2b(1),dr2b(2),dr2b(3)}));
o_e = simplify(subs(o_e,{dax,day,daz,da,db,dc},{da2b(1),da2b(2),da2b(3),dr2b(1),dr2b(2),dr2b(3)}));
o_d = simplify(subs(o_d,{dax,day,daz,da,db,dc},{da2b(1),da2b(2),da2b(3),dr2b(1),dr2b(2),dr2b(3)}));
varphi_n = simplify(subs(varphi_n,{dax,day,daz,da,db,dc},{da2b(1),da2b(2),da2b(3),dr2b(1),dr2b(2),dr2b(3)}));
varphi_e = simplify(subs(varphi_e,{dax,day,daz,da,db,dc},{da2b(1),da2b(2),da2b(3),dr2b(1),dr2b(2),dr2b(3)}));
varphi_d = simplify(subs(varphi_d,{dax,day,daz,da,db,dc},{da2b(1),da2b(2),da2b(3),dr2b(1),dr2b(2),dr2b(3)}));

% pretty(simplify(eta_n))
% pretty(simplify(eta_e))
% pretty(simplify(eta_d))
% pretty(simplify(o_n))
% pretty(simplify(o_e))
% pretty(simplify(o_d))
% pretty(simplify(varphi_n))
% pretty(simplify(varphi_e))
% pretty(simplify(varphi_d))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Detailed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = eta_n;
% x = eta_e;
% x = eta_d;
% x = o_n;
% x = o_e; 
% x = o_d;
% x = varphi_n;
% x = varphi_e;
% x = varphi_d;

% dxdax = collect(simplify(subs(x,{day,daz,da,db,dc,dgp},{0,0,0,0,0,0})),dax); % termos de 'x' dependentes de dax
% dxday = collect(simplify(subs(x,{dax,daz,da,db,dc,dgp},{0,0,0,0,0,0})),day); % termos de 'x' dependentes de day
% dxdaz = collect(simplify(subs(x,{dax,day,da,db,dc,dgp},{0,0,0,0,0,0})),daz); % termos de 'x' dependentes de daz
% dxda = collect(simplify(subs(x,{dax,day,daz,db,dc,dgp},{0,0,0,0,0,0})),da); % termos de 'x' dependentes de da
% dxdb = collect(simplify(subs(x,{dax,day,daz,da,dc,dgp},{0,0,0,0,0,0})),db); % termos de 'x' dependentes de db
% dxdc = collect(simplify(subs(x,{dax,day,daz,da,db,dgp},{0,0,0,0,0,0})),dc); % termos de 'x' dependentes de dc
% dxdgp = collect(simplify(subs(x,{dax,day,daz,da,db,dc},{0,0,0,0,0,0})),dgp); % termos de 'x' dependentes de dgp

% pretty(simplify(dxdax))
% pretty(simplify(dxday))
% pretty(simplify(dxdaz))
% pretty(simplify(dxda))
% pretty(simplify(dxdb))
% pretty(simplify(dxdc))
% pretty(simplify(dxdgp))

% eta_n = dxdax + dxday + dxdaz + dxda + dxdb + dxdc + dxdgp;
% eta_e = dxdax + dxday + dxdaz + dxda + dxdb + dxdc + dxdgp;
% eta_d = dxdax + dxday + dxdaz + dxda + dxdb + dxdc + dxdgp;
% o_n = dxdax + dxday + dxdaz + dxda + dxdb + dxdc + dxdgp;
% o_e = dxdax + dxday + dxdaz + dxda + dxdb + dxdc + dxdgp;
% o_d = dxdax + dxday + dxdaz + dxda + dxdb + dxdc + dxdgp;
% varphi_n = dxdax + dxday + dxdaz + dxda + dxdb + dxdc + dxdgp;
% varphi_e = dxdax + dxday + dxdaz + dxda + dxdb + dxdc + dxdgp;
% varphi_d = dxdax + dxday + dxdaz + dxda + dxdb + dxdc + dxdgp;

dxdaN = collect(simplify(subs(x,{daE,daD,drN,drE,drD,dgp},{0,0,0,0,0,0})),daN); % termos de 'x' dependentes de daN
dxdaE = collect(simplify(subs(x,{daN,daD,drN,drE,drD,dgp},{0,0,0,0,0,0})),daE); % termos de 'x' dependentes de daE
dxdaD = collect(simplify(subs(x,{daN,daE,drN,drE,drD,dgp},{0,0,0,0,0,0})),daD); % termos de 'x' dependentes de daD
dxdrN = collect(simplify(subs(x,{daN,daE,daD,drE,drD,dgp},{0,0,0,0,0,0})),drN); % termos de 'x' dependentes de drN
dxdrE = collect(simplify(subs(x,{daN,daE,daD,drN,drD,dgp},{0,0,0,0,0,0})),drE); % termos de 'x' dependentes de drE
dxdrD = collect(simplify(subs(x,{daN,daE,daD,drN,drE,dgp},{0,0,0,0,0,0})),drD); % termos de 'x' dependentes de drD
dxdgp = collect(simplify(subs(x,{daN,daE,daD,drN,drE,drD},{0,0,0,0,0,0})),dgp); % termos de 'x' dependentes de dgp

dxdaN = subs(dxdaN,{ax^2 + ay^2 + az^2},{gp^2});
dxdaE = subs(dxdaE,{ax^2 + ay^2 + az^2},{gp^2});
dxdaD = subs(dxdaD,{ax^2 + ay^2 + az^2},{gp^2});
dxdrN = subs(dxdrN,{ax^2 + ay^2 + az^2},{gp^2});
dxdrE = subs(dxdrE,{ax^2 + ay^2 + az^2},{gp^2});
dxdrD = subs(dxdrD,{ax^2 + ay^2 + az^2},{gp^2});
dxdgp = subs(dxdgp,{ax^2 + ay^2 + az^2},{gp^2});

dxdaN = subs(dxdaN,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdaE = subs(dxdaE,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdaD = subs(dxdaD,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdrN = subs(dxdrN,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdrE = subs(dxdrE,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdrD = subs(dxdrD,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdgp = subs(dxdgp,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});

dxdaN = subs(dxdaN,{a*ax + b*ay + c*az},{-gp*rD});
dxdaE = subs(dxdaE,{a*ax + b*ay + c*az},{-gp*rD});
dxdaD = subs(dxdaD,{a*ax + b*ay + c*az},{-gp*rD});
dxdrN = subs(dxdrN,{a*ax + b*ay + c*az},{-gp*rD});
dxdrE = subs(dxdrE,{a*ax + b*ay + c*az},{-gp*rD});
dxdrD = subs(dxdrD,{a*ax + b*ay + c*az},{-gp*rD});
dxdgp = subs(dxdgp,{a*ax + b*ay + c*az},{-gp*rD});

dxdaN = subs(dxdaN,{a*ax*rE^2 - gp*rD*rN^2 + ay*b*rE^2 + az*c*rE^2},{-gp*rD*(rN^2 + rE^2)});
dxdaE = subs(dxdaE,{a*ax*rE^2 - gp*rD*rN^2 + ay*b*rE^2 + az*c*rE^2},{-gp*rD*(rN^2 + rE^2)});
dxdaD = subs(dxdaD,{a*ax*rE^2 - gp*rD*rN^2 + ay*b*rE^2 + az*c*rE^2},{-gp*rD*(rN^2 + rE^2)});
dxdrN = subs(dxdrN,{a*ax*rE^2 - gp*rD*rN^2 + ay*b*rE^2 + az*c*rE^2},{-gp*rD*(rN^2 + rE^2)});
dxdrE = subs(dxdrE,{a*ax*rE^2 - gp*rD*rN^2 + ay*b*rE^2 + az*c*rE^2},{-gp*rD*(rN^2 + rE^2)});
dxdrD = subs(dxdrD,{a*ax*rE^2 - gp*rD*rN^2 + ay*b*rE^2 + az*c*rE^2},{-gp*rD*(rN^2 + rE^2)});
dxdgp = subs(dxdgp,{a*ax*rE^2 - gp*rD*rN^2 + ay*b*rE^2 + az*c*rE^2},{-gp*rD*(rN^2 + rE^2)});

dxdaN = subs(dxdaN,{a*ax*rN^2 - gp*rD*rE^2 + ay*b*rN^2 + az*c*rN^2},{-gp*rD*(rN^2 + rE^2)});
dxdaE = subs(dxdaE,{a*ax*rN^2 - gp*rD*rE^2 + ay*b*rN^2 + az*c*rN^2},{-gp*rD*(rN^2 + rE^2)});
dxdaD = subs(dxdaD,{a*ax*rN^2 - gp*rD*rE^2 + ay*b*rN^2 + az*c*rN^2},{-gp*rD*(rN^2 + rE^2)});
dxdrN = subs(dxdrN,{a*ax*rN^2 - gp*rD*rE^2 + ay*b*rN^2 + az*c*rN^2},{-gp*rD*(rN^2 + rE^2)});
dxdrE = subs(dxdrE,{a*ax*rN^2 - gp*rD*rE^2 + ay*b*rN^2 + az*c*rN^2},{-gp*rD*(rN^2 + rE^2)});
dxdrD = subs(dxdrD,{a*ax*rN^2 - gp*rD*rE^2 + ay*b*rN^2 + az*c*rN^2},{-gp*rD*(rN^2 + rE^2)});
dxdgp = subs(dxdgp,{a*ax*rN^2 - gp*rD*rE^2 + ay*b*rN^2 + az*c*rN^2},{-gp*rD*(rN^2 + rE^2)});

dxdaN = subs(dxdaN,{gp*a^2 + ax*rD*a + gp*b^2 + ay*rD*b + gp*c^2 + az*rD*c},{gp*(rE^2 + rN^2)});
dxdaE = subs(dxdaE,{gp*a^2 + ax*rD*a + gp*b^2 + ay*rD*b + gp*c^2 + az*rD*c},{gp*(rE^2 + rN^2)});
dxdaD = subs(dxdaD,{gp*a^2 + ax*rD*a + gp*b^2 + ay*rD*b + gp*c^2 + az*rD*c},{gp*(rE^2 + rN^2)});
dxdrN = subs(dxdrN,{gp*a^2 + ax*rD*a + gp*b^2 + ay*rD*b + gp*c^2 + az*rD*c},{gp*(rE^2 + rN^2)});
dxdrE = subs(dxdrE,{gp*a^2 + ax*rD*a + gp*b^2 + ay*rD*b + gp*c^2 + az*rD*c},{gp*(rE^2 + rN^2)});
dxdrD = subs(dxdrD,{gp*a^2 + ax*rD*a + gp*b^2 + ay*rD*b + gp*c^2 + az*rD*c},{gp*(rE^2 + rN^2)});
dxdgp = subs(dxdgp,{gp*a^2 + ax*rD*a + gp*b^2 + ay*rD*b + gp*c^2 + az*rD*c},{gp*(rE^2 + rN^2)});

dxdaN = subs(dxdaN,{ax^2*rE^2 + ay^2*rE^2 + az^2*rE^2 + gp^2*rN^2},{gp^2*(rE^2 + rN^2)});
dxdaE = subs(dxdaE,{ax^2*rE^2 + ay^2*rE^2 + az^2*rE^2 + gp^2*rN^2},{gp^2*(rE^2 + rN^2)});
dxdaD = subs(dxdaD,{ax^2*rE^2 + ay^2*rE^2 + az^2*rE^2 + gp^2*rN^2},{gp^2*(rE^2 + rN^2)});
dxdrN = subs(dxdrN,{ax^2*rE^2 + ay^2*rE^2 + az^2*rE^2 + gp^2*rN^2},{gp^2*(rE^2 + rN^2)});
dxdrE = subs(dxdrE,{ax^2*rE^2 + ay^2*rE^2 + az^2*rE^2 + gp^2*rN^2},{gp^2*(rE^2 + rN^2)});
dxdrD = subs(dxdrD,{ax^2*rE^2 + ay^2*rE^2 + az^2*rE^2 + gp^2*rN^2},{gp^2*(rE^2 + rN^2)});
dxdgp = subs(dxdgp,{ax^2*rE^2 + ay^2*rE^2 + az^2*rE^2 + gp^2*rN^2},{gp^2*(rE^2 + rN^2)});

dxdaN = subs(dxdaN,{ax^2*rN^2 + ay^2*rN^2 + az^2*rN^2 + gp^2*rE^2},{gp^2*(rE^2 + rN^2)});
dxdaE = subs(dxdaE,{ax^2*rN^2 + ay^2*rN^2 + az^2*rN^2 + gp^2*rE^2},{gp^2*(rE^2 + rN^2)});
dxdaD = subs(dxdaD,{ax^2*rN^2 + ay^2*rN^2 + az^2*rN^2 + gp^2*rE^2},{gp^2*(rE^2 + rN^2)});
dxdrN = subs(dxdrN,{ax^2*rN^2 + ay^2*rN^2 + az^2*rN^2 + gp^2*rE^2},{gp^2*(rE^2 + rN^2)});
dxdrE = subs(dxdrE,{ax^2*rN^2 + ay^2*rN^2 + az^2*rN^2 + gp^2*rE^2},{gp^2*(rE^2 + rN^2)});
dxdrD = subs(dxdrD,{ax^2*rN^2 + ay^2*rN^2 + az^2*rN^2 + gp^2*rE^2},{gp^2*(rE^2 + rN^2)});
dxdgp = subs(dxdgp,{ax^2*rN^2 + ay^2*rN^2 + az^2*rN^2 + gp^2*rE^2},{gp^2*(rE^2 + rN^2)});

dxdaN = subs(dxdaN,{rD*ax^2 + a*gp*ax + rD*ay^2 + b*gp*ay + rD*az^2 + c*gp*az},{0});
dxdaE = subs(dxdaE,{rD*ax^2 + a*gp*ax + rD*ay^2 + b*gp*ay + rD*az^2 + c*gp*az},{0});
dxdaD = subs(dxdaD,{rD*ax^2 + a*gp*ax + rD*ay^2 + b*gp*ay + rD*az^2 + c*gp*az},{0});
dxdrN = subs(dxdrN,{rD*ax^2 + a*gp*ax + rD*ay^2 + b*gp*ay + rD*az^2 + c*gp*az},{0});
dxdrE = subs(dxdrE,{rD*ax^2 + a*gp*ax + rD*ay^2 + b*gp*ay + rD*az^2 + c*gp*az},{0});
dxdrD = subs(dxdrD,{rD*ax^2 + a*gp*ax + rD*ay^2 + b*gp*ay + rD*az^2 + c*gp*az},{0});
dxdgp = subs(dxdgp,{rD*ax^2 + a*gp*ax + rD*ay^2 + b*gp*ay + rD*az^2 + c*gp*az},{0});

dxdaN = subs(dxdaN,{ax^2 + ay^2 + az^2},{gp^2});
dxdaE = subs(dxdaE,{ax^2 + ay^2 + az^2},{gp^2});
dxdaD = subs(dxdaD,{ax^2 + ay^2 + az^2},{gp^2});
dxdrN = subs(dxdrN,{ax^2 + ay^2 + az^2},{gp^2});
dxdrE = subs(dxdrE,{ax^2 + ay^2 + az^2},{gp^2});
dxdrD = subs(dxdrD,{ax^2 + ay^2 + az^2},{gp^2});
dxdgp = subs(dxdgp,{ax^2 + ay^2 + az^2},{gp^2});

dxdaN = subs(dxdaN,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdaE = subs(dxdaE,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdaD = subs(dxdaD,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdrN = subs(dxdrN,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdrE = subs(dxdrE,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdrD = subs(dxdrD,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});
dxdgp = subs(dxdgp,{a^2 + b^2 + c^2},{rN^2 + rE^2 + rD^2});

dxdaN = subs(dxdaN,{a*ax + b*ay + c*az},{-gp*rD});
dxdaE = subs(dxdaE,{a*ax + b*ay + c*az},{-gp*rD});
dxdaD = subs(dxdaD,{a*ax + b*ay + c*az},{-gp*rD});
dxdrN = subs(dxdrN,{a*ax + b*ay + c*az},{-gp*rD});
dxdrE = subs(dxdrE,{a*ax + b*ay + c*az},{-gp*rD});
dxdrD = subs(dxdrD,{a*ax + b*ay + c*az},{-gp*rD});
dxdgp = subs(dxdgp,{a*ax + b*ay + c*az},{-gp*rD});

pretty(simplify(dxdaN))
pretty(simplify(dxdaE))
pretty(simplify(dxdaD))
pretty(simplify(dxdrN))
pretty(simplify(dxdrE))
pretty(simplify(dxdrD))
pretty(simplify(dxdgp))

% Analysis as a function of da, db and dc

dr2l = Cb2l*[da;db;dc];

dxdrNabc = subs(dxdrN,drN,dr2l(1));
dxdrEabc = subs(dxdrE,drE,dr2l(2));
dxdrDabc = subs(dxdrD,drD,dr2l(3));

dxdrNa = collect(simplify(subs(dxdrNabc,{db,dc},{0,0})),da);                % termos de 'dxdrN' dependentes de da
dxdrNb = collect(simplify(subs(dxdrNabc,{da,dc},{0,0})),db);                % termos de 'dxdrN' dependentes de db
dxdrNc = collect(simplify(subs(dxdrNabc,{da,db},{0,0})),dc);                % termos de 'dxdrN' dependentes de dc

dxdrEa = collect(simplify(subs(dxdrEabc,{db,dc},{0,0})),da);                % termos de 'dxdrE' dependentes de da
dxdrEb = collect(simplify(subs(dxdrEabc,{da,dc},{0,0})),db);                % termos de 'dxdrE' dependentes de db
dxdrEc = collect(simplify(subs(dxdrEabc,{da,db},{0,0})),dc);                % termos de 'dxdrE' dependentes de dc

dxdrDa = collect(simplify(subs(dxdrDabc,{db,dc},{0,0})),da);                % termos de 'dxdrD' dependentes de da
dxdrDb = collect(simplify(subs(dxdrDabc,{da,dc},{0,0})),db);                % termos de 'dxdrD' dependentes de db
dxdrDc = collect(simplify(subs(dxdrDabc,{da,db},{0,0})),dc);                % termos de 'dxdrD' dependentes de dc

dxdrNEDa = collect(simplify(dxdrNa + dxdrEa + dxdrDa),da);                  % termos de 'dxdrN + dxdrE + dxdrD' dependentes de da
dxdrNEDb = collect(simplify(dxdrNb + dxdrEb + dxdrDb),db);                  % termos de 'dxdrN + dxdrE + dxdrD' dependentes de db
dxdrNEDc = collect(simplify(dxdrNc + dxdrEc + dxdrDc),dc);                  % termos de 'dxdrN + dxdrE + dxdrD' dependentes de dc

pretty(simplify(dxdrNEDa))
pretty(simplify(dxdrNEDb))
pretty(simplify(dxdrNEDc))

dxdrNEDa = subs(dxdrNEDa,rE^2 + rN^2,a^2 + b^2 + c^2 - rD^2);
dxdrNEDb = subs(dxdrNEDb,rE^2 + rN^2,a^2 + b^2 + c^2 - rD^2);
dxdrNEDc = subs(dxdrNEDc,rE^2 + rN^2,a^2 + b^2 + c^2 - rD^2);

dxdrNEDa = subs(dxdrNEDa,rD,-(a*ax + b*ay + c*az)/gp);
dxdrNEDb = subs(dxdrNEDb,rD,-(a*ax + b*ay + c*az)/gp);
dxdrNEDc = subs(dxdrNEDc,rD,-(a*ax + b*ay + c*az)/gp);

dxdrNEDa = subs(dxdrNEDa,gp,sqrt(ax^2 + ay^2 + az^2));
dxdrNEDb = subs(dxdrNEDb,gp,sqrt(ax^2 + ay^2 + az^2));
dxdrNEDc = subs(dxdrNEDc,gp,sqrt(ax^2 + ay^2 + az^2));

pretty(simplify(dxdrNEDa))
pretty(simplify(dxdrNEDb))
pretty(simplify(dxdrNEDc))

% b = c = 0

dxdrNEDa_bc0 = subs(dxdrNEDa,{b,c},{0,0});
dxdrNEDb_bc0 = subs(dxdrNEDb,{b,c},{0,0});
dxdrNEDc_bc0 = subs(dxdrNEDc,{b,c},{0,0});

pretty(simplify(dxdrNEDa_bc0))
pretty(simplify(dxdrNEDb_bc0))
pretty(simplify(dxdrNEDc_bc0))

% a = c = 0

dxdrNEDa_ac0 = subs(dxdrNEDa,{a,c},{0,0});
dxdrNEDb_ac0 = subs(dxdrNEDb,{a,c},{0,0});
dxdrNEDc_ac0 = subs(dxdrNEDc,{a,c},{0,0});

pretty(simplify(dxdrNEDa_ac0))
pretty(simplify(dxdrNEDb_ac0))
pretty(simplify(dxdrNEDc_ac0))

% a = b = 0

dxdrNEDa_ab0 = subs(dxdrNEDa,{a,b},{0,0});
dxdrNEDb_ab0 = subs(dxdrNEDb,{a,b},{0,0});
dxdrNEDc_ab0 = subs(dxdrNEDc,{a,b},{0,0});

pretty(simplify(dxdrNEDa_ab0))
pretty(simplify(dxdrNEDb_ab0))
pretty(simplify(dxdrNEDc_ab0))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Bias Estimation Algorithms %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms etan etae etad on oe od varphin varphie varphid real

eqn1 = simplify(eta_n) == etan;
eqn2 = simplify(eta_e) == etae;
eqn3 = simplify(eta_d) == etad;
eqn4 = simplify(o_n) == on;
eqn5 = simplify(o_e) == oe;
eqn6 = simplify(o_d) == od;
eqn7 = simplify(varphi_n) == varphin;
eqn8 = simplify(varphi_e) == varphie;
eqn9 = simplify(varphi_d) == varphid;

[H,Y] = equationsToMatrix([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6],[dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp]);
rank(H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% daz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% daz_eqn3 = solve(eqn3,daz);
% daz_eqn3_est = subs(daz_eqn3,{dgp},{0});
% Ddaz_eqn3 = daz_eqn3_est - daz_eqn3;
% 
% pretty(simplify(daz_eqn3))
% pretty(simplify(daz_eqn3_est))
% pretty(simplify(Ddaz_eqn3))
%
% daz_eqn6 = solve(eqn6,daz);
% daz_eqn6_est = subs(daz_eqn6,{dgp},{0});
% Ddaz_eqn6 = daz_eqn6_est - daz_eqn6;
% 
% pretty(simplify(daz_eqn6))
% pretty(simplify(daz_eqn6_est))
% pretty(simplify(Ddaz_eqn6))

% eqn10 = subs(eqn3,dgp,0);
% eqn11 = subs(eqn6,dgp,0);
% 
% [H1,Y1] = equationsToMatrix([eqn10,eqn11],daz);
% daz_ls = inv(transpose(H1)*H1)*transpose(H1)*Y1;
% 
% pretty(simplify(daz_ls))

%%%%%%%%%%%%%%%%%%%%% dDDADR_12_ik and dDDADR_12_il %%%%%%%%%%%%%%%%%%%%%%%

clc

% dDDADR_12_ik_eqn1 = solve(eqn1,dDDADR_12_ik);
% dDDADR_12_ik_eqn1_3 = simplify(subs(dDDADR_12_ik_eqn1,daz,daz_eqn3));
% dDDADR_12_ik_eqn1_3_est = subs(dDDADR_12_ik_eqn1_3,{dDDADR_12_ij,dDDADR_12_il,dax,day,da,db},{0,0,0,0,0,0});
% DdDDADR_12_ik_eqn1_3 = dDDADR_12_ik_eqn1_3_est - dDDADR_12_ik_eqn1_3;

% pretty(simplify(dDDADR_12_ik_eqn1_3))
% pretty(simplify(dDDADR_12_ik_eqn1_3_est))
% pretty(simplify(DdDDADR_12_ik_eqn1_3))

% eqn12 = subs(eqn1,{dDDADR_12_ij,dax,day,da,db},{0,0,0,0,0});
% eqn13 = subs(eqn2,{dDDADR_12_ij,dax,day,da,db},{0,0,0,0,0});
% eqn14 = subs(eqn4,{dDDADR_12_ij,dax,day,dc},{0,0,0,0});
% eqn15 = subs(eqn5,{dDDADR_12_ij,dax,day,dc},{0,0,0,0});
% 
% eqn16 = eqn12 - eqn3*b^2/(a^2 + b^2);
% eqn17 = eqn12 + eqn6*b/a;
% eqn18 = eqn13 - eqn3*a^2/(a^2 + b^2);
% eqn19 = eqn13 + eqn6*a/b;
% eqn20 = eqn14 + eqn3*b*c/(2*(a^2 + b^2));
% eqn21 = eqn14 - eqn6*c/(2*a);
% eqn22 = eqn15 + eqn3*a*c/(2*(a^2 + b^2));
% eqn23 = eqn15 - eqn6*c/(2*b);
% 
% [H2,Y2] = equationsToMatrix([eqn16,eqn17,eqn18,eqn19,eqn20,eqn21,eqn22,eqn23],[dDDADR_12_ik,dDDADR_12_il]);
% H2 = simplify(H2);
% Y2 = simplify(Y2);
% 
% X2_ls = inv(transpose(H2)*H2)*transpose(H2)*Y2;
% X2_ls = simplify(X2_ls);
% dDDADR_12_ik_ls = X2_ls(1);
% dDDADR_12_il_ls = X2_ls(2);
% 
% pretty(simplify(dDDADR_12_ik_ls))
% pretty(simplify(dDDADR_12_il_ls))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ambiguities resolved and a = 0 %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc
% 
% eqn24 = subs(eqn1,{a,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,day,dgp},{0,0,0,0,0,0});
% eqn25 = subs(eqn2,{a,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,day,dgp},{0,0,0,0,0,0});
% eqn26 = subs(eqn3,{a,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,day,dgp},{0,0,0,0,0,0});
% eqn27 = subs(eqn4,{a,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,day,dgp},{0,0,0,0,0,0});
% eqn28 = subs(eqn5,{a,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,day,dgp},{0,0,0,0,0,0});
% eqn29 = subs(eqn6,{a,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,day,dgp},{0,0,0,0,0,0});

%%%%%%%%%%%%%%%%%%%%% LS solution to daz, db, and dc %%%%%%%%%%%%%%%%%%%%%%

% eqn30 = eqn24 - eqn26;
% eqn31 = eqn27 + eqn26*c/(2*b);
% 
% [H3,Y3] = equationsToMatrix([eqn30,eqn25,eqn26,eqn31],[daz,db,dc]);
% H3 = simplify(H3);
% Y3 = simplify(Y3);
% 
% X3_ls = inv(transpose(H3)*H3)*transpose(H3)*Y3;
% X3_ls = simplify(X3_ls);
% daz_ls = X3_ls(1);
% db_ls = X3_ls(2);
% dc_ls = X3_ls(3);
% 
% pretty(simplify(daz_ls))
% pretty(simplify(db_ls))
% pretty(simplify(dc_ls))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ambiguities resolved and b = 0 %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

eqn32 = subs(eqn1,{b,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il},{0,0,0,0});
eqn33 = subs(eqn2,{b,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il},{0,0,0,0});
eqn34 = subs(eqn3,{b,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il},{0,0,0,0});
eqn35 = subs(eqn4,{b,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il},{0,0,0,0});
eqn36 = subs(eqn5,{b,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il},{0,0,0,0});
eqn37 = subs(eqn6,{b,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il},{0,0,0,0});
eqn38 = subs(eqn7,{b,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il},{0,0,0,0});
eqn39 = subs(eqn8,{b,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il},{0,0,0,0});
eqn40 = subs(eqn9,{b,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il},{0,0,0,0});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% daz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

daz_eqn34 = solve(eqn34,daz);
daz_eqn34_est = subs(daz_eqn34,{dgp},{0});
Ddaz_eqn34 = daz_eqn34_est - daz_eqn34;

pretty(simplify(daz_eqn34))
pretty(simplify(daz_eqn34_est))
pretty(simplify(Ddaz_eqn34))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% da %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

da_eqn32 = solve(eqn32,da);
da_eqn32_est = subs(da_eqn32,{dax},{0});
Dda_eqn32 = da_eqn32_est - da_eqn32;

pretty(simplify(da_eqn32))
pretty(simplify(da_eqn32_est))
pretty(simplify(Dda_eqn32))

eqn41 = eqn33 - eqn34;

da_eqn41 = solve(eqn41,da);
da_eqn41_est = subs(da_eqn41,{dax},{0});
Dda_eqn41 = da_eqn41_est - da_eqn41;

pretty(simplify(da_eqn41))
pretty(simplify(da_eqn41_est))
pretty(simplify(Dda_eqn41))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eqn42 = eqn36 + eqn34*c/(2*a);

dc_eqn42 = solve(eqn42,dc);
dc_eqn42_est = subs(dc_eqn42,{dax},{0});
Ddc_eqn42 = dc_eqn42_est - dc_eqn42;

pretty(simplify(dc_eqn42))
pretty(simplify(dc_eqn42_est))
pretty(simplify(Ddc_eqn42))

%%%%%%%%%%%%%%%%%%%%% LS solution to daz, da, and dc %%%%%%%%%%%%%%%%%%%%%%

eqn43 = subs(eqn32,{dax,dgp},{0,0});
eqn44 = subs(eqn34,{dax,dgp},{0,0});
eqn45 = subs(eqn41,{dax,dgp},{0,0});
eqn46 = subs(eqn42,{dax,dgp},{0,0});

[H4,Y4] = equationsToMatrix([eqn43,eqn45,eqn44,eqn46],[daz,da,dc]);
H4 = simplify(H4);
Y4 = simplify(Y4);

X4_ls = inv(transpose(H4)*H4)*transpose(H4)*Y4;
X4_ls = simplify(X4_ls);
daz_ls = X4_ls(1);
da_ls = X4_ls(2);
dc_ls = X4_ls(3);

pretty(simplify(daz_ls))
pretty(simplify(da_ls))
pretty(simplify(dc_ls))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% varphi's %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varphi_n_eqn38 = solve(eqn38,varphin);
varphi_d_eqn40 = solve(eqn40,varphid);

eqn47 = eqn39 + eqn34*c/(2*a);
eqn48 = subs(eqn47,dc,dc_eqn42);

varphi_e_eqn48 = solve(eqn48,varphie);

pretty(simplify(varphi_n_eqn38))
pretty(simplify(varphi_e_eqn48))
pretty(simplify(varphi_d_eqn40))

clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% teste numérico %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

gSTD = 9.80665;                                                             % standard gravity - para conversao entre m/s² e g

gP = 9.786466787716226;                                                     % gravidade LINCS
A = 2;
B = 0;
C = 0;                                                                      % baseline coordinates in body-frame

bDDADR_12_ij = 0;                                                           % double-differenced ADR bias [m]
bDDADR_12_ik = 0;                                                           % double-differenced ADR bias [m]
bDDADR_12_il = 0;                                                           % double-differenced ADR bias [m]
bax = 1*gSTD/1000;                                                          % bias dos acelerômetros (mg)
bay = -1*gSTD/1000;                                                         % bias dos acelerômetros (mg)
baz = 1*gSTD/1000;                                                          % bias dos acelerômetros (mg)
ba = -0.02;                                                                 % bias baseline no eixo x (m)
bb = 0.02;                                                                  % bias baseline no eixo y (m)
bc = 0.02;                                                                 % bias baseline no eixo z (m)
bgrav = -0.005*gSTD/1000;                                                   % bias na gravidade (mg) - Savage

eta_n_teste = double(subs(eta_n,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));
eta_e_teste = double(subs(eta_e,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));
eta_d_teste = double(subs(eta_d,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));
o_n_teste = double(subs(o_n,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));
o_e_teste = double(subs(o_e,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));
o_d_teste = double(subs(o_d,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));
varphi_n_teste = double(subs(varphi_n,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));
varphi_e_teste = double(subs(varphi_e,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));
varphi_d_teste = double(subs(varphi_d,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));

Ddaz_teste = double(subs(Ddaz_eqn34,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));
Dda_teste = double(subs(Dda_eqn32,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));
Ddc_teste = double(subs(Ddc_eqn42,{a,b,c,gp,dDDADR_12_ij,dDDADR_12_ik,dDDADR_12_il,dax,day,daz,da,db,dc,dgp},{A,B,C,gP,bDDADR_12_ij,bDDADR_12_ik,bDDADR_12_il,bax,bay,baz,ba,bb,bc,bgrav}));

format short

eta_n_teste*180/pi
eta_e_teste*180/pi
eta_d_teste*180/pi
o_n_teste*180/pi
o_e_teste*180/pi
o_d_teste*180/pi
varphi_n_teste*180/pi
varphi_e_teste*180/pi
varphi_d_teste*180/pi

Ddaz_teste*1000/gSTD
Dda_teste
Ddc_teste

