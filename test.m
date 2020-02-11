f = 1:.2:100;
SF1 = 1./f+2;
SF2 = 1./f+3;

figure;
plot(f,SF1)
hold on;
plot(f,SF2)

figure;
plot(f,2./SF1); hold on
plot(f,2./SF2);

plot(f,(2./SF1)-(2./SF2))