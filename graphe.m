data_nf = load('nf.txt');
data_er = load('er.txt');
data_li = load('li.txt');
data_reset = load('reset.txt');

plot(data_nf(:,2),log10(data_nf(:,1)),'-k','LineWidth',3);
hold on 
my_colour = [255 212 59] ./ 255;
plot(data_er(:,2),log10(data_er(:,1)),'-','LineWidth',3,'Color',my_colour);
hold on
plot(data_li(:,2),log10(data_li(:,1)),'-g','LineWidth',3);
hold on 
plot(data_reset(:,2),log10(data_reset(:,1)),'-r','LineWidth',3);
legend("NF","ER","LI","RESET");
xlabel("Number of iterations")
ylabel("convergence")
title("Numerical behavior of recovery strategies with ULFM")

print -djpg image2.jpg