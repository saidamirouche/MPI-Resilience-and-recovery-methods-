
# LOAD different data files 

data_nf = load('nf_time.txt');
data_er = load('er.txt');
data_li = load('li.txt');
data_reset = load('reset.txt');

data_er_1 = load('er_1.txt');
data_li_1 = load('li_1.txt');
data_reset_1 = load('reset_1.txt');

data_er_2 = load('er_2.txt');
data_li_2 = load('li_2.txt');
data_reset_2 = load('reset_2.txt');

data_er_3 = load('er_3.txt');
data_li_3 = load('li_3.txt');
data_reset_3 = load('reset_3.txt');



# PLOTING the iterations data files 
subplot(221);

plot(data_nf(:,2),(data_nf(:,1)),'-k','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
#my_colour = [255 212 59] ./ 255;
my_colour = [251 167 18] ./ 255;
plot(data_er_1(:,2),(data_er_1(:,1)),'-','LineWidth',2,'Color',my_colour);
set(gca, 'YScale', 'log')
hold on
plot(data_li_1(:,2),(data_li_1(:,1)),'-g','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
plot(data_reset_1(:,2),(data_reset_1(:,1)),'-r','LineWidth',2);
set(gca, 'YScale', 'log')
legend("NF","ER","LI","RESET");
xlabel({"Iterations"},'FontSize',16,'FontWeight','bold','Color','b')
ylabel({" ||(b-Ax)||/||b||"},'FontSize',16,'FontWeight','bold','Color','b')
title("4 processes killed between [100, 160] iterations",'FontSize',16,'FontWeight','bold','Color','k')

subplot(222);



plot(data_nf(:,2),(data_nf(:,1)),'-k','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
#my_colour = [255 212 59] ./ 255;
my_colour = [251 167 18] ./ 255;
plot(data_er_2(:,2),(data_er_2(:,1)),'-','LineWidth',2,'Color',my_colour);
set(gca, 'YScale', 'log')
hold on
plot(data_li_2(:,2),(data_li_2(:,1)),'-g','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
plot(data_reset_2(:,2),(data_reset_2(:,1)),'-r','LineWidth',2);
set(gca, 'YScale', 'log')
legend("NF","ER","LI","RESET");
xlabel({"Iterations"},'FontSize',16,'FontWeight','bold','Color','b')
ylabel({" ||(b-Ax)||/||b||"},'FontSize',16,'FontWeight','bold','Color','b')
title("4 processes killed between [300, 360] iterations",'FontSize',16,'FontWeight','bold','Color','k')

subplot(223);


plot(data_nf(:,2),(data_nf(:,1)),'-k','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
#my_colour = [255 212 59] ./ 255;
my_colour = [251 167 18] ./ 255;
plot(data_er_3(:,2),(data_er_3(:,1)),'-','LineWidth',2,'Color',my_colour);
set(gca, 'YScale', 'log')
hold on
plot(data_li_3(:,2),(data_li_3(:,1)),'-g','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
plot(data_reset_3(:,2),(data_reset_3(:,1)),'-r','LineWidth',2);
set(gca, 'YScale', 'log')
legend("NF","ER","LI","RESET");
xlabel({"Iterations"},'FontSize',16,'FontWeight','bold','Color','b')
ylabel({" ||(b-Ax)||/||b||"},'FontSize',16,'FontWeight','bold','Color','b')
title("4 processes killed between [600, 660] iterations",'FontSize',16,'FontWeight','bold','Color','k')

subplot(224);

plot(data_nf(:,2),(data_nf(:,1)),'-k','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
#my_colour = [255 212 59] ./ 255;
my_colour = [251 167 18] ./ 255;
plot(data_er(:,2),(data_er(:,1)),'-','LineWidth',2,'Color',my_colour);

hold on
plot(data_li(:,2),(data_li(:,1)),'-g','LineWidth',2);

hold on 
plot(data_reset(:,2),(data_reset(:,1)),'-r','LineWidth',2);

legend("NF","ER","LI","RESET");
xlabel({"Iterations"},'FontSize',16,'FontWeight','bold','Color','b')
ylabel({" ||(b-Ax)||/||b||"},'FontSize',16,'FontWeight','bold','Color','b')
title("4 processes killed randomly in any moment",'FontSize',16,'FontWeight','bold','Color','k')


#Saving the plot in a photo
set(gcf,'PaperUnits','inches','PaperPosition',[1 1 19 12])

print -djpeg nbr_iter.jpg -r100

figure 


# PLOTING the excution time data files 


subplot(221);


plot((data_nf(:,3)),(data_nf(:,1)),'-k','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
#my_colour = [255 212 59] ./ 255;
my_colour = [251 167 18] ./ 255;
plot((data_er_1(:,3)),(data_er_1(:,1)),'-','LineWidth',2,'Color',my_colour);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

hold on
plot((data_li_1(:,3)),(data_li_1(:,1)),'-g','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
plot((data_reset_1(:,3)),(data_reset_1(:,1)),'-r','LineWidth',2);
set(gca, 'YScale', 'log')
legend("NF","ER","LI","RESET",'location','southwest' );

xlabel({"Time"},'FontSize',16,'FontWeight','bold','Color','b')
ylabel({" ||(b-Ax)||/||b||"},'FontSize',16,'FontWeight','bold','Color','b')
title("4 processes killed between [100, 160] iterations",'FontSize',16,'FontWeight','bold','Color','k')

subplot(222);
plot((data_nf(:,3)),(data_nf(:,1)),'-k','LineWidth',2);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

hold on 
#my_colour = [255 212 59] ./ 255;
my_colour = [251 167 18] ./ 255;
plot((data_er_2(:,3)),(data_er_2(:,1)),'-','LineWidth',2,'Color',my_colour);
set(gca, 'YScale', 'log')
hold on
plot((data_li_2(:,3)),(data_li_2(:,1)),'-g','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
plot((data_reset_2(:,3)),(data_reset_2(:,1)),'-r','LineWidth',2);
set(gca, 'YScale', 'log')
legend("NF","ER","LI","RESET",'location','southwest' );
xlabel({"Time"},'FontSize',16,'FontWeight','bold','Color','b')
ylabel({" ||(b-Ax)||/||b||"},'FontSize',16,'FontWeight','bold','Color','b')
title("4 processes killed between [300, 360] iterations",'FontSize',16,'FontWeight','bold','Color','k')

subplot(223);
plot((data_nf(:,3)),(data_nf(:,1)),'-k','LineWidth',2);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

hold on 
#my_colour = [255 212 59] ./ 255;
my_colour = [251 167 18] ./ 255;
plot((data_er_3(:,3)),(data_er_3(:,1)),'-','LineWidth',2,'Color',my_colour);
set(gca, 'YScale', 'log')
hold on
plot((data_li_3(:,3)),(data_li_3(:,1)),'-g','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
plot((data_reset_3(:,3)),(data_reset_3(:,1)),'-r','LineWidth',2);
set(gca, 'YScale', 'log')
legend("NF","ER","LI","RESET",'location','southwest' );
xlabel({"Time"},'FontSize',16,'FontWeight','bold','Color','b')
ylabel({" ||(b-Ax)||/||b||"},'FontSize',16,'FontWeight','bold','Color','b')
title("4 processes killed between [600, 660] iterations",'FontSize',16,'FontWeight','bold','Color','k')

subplot(224);

plot((data_nf(:,3)),(data_nf(:,1)),'-k','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
#my_colour = [255 212 59] ./ 255;
my_colour = [251 167 18] ./ 255;
plot((data_er(:,3)),(data_er(:,1)),'-','LineWidth',2,'Color',my_colour);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

hold on
plot((data_li(:,3)),(data_li(:,1)),'-g','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
plot((data_reset(:,3)),(data_reset(:,1)),'-r','LineWidth',2);
set(gca, 'YScale', 'log')
legend("NF","ER","LI","RESET",'location','southwest' );
xlabel({"Time"},'FontSize',16,'FontWeight','bold','Color','b')
ylabel({" ||(b-Ax)||/||b||"},'FontSize',16,'FontWeight','bold','Color','b')
title("4 processes killed randomly in any moment",'FontSize',16,'FontWeight','bold','Color','k')

set(gcf,'PaperUnits','inches','PaperPosition',[1 1 19 12])

print -djpeg timees.jpg -r100


