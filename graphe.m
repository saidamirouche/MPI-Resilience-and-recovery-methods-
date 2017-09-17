
# LOAD different data files 


data_er_1 = load('er_1_no_agree.txt');
data_li_1 = load('li_1_no_agree.txt');
data_reset_1 = load('reset_1_no_agree.txt');
data_er_1_agree = load('er_1_agree.txt');
data_li_1__agree = load('li_1_agree.txt');
data_reset_1__agree = load('reset_1_agree.txt');


data_er_2 = load('er_2_no_agree.txt');
data_li_2 = load('li_2_no_agree.txt');
data_reset_2 = load('reset_2_no_agree.txt');
data_er_2_agree  = load('er_2_agree.txt');
data_li_2_agree  = load('li_2_agree.txt');
data_reset_2_agree  = load('reset_2_agree.txt');


data_er_3 = load('er_3_no_agree.txt');
data_li_3 = load('li_3_no_agree.txt');
data_reset_3 = load('reset_3_no_agree.txt');
data_er_3_agree = load('er_3_agree.txt');
data_li_3_agree = load('li_3_agree.txt');
data_reset_3_agree = load('reset_3_agree.txt');



data_nf = load('nf_rand.txt');
data_er = load('er_rand_no_agree.txt');
data_li = load('LI_rand_no_agree.txt');
data_reset = load('reset_rand_no_agree.txt');
data_er_agree = load('er_rand_agree.txt');
data_li_agree = load('li_rand_agree.txt');
data_reset_agree = load('reset_rand_agree.txt');

no_agree = load('nf_non_agree.txt');
relax = load('er_sync_3.txt');



# PLOTING THE DIFFERENCE BETWEEN USING AGREEMENT AND NOT USING IT
plot(data_nf(:,3),(data_nf(:,1)),'-k','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
plot(no_agree(:,3),(no_agree(:,1)),'--r','LineWidth',2);
ax = legend("With Agreement","Without Agreement");
set(ax, 'fontsize', 15,'FontWeight','bold');
title("Time execution between using Agreement and not using it",'fontsize', 12)
set(gcf,'PaperUnits','inches','PaperPosition',[1 1 19 12])

print -djpeg agree.jpg -r100


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
title("4 processes killed between [200, 260] iterations",'FontSize',16,'FontWeight','bold','Color','k')

subplot(223);


plot(data_nf(:,2),(data_nf(:,1)),'-k','LineWidth',2);
set(gca, 'YScale', 'log')
hold on 
#my_colour = [255 212 59] ./ 255;
my_colour = [251 167 18] ./ 255;
plot(data_er_2(:,2),(data_er_2(:,1)),'-','LineWidth',2,'Color',my_colour);
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
title("4 processes killed between [300, 360] iterations",'FontSize',16,'FontWeight','bold','Color','k')

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

hold on
plot(relax(:,2),(relax(:,1)),'--r','LineWidth',2);

ax = legend("NF","ER","LI","RESET","Relaxed");
#set(ax, 'fontsize', 16,'FontWeight','bold');
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
plot((data_er_2(:,3)),(data_er_2(:,1)),'-','LineWidth',2,'Color',my_colour);
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

plot((data_nf(:,3)),(data_nf(:,1)),'-k','LineWidth',3);
set(gca, 'YScale', 'log')
hold on 
#my_colour = [255 212 59] ./ 255;
my_colour = [251 167 18] ./ 255;
plot((data_er(:,3)),(data_er(:,1)),'-','LineWidth',3,'Color',my_colour);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

hold on
plot((data_li(:,3)),(data_li(:,1)),'-g','LineWidth',3);
set(gca, 'YScale', 'log')
hold on 
plot((data_reset(:,3)),(data_reset(:,1)),'-r','LineWidth',3);
set(gca, 'YScale', 'log')
hold on
plot(no_agree(:,3),(no_agree(:,1)),'--k','LineWidth',3);
set(gca, 'YScale', 'log')

hold on
plot((data_er_agree(:,3)),(data_er_agree(:,1)),'--','LineWidth',3,'Color',my_colour);

hold on
plot((data_li_agree(:,3)),(data_li_agree(:,1)),'--g','LineWidth',3);

hold on
plot((data_reset_agree(:,3)),(data_reset_agree(:,1)),'--r','LineWidth',3);


hold on 
plot(relax(:,3),(relax(:,1)),'-b','LineWidth',3);
set(gca, 'YScale', 'log')

 axis tight;
 
ax = legend("NF ","ER Without Agree","LI Without Agree","RESET Without Agree","NF without Agree","ER","LI","Reset","Relaxed ",'location','southwest' );
set(ax, 'fontsize', 16,'FontWeight','bold');
xlabel({"Time"},'FontSize',16,'FontWeight','bold','Color','b')
ylabel({" ||(b-Ax)||/||b||"},'FontSize',16,'FontWeight','bold','Color','b')
title("4 processes killed randomly in any moment",'FontSize',16,'FontWeight','bold','Color','k')

set(gcf,'PaperUnits','inches','PaperPosition',[1 1 19 12])

print -djpeg timee_all.jpg -r100



