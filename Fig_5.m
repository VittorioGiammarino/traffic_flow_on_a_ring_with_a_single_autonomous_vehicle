load("Fig_5_stability_in_the_a-b_plane.mat")

figure()
imagesc(b_bar,a_bar*h_eq*h_eq,Stability_matrix)
%hold on
%plot(0.1,140,'ro','MarkerSize',20)
hold on
plot(0.5,20,'rx','MarkerSize',20)
hold on
plot(b_bar_bound,a_bar_bound,'w-o','MarkerSize',1.2)
caxis([1 9])
colorbar('Ticks',[1,2,3,4,5,6,7,8,9],...
         'TickLabels',{'3-vehicles','5-vehicles','10-vehicles',...
         '20-vehicles','60-vehicles','100-vehicles','200-vehicles','500-vehicles',...
         'Stability'})
colormap('jet')
legend({'a=20, b=0.5', 'Bound for ||\Gamma||_{\infty}\leq 1'},'TextColor','w','FontSize',19)
legend('boxoff')
xlabel('b')
ylabel('a')
set(gcf, 'Position',  [100, 100, 800, 400])