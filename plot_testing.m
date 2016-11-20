
     fig_final = figure('units','normalized','outerposition',[0 0 1 1]);

               
        count = 1;
        step = floor(tn/3);
        for i = step:step:tn
        
            
            subplot(18,10,60*(count-1)+[2:5 12:15 22:25])
                        
                contourf(x,m,y(:,:,i),'edgecolor','none')
                hold on
                title(['Example 3, u(t = '  num2str(round(t(i))) ',x,m)'])
                axis([0 25 0 1 0 umax])
                caxis([0,5])
                plot3([0 30],[.5 .5],[5 5],'color',[1 1 1])
                %colorbar
                view(2)
                set(gca,'ytick',[])
                set(gca,'xtick',[])
            
            subplot(18,10,60*(count-1)+[1 11 21])
          
                plot([IC_1_d_m(1) Soln(t(i),m_fine(2:end-1)) IC_1_d_m(end)]/max(Soln(t(i),m_fine(2:end-1))),m_fine,'linewidth',3)
                axis([0 1.1 0 1])

                set(gca,'xdir','reverse')
                ylabel('m')
                xlabel('u(t,m)')
                
                %label for plots
                text(1.5,1.1,['(' char(96+count) ')'],'fontsize',15)
            
          
            subplot(18,10,60*(count-1)+[32:35 42:45])
            
                plot(x,z_tx(:,i),'linewidth',3)
                hold on
                plot(x,z_nonaut(i,:),'color',[0 .5 0],'linewidth',3)
                axis([0 13 0 1.11])
                
                if count == 3
                    xlabel('x')
                end
                
                ylabel('w(t,x)')
%                 set(gca)

                if count == 1
                    h=legend('Structured simulation','Nonautonomous simulation','location','northeast');
                end
                
                          
                
            
            count = count + 1;
        end
        
        subplot(18,10,[7:10,17:20,27:30,37:40,47:50,57:60,67:70,77:80])
            hold on

            tcont = linspace(0,15,300);
            scont = linspace(0,1,400);

            [Tcont,Scont] = meshgrid(tcont,scont);

            A = Soln(Tcont,Scont);

            contourf(Tcont,Scont,log(A),'edgecolor','none')        


            plot(tcont,sigma_inv(int_f_s(tcont),0.35),'color',[0.5 0.5 0.5],'linewidth',3)

            % h=legend('Distribution of m','$h(t;\underline{m})$');
            % set(h,'interpreter','latex')

            plot(tcont,sigma_inv(int_f_s(tcont),0.15),'color',[0.5 0.5 0.5],'linewidth',3)
            plot(tcont,sigma_inv(int_f_s(tcont),0.05),'color',[0.5 0.5 0.5],'linewidth',3)


            title('Distribution along $m$ over time in Example 3','interpreter','latex')
            xlabel('t')
            ylabel('m')

            caxis([min(min(log(A(A~=0)))) max(max(log(A)))]/1.5)
            colorbar
            text(-1.15,1.05,'(d)','fontsize',15)
            
        
        
        
        s3 = subplot(18,10,[107:110 117:120 127:130 137:140 147:150 157:160 167:170 177:180]);
            
        h1 = openfig('Ex3_range_gamma.fig','reuse');
        ax1 = gca;
        fig3 = get(ax1,'children');
        
        copyobj(fig3,s3);
        
        figure(fig_final)
        subplot(18,10,[107:110 117:120 127:130 137:140 147:150 157:160 167:170 177:180]);
        
        xlabel('x')
        ylabel('w(t=15,x)')
        title('Various values of $\gamma$','interpreter','latex')
        axis([0 25 0 1.2])
        text(-1.25,1.05,'(e)','fontsize',15)
%         legend('2.4','2.55','3.594','4.537','5.68','6','6.2','location','northeast')
        
        
    set(fig_final,'color',[1 1 1])
    export_fig(fig_final,'ex_3_full.eps')
    saveas(fig_final,'ex_3_full.fig')
        