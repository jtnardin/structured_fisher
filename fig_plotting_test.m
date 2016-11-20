
        figure('units','normalized','outerposition',[0 0 1 1])

               
        count = 1;
        step = floor(tn/3);
        for i = step:step:tn
        
            
            subplot(18,10,60*(count-1)+[2:5 12:15 22:25])
                        
                contourf(x,m,y(:,:,i),'edgecolor','none')
                hold on
                title(['Example 1, u(t = '  num2str(round(t(i))) ',x,m)'])
    %             xlabel('x','fontsize',30)
                axis([0 13 0 1 0 umax])
                caxis([0,5])
                plot3([0 30],[.5 .5],[5 5],'color',[1 1 1])
                %colorbar
                view(2)
                set(gca,'ytick',[])
                set(gca,'xtick',[])
            
            subplot(18,10,60*(count-1)+[1 11 21])
          
                plot([IC_1_d_m(1) Soln(t(i),m_fine(2:end-1)) IC_1_d_m(end)]/max(Soln(t(i),m_fine(2:end-1))),m_fine,'linewidth',6)
                axis([0 1.1 0 1])

                set(gca,'xdir','reverse')
                ylabel('m')
                xlabel('u(t,m)')
                
                %label for plots
                text(1.5,1.1,['(' char(96+count) ')'])
            
          
            subplot(18,10,60*(count-1)+[32:35 42:45])
            
                plot(x,z(:,i),'linewidth',6)
                hold on
                plot(x,z_nonaut(i,:),'color',[0 .5 0],'linewidth',6)
                axis([0 13 0 1.11])
                xlabel('x')
                ylabel('w(t,x)')
                set(gca)

                if count == 1
                    h=legend('Structured simulation','Nonautonomous simulation','location','northeast');
                end
                
                          
                
            
            count = count + 1;
        end
        
        subplot(18,10,[6:10,16:20,26:30,36:40,46:50,56:60, 66:70])
            hold on

            tcont = linspace(0,10,300);
            scont = linspace(0,1,400);

            [Tcont,Scont] = meshgrid(tcont,scont);

            A = Soln(Tcont,Scont);

            contourf(Tcont,Scont,log(A),'edgecolor','none')        


            plot(tcont,sigma_inv(int_f_s(tcont),0.35),'color',[0.5 0.5 0.5],'linewidth',1.5)

            % h=legend('Distribution of m','$h(t;\underline{m})$');
            % set(h,'interpreter','latex')

            plot(tcont,sigma_inv(int_f_s(tcont),0.15),'color',[0.5 0.5 0.5],'linewidth',1.5)
            plot(tcont,sigma_inv(int_f_s(tcont),0.05),'color',[0.5 0.5 0.5],'linewidth',1.5)


            title('Distribution along $m$ over time in Example 1','interpreter','latex')
            xlabel('t')
            ylabel('m')

            caxis([min(min(log(A(A~=0)))) max(max(log(A)))]/1.5)
            colorbar
            view(2)
        
        cd janus_code/
        subplot(18,10,[86:90 96:100 106:110 116:120 126:130 136:140 146:150])
        
            load('ex1_range_alpha')


            hold on

            for i = 1:length(alpha)

                plot(x,z{i}(:,end),'linewidth',.25)

            end


            arrow([10,0.4],[30,0.4])
            text(10,0.45,'$\alpha$ increasing','interpreter','latex')

            axis([0 40 0 1])

            xlabel('x')
            ylabel('w(t,x)')

            title('Example 1, w(t=30,x) for various $\alpha$ values','interpreter','latex')

        cd ..

        
    set(gcf,'color',[1 1 1])
    export_fig(gcf,'ex_1_testing.eps')
        