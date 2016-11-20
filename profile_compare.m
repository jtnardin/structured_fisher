%profile_compare.m written 9-21-16 by JTN

height = [1 2 4 8];

decay = [.1 .25 .5 1];

y = cell(4,4);
z = cell(4,4);

for i = 1:4
    for j = 1:4
        [i,j]
        [soln_mx, soln_x] = fisher_diffusion_struct_Dm_f([height(i),decay(j)]);
        y{i,j} = soln_mx;
        z{i,j} = soln_x;
    end
end

xn = 151; %number of x points
x = linspace(0,30,xn);

figure
count = 1;

for i = 1:4
    for j = 1:4
        
        subplot(4,4,count)
        
        plot(x,z{i,j});
        count = count + 1;
        title(['height = ' num2str(height(i)) ' decay = ' num2str(decay(j))])
    end
end

save('ex2_take1.mat')