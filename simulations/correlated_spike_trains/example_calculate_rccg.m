sp = synthetic_spikes(0.05,30e-3,5);
% sp = zeros(16e4,30);
% sp(randi(numel(sp),30*10,1)) = 1;

dt = 1/16e3;


count = 0;
for i = 1:5
    [ccg(:,:,i),I,J,ccg0(:,:,i)] = calculate_CCG(sp,dt,i,i);
    for j = i+1:5
        count = count+1;
        [ccg_null(:,:,count),~,~,ccg0_null(:,:,count)] = calculate_CCG(sp,dt,i,j);
    end
end

ccg = mean(ccg,3)-mean(ccg_null,3);
ccg0 = mean(ccg0,3)-mean(ccg0_null,3);


kernl = normpdf(linspace(-2,2,5e-3/dt),0,1);
kernl = kernl/sum(kernl);

% ccg_filter = filter(kernl,1,ccg);
ccg_filter = ccg;
% ccg0_filter = filter(kernl,1,ccg0);
ccg0_filter = ccg0;


m = (size(ccg,1)-1)/2;
t_lag = [10.^linspace(0,2,20),150,200,300,400];
r_ccg = zeros(length(t_lag),size(ccg_filter,2));
for i = 1:length(t_lag)
    it = floor(t_lag(i)/dt*1e-3);
    idcs = m-it:m+it;
    r_ccg(i,:) = sum(ccg_filter(idcs,:))./sqrt(sum(ccg0_filter(idcs,I(:))).*sum(ccg0_filter(idcs,J(:))));
end

figureNB;
    plot(t_lag,mean(r_ccg,2),'.-k','LineWidth',2,'MarkerSize',20)
    set(gca,'xscale','log')


function [ccg,I,J,ccg0] = calculate_CCG(sp,dt,t1,t2)
    h = waitbar(0);
    N = size(sp,2);
    count = 0;
    ccg = zeros(1/dt+1,N*(N-1)/2);
    ccg0 = zeros(1/dt+1,N);
    I = nan(N*(N-1)/2,1);
    J = nan(N*(N-1)/2,1);
    for i = 1:N
        ccg0(:,i) = CCG(sp(:,i,t1),sp(:,i,t2),0.5,dt);
    end
    for i = 1:N
        for j = i+1:N
            count = count+1;
            update_waitbar(h,count,N*(N-1)/2);
            ccg(:,count) = CCG(sp(:,i,t1),sp(:,j,t2),0.5,dt);
            I(count) = i;
            J(count) = j;
        end
    end
    delete(h);
end