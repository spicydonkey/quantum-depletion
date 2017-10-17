% TF fit on existing data

% config fit
fit_tf.eqn=@(p,x1) max(p(1)*(1-(x1/p(2)).^2),0);
fit_tf.coefname={'n0','R'};
fit_tf.param0=[3e4,2];
fit_tf.opt=statset('TolFun',1e-10,'TolX',1e-10,'MaxIter',1e3);

%
xx=k*1e-6;
yy=nk*1e18;

fit_tf.fit=fitnlm(xx,yy,fit_tf.eqn,fit_tf.param0,...
    'CoefficientNames',fit_tf.coefname,...
    'Options',fit_tf.opt);

xf=linspace(min(xx),max(xx),1e3);
yf=feval(fit_tf.fit,xf);