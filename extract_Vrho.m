clear
load('staged377_invc_adj.mat','muidx','rot','dMx','muidx')
numer=rot.numer;
denom=rot.denom;
labA=10;
chart_p = coco_read_solution('stagedbc377', labA,'chart');
x=chart_p.x(muidx(1:2,:));
[r2,b]=deal(chart_p.x(5),chart_p.x(6));
perm=diag(ones(numer,1),denom-numer)+diag(ones(denom-numer,1),-numer);
sigma=kron(perm,eye(2));
xrot=num2cell(reshape(sigma*x(:),2,denom),1);
dMxval=cellfun(@(x)dMx(x,r2,b),xrot,'uniformoutput',false);
V_rho=blkdiag(dMxval{:})*sigma;
evs=eig(V_rho);
plot(real(evs),imag(evs),'o')
axis equal