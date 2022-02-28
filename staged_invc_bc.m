%% Demonstration for staged construction of invariant curve and the adjoint problem for phase response
%
% Part of H. Dankowicz, J. Sieber, "Sensitivity Analysis for Periodic
% Orbits and Quasiperiodic Invariant Tori Using the Adjoint Method"
%
%%
clear
%% approximate golden ratio rotation number by continued fraction
gr=(sqrt(sym(5))-1)/2;
[ngr,dgr]=numden(str2sym(rat(double(gr),1e-5)));
rot=struct('numer',double(ngr),'denom',double(dgr));
%% System parameters, map and its derivatives
[theta,alpha,a,r1]=deal(2*pi*rot.numer/rot.denom,1/4,-1/4,0);
mrot=[cos(theta),-sin(theta);sin(theta),cos(theta)];
M=   @(x,r2,b) (1+alpha)*mrot*x+(x.'*x)*mrot*[a,-b;b,a]*x+[r1;r2].*x;
dMx= @(x,r2,b) (1+alpha)*mrot+(x.'*x)*mrot*[a,-b;b,a]+diag([r1;r2])+2*mrot*[a,-b;b,a]*x*x.';
dMr2=@(x,r2,b) [0;x(2)];
dMb= @(x,r2,b) mrot*[0,-1;1,0]*(x.'*x)*x;
%% convert to M(x,p) format, parameters are r2 and b
Mxp=  @(x,p)M(x,p(1),p(2));
dMxpx=@(x,p)dMx(x,p(1),p(2));
dMxpp=@(x,p)[dMr2(x,p(1),p(2)),dMb(x,p(1),p(2))];
%% Coco function is residual of M
mx = 1:2; my = 3:4; mp= 5:6; mdrho=7;
Mres=@(u)Mxp(u(mx),u(mp))-u(my);
dMres=@(u)[dMxpx(u(mx),u(mp)),-eye(2),dMxpp(u(mx),u(mp)),zeros(2,1)];
bx=1:2; bxnext=3:4; by=5:6; bdrho=7;
fbc=@(u)u(bx)+u(bdrho)*rot.denom*(u(bxnext)-u(bx))-u(by);
dbc=@(u)[eye(2)*(1-rot.denom*u(bdrho)),eye(2)*rot.denom*u(bdrho),-eye(2),...
    rot.denom*(u(bxnext)-u(bx))];
fcn=@(f) @(p,d,u) deal(d, f(u)); % convert plain fcn to coco fcn
fid=@(s,i)[s,num2str(i)];        % enumerated names for function identifiers
%% initial guess  (circle of radius sqrt(-alpha/a)) and parameter values
p0=[0;0];
drho0=0;
phi=2*pi*linspace(0,1,rot.denom+1);
phi=phi(1:end-1);
x0=[cos(phi);sin(phi)]*sqrt(-alpha/a);
%% Collect rot.denom copies of the map to the continuation problem
prob = coco_prob;
for i=1:rot.denom
  irot = mod(i+rot.numer-1,rot.denom)+1;
  prob = coco_add_func(prob, fid('M',i), fcn(Mres),fcn(dMres), [], 'zero', ...
      'u0', [x0(:,i);x0(:,irot);p0;drho0]);
  prob=coco_add_adjt(prob,fid('M',i));
  muidx(:,i)=coco_get_func_data(prob,fid('M',i),'uidx');
  maidx(:,i)=coco_get_adjt_data(prob,fid('M',i),'axidx');
end
prob=coco_add_pars(prob,'pars',muidx([mp,mdrho],1),{'r2','b','drho'});
prob=coco_add_adjt(prob,'pars',{'e.r2','e.b','e.drho'},...
    'aidx',maidx([mp,mdrho],1),'l0',[0;0;1]);
%% Add boundary conditions
for i=1:rot.denom
  irot = mod(i+rot.numer-1,rot.denom)+1;
  inext = mod(i+rot.numer,rot.denom)+1;
  prob = coco_add_func(prob,fid('bc',i), fcn(fbc),fcn(dbc), [], 'zero', ...
      'uidx', [muidx(mx,irot);muidx(mx,inext);muidx(my,i);muidx(mdrho,i)]);
  prob=coco_add_adjt(prob,fid('bc',i),...
      'aidx',[maidx(mx,irot);maidx(mx,inext);maidx(my,i);maidx(mdrho,i)]);
end
%% Add glue between parameters
for i=2:rot.denom
  prob = coco_add_glue(prob,fid('pglue',i),muidx([mp,mdrho],1),muidx([mp,mdrho],i));
  prob=coco_add_adjt(prob,fid('pglue',i),'aidx',[maidx([mp,mdrho],1);maidx([mp,mdrho],i)]);
end
%% Add phase condition
phascond=@(p,d,u)deal(d,d.dx(:)'*(u(:)-d.xref(:)));
dphascond=@(p,d,u)deal(d,d.dx(:)');
data=struct('dx',rot.denom*(x0(:,[2:end,1])-x0),'xref',x0);
prob = coco_add_func(prob, 'phasecond',phascond,dphascond,data,'zero', ...
    'uidx',muidx(mx,:));
prob=coco_add_adjt(prob,'phasecond','aidx',maidx(mx,:));
%% continuation settings
prob=coco_add_event(prob,'UZ','r2',-0.9:0.02:0);
prob=coco_set(prob,'cont','NPR',inf,'ItMX',100, 'NAdapt', 100);
runid=['stagedbc',num2str(rot.denom)];
bd=coco(prob,runid,[],1,{'r2','b','e.r2','e.b'},{[-0.9,0]});
%%
xidx=muidx(mx,:);
figure(1);clf;
pplane=coco_bd_col(bd,{'r2','b'});
plot(pplane(1,:),pplane(2,:),'o-');
labs=sort(coco_bd_labs(bd,'UZ'));
figure(2);clf;hold on
ytor=NaN(numel(xidx),length(labs));
for i=1:length(labs)
    chart_p = coco_read_solution(runid, labs(i),'chart');
    ytor(:,i)=chart_p.x(xidx(:));
    tor=reshape(chart_p.x(xidx(:)),2,[]);
    plot(tor(1,[1:end,1]),tor(2,[1:end,1]),'b-')
    axis equal;grid on
    title(sprintf('i=%d,b=%g,r2=%g',labs(i),...
        chart_p.x(muidx(mp(2),1)),chart_p.x(muidx(mp(1),1))));
    drawnow
    %pause
end
%%
save(['staged',num2str(rot.denom),'_invc_adj.mat'])
