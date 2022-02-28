%% Demonstration for staged construction of invariant curve and the adjoint problem for phase response
%
% Script regenerating Figure 1 in H. Dankowicz, J. Sieber, "Sensitivity
% Analysis for Periodic Orbits and Quasiperiodic Invariant Tori Using the
% Adjoint Method"
%
%%
clear
load('staged377_invc_adj.mat');
%%
dim=2;
nx=rot.denom;
%% get results from solution labels
xidx=muidx(mx,:);
r2idx=muidx(mp(1));
bidx=muidx(mp(2));
dridx=muidx(mdrho);
pplane=coco_bd_col(bd,{'r2','b'});
labs=sort(coco_bd_labs(bd,'UZ'));
ytor=NaN(numel(xidx),length(labs));
for i=length(labs):-1:1
    fprintf('label %d of %d\n',i,length(labs));
    chart_p = coco_read_solution(runid, labs(i),'chart');
    mutor(:,i)=chart_p.x(xidx(:));
    [r2(i),b(i),drho(i)]=deal(chart_p.x(r2idx),chart_p.x(bidx),chart_p.x(dridx));
end
%% get values of adjoint variables for residual M
labsel=[10,length(labs)-1];
lnames={'A','B'};
ator=NaN(dim,rot.denom,length(labsel));
for i=1:length(labsel)
    for k=rot.denom:-1:1
        chart_a=coco_read_adjoint(['M',num2str(k)],runid,labs(labsel(i)),'chart');
        ator(:,k,i)=chart_a.x;
    end
end
%%
ir=@(i)mod(i-rot.numer-1,rot.denom)+1;
xall=reshape(mutor,2,rot.denom,[]);
fac=0.5;nst=1:30:rot.denom;
ndev=5;
projmat=[0,-1;1,0];
devnorm=1e-3;
nt=6;
b_off=[0,-0.02];
lw={'linewidth',2};
ltx={'Interpreter','LaTeX','fontsize',18};
txt={'fontsize',14};
figure(1);clf;
subplot(2,3,1);ax=gca;
plot(ax,r2,b,'.-',lw{:});
hold(ax,'on');
plot(ax,r2(labsel),b(labsel),'ko','markersize',8,'MarkerFaceColor','k');
text(ax,r2(labsel)+b_off(1),b(labsel)+b_off(2),lnames,ltx{:},txt{:},...
    'HorizontalAlignment','center');
set(ax,'box','on',lw{:},txt{:});
grid(ax,'on');
xlabel(ax,'$r_2$',ltx{:});
ylabel(ax,'$b$',ltx{:});
ax2=gca;hold(ax2,'on');
clrs=lines();


%%
for i0=1:length(labsel)
    ind=labsel(i0);
    x2=xall(:,:,ind);
    GpI{i0}=Gamma(x2,@(x)dMx(x,r2(ind),b(ind)),rot);
    dx2=(rotmat(nx,1,dim)*x2(:)-x2(:))*nx;
    dx2=reshape(dx2,dim,[]);
    lam2=ator(:,ir(1:nx),i0);
    fprintf('lam in left nullspace? orthogonality in point %s: %g\n',...
        lnames{i0},norm(lam2(:).'*(GpI{i0}-eye(size(GpI{i0})))/norm(lam2)))
    lamor=projmat*lam2;
    lamor=lamor./sqrt(sum(lamor.^2,1))*fac;
    subplot(2,3,3*(i0-1)+2);ax2=gca;
    plot(ax2,[x2(1,nst)-lamor(1,nst)*fac;x2(1,nst)+lamor(1,nst)*fac],...
        [x2(2,nst)-lamor(2,nst)*fac;x2(2,nst)+lamor(2,nst)*fac],'-',lw{:},'color',clrs(2,:))
    hold(ax2,'on');
    plot(ax2,x2(1,:),x2(2,:),'r.',...
        [x2(1,:);x2(1,:)+dx2(1,:)/nx],[x2(2,:);x2(2,:)+dx2(2,:)/nx],'.-','color',clrs(1,:));
    ax2.DataAspectRatio=[1,1,1];
    grid(ax2,'on');
    set(ax2,'box','on',lw{:},txt{:});
    xlim(ax2,[-2,2]);
    ylim(ax2,[-2,2]);
    xlabel(ax2,'$x_1$',ltx{:});
    ylabel(ax2,'$x_2$',ltx{:});
    text(ax2,1.8,1.7,lnames{i0},ltx{:},txt{:},...
    'HorizontalAlignment','center')
    % test asymptotic phase
    [mx,ix]=max(abs(projmat*lam2),[],2);
    isel=ix(1);
    dev=2*(rand(dim,ndev)-0.5)*devnorm;
    proj=lam2(:,isel)'*dev;
    xtini=repmat(x2(:,isel),1,ndev,2);
    xtini(:,:,1)=xtini(:,:,1)+dev;
    xtini(:,:,2)=xtini(:,:,2)+dx2(:,isel*ones(1,ndev)).*proj(ones(dim,1),:)*nx;
    plot(ax2,xtini(1,:),xtini(2,:),'ko','MarkerFaceColor','k');
    clear xt
    xt(:,:,:,1)=xtini;
    for k=1:nt
        xt(:,:,:,k+1)=Mvec(@(x)M(x,r2(ind),b(ind)),xt(:,:,:,k));
    end
    xdev=reshape(xt(:,:,2,:)-xt(:,:,1,:),dim,ndev,nt+1);
    %
    subplot(2,3,3*(i0-1)+3);ax4=gca;
    semilogy(ax4,1:nt+1,squeeze(sqrt(sum(xdev.^2,1))),'+-','color',clrs(1,:))
    grid(ax4,'on');
    ax4.YMinorGrid='off';
    ax4.YMinorTick='off';
    set(ax4,'box','on',lw{:},txt{:});
    %ylim(ax4,[1e-8*10^(2*i0-2),1e-3]);
    xlabel(ax4,'iterate',ltx{:});
    ylabel(ax4,'$q_\mathrm{tr}$ component',ltx{:});
    text(ax4,6.5,ax4.YLim(2)/1.8/((3-i0)^0.5),lnames{i0},ltx{:},txt{:},...
        'HorizontalAlignment','center')
end
%% Spectrum in point B for Gamma and Gamma_tr
dxmat=num2cell(reshape(dx2*nx,2,1,[]),[1,2]);
lam2mat=num2cell(reshape(lam2,1,2,[]),[1,2]);
qtg=blkdiag(dxmat{:})*blkdiag(lam2mat{:});
qtr=speye(size(qtg))-qtg;
prop=GpI{2};
evproptr=eig(full(prop*qtr));
evprop=eig(full(prop));
%%
subplot(2,3,4);ax3=gca;cla(ax3);hold(ax3,'on');
pp=plot(ax3,real(evprop)-1,imag(evprop),'o','color',clrs(1,:),lw{:},...
    'DisplayName','spec $\Gamma_\rho$');
pptr=plot(ax3,real(evproptr)-1,imag(evproptr),'k.',lw{:},'Markersize',10,...
    'DisplayName','spec $\Gamma_\mathrm{tr}$');
axis(ax3,'equal');
set(ax3,'box','on',lw{:},txt{:});
legend(ax3,[pp,pptr],'Location','best',ltx{:},txt{:});
xlabel(ax3,'Re',ltx{:});
ylabel(ax3,'Im',ltx{:});
grid(ax3,'on');
%%
function xt=Mvec(M,x)
dim=size(x,1);
xv=reshape(x,dim,[]);
for i=size(xv,2):-1:1
    xt(:,i)=M(xv(:,i));
end
xt=reshape(xt,size(x));
end
function R=rotmat(nx,nrot,dim)
mrot1=diag(sparse(ones(nx-nrot,1)),nrot)+diag(sparse(ones(nrot,1)),nrot-nx);
R=kron(mrot1,speye(dim));
end
function G=Gamma(x,dxM,rot)
[dim,nx]=size(x);
assert(nx==rot.denom);
mrot=rotmat(nx,nx-rot.numer,dim);
xrot=reshape(mrot*x(:),dim,nx);
for i=size(xrot,2):-1:1
    dxMval{i}=dxM(xrot(:,i));
end
G=blkdiag(dxMval{:})*mrot;
end
