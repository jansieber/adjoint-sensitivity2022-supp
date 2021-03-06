function varargout=sym_bistable(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=3;
   return
  case 'nout'
   varargout{1}=2;
   return
  case 'argrange'
   varargout{1}=struct('t',1:1,'x',2:3,'p',4:6);
   return
  case 'argsize'
   varargout{1}=struct('t',1,'x',2,'p',3);
   return
  case 'vector'
   varargout{1}=struct('t',0,'x',1,'p',1);
   return
  case 'extension'
   varargout{1}='rhs';
   return
  case 'maxorder'
   varargout{1}=2;
   return
end
nout=2;
order=varargin{1};
f=str2func(sprintf('sym_bistable_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});




function [out1,out2] = sym_bistable_rhs_0(t,x,v,T,a,gam,t_dev,x_dev,v_dev,T_dev,a_dev,gam_dev)
%SYM_BISTABLE_RHS_0
%    [OUT1,OUT2] = SYM_BISTABLE_RHS_0(T,X,V,T,A,GAM,T_DEV,X_DEV,V_DEV,T_DEV,A_DEV,GAM_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    01-Mar-2021 20:31:16

out1 = v;
if nargout > 1
    t2 = -x;
    out2 = t2-gam.*v+t2.*x.^2+a.*cos((t.*pi.*2.0)./T);
end


function [out1,out2] = sym_bistable_rhs_1(t,x,v,T,a,gam,t_dev,x_dev,v_dev,T_dev,a_dev,gam_dev)
%SYM_BISTABLE_RHS_1
%    [OUT1,OUT2] = SYM_BISTABLE_RHS_1(T,X,V,T,A,GAM,T_DEV,X_DEV,V_DEV,T_DEV,A_DEV,GAM_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    01-Mar-2021 20:31:16

out1 = v_dev;
if nargout > 1
    t2 = 1.0./T;
    t3 = t.*t2.*pi.*2.0;
    out2 = -x_dev-gam_dev.*v-gam.*v_dev+a_dev.*cos(t3)-x.^2.*x_dev.*3.0-a.*sin(t3).*(t2.*t_dev.*pi.*2.0-T_dev.*t.*t2.^2.*pi.*2.0);
end


function [out1,out2] = sym_bistable_rhs_2(t,x,v,T,a,gam,t_dev,x_dev,v_dev,T_dev,a_dev,gam_dev)
%SYM_BISTABLE_RHS_2
%    [OUT1,OUT2] = SYM_BISTABLE_RHS_2(T,X,V,T,A,GAM,T_DEV,X_DEV,V_DEV,T_DEV,A_DEV,GAM_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    01-Mar-2021 20:31:16

out1 = 0.0;
if nargout > 1
    t2 = 1.0./T;
    t3 = t2.^2;
    t4 = t.*t2.*pi.*2.0;
    t5 = t2.*t_dev.*pi.*2.0;
    t6 = T_dev.*t.*t3.*pi.*2.0;
    t7 = sin(t4);
    t8 = -t6;
    t9 = t5+t8;
    out2 = gam_dev.*v_dev.*-2.0-x.*x_dev.^2.*6.0-a.*t7.*(T_dev.^2.*t.*t2.^3.*pi.*4.0-T_dev.*t3.*t_dev.*pi.*4.0)-a.*t9.^2.*cos(t4)-a_dev.*t7.*t9.*2.0;
end

