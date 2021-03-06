function varargout=sym_g(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=1;
   return
  case 'nout'
   varargout{1}=1;
   return
  case 'argrange'
   varargout{1}=struct('x',1:3);
   return
  case 'argsize'
   varargout{1}=struct('x',3);
   return
  case 'vector'
   varargout{1}=struct('x',1);
   return
  case 'extension'
   varargout{1}='rhs';
   return
  case 'maxorder'
   varargout{1}=2;
   return
end
nout=1;
order=varargin{1};
f=str2func(sprintf('sym_g_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});




function out1 = sym_g_rhs_0(x1,x2,x3,x1_dev,x2_dev,x3_dev)
%SYM_G_RHS_0
%    OUT1 = SYM_G_RHS_0(X1,X2,X3,X1_DEV,X2_DEV,X3_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    01-Mar-2021 18:53:36

out1 = x1./(x2.^2+1.0);


function out1 = sym_g_rhs_1(x1,x2,x3,x1_dev,x2_dev,x3_dev)
%SYM_G_RHS_1
%    OUT1 = SYM_G_RHS_1(X1,X2,X3,X1_DEV,X2_DEV,X3_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    01-Mar-2021 18:53:36

t2 = x2.^2;
t3 = t2+1.0;
out1 = x1_dev./t3-1.0./t3.^2.*x1.*x2.*x2_dev.*2.0;


function out1 = sym_g_rhs_2(x1,x2,x3,x1_dev,x2_dev,x3_dev)
%SYM_G_RHS_2
%    OUT1 = SYM_G_RHS_2(X1,X2,X3,X1_DEV,X2_DEV,X3_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    01-Mar-2021 18:53:36

t2 = x2.^2;
t3 = x2_dev.^2;
t4 = t2+1.0;
t5 = 1.0./t4.^2;
out1 = t3.*t5.*x1.*-2.0-t5.*x2.*x1_dev.*x2_dev.*4.0+t2.*t3.*1.0./t4.^3.*x1.*8.0;

