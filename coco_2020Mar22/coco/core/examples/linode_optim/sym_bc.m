function varargout=sym_bc(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=1;
   return
  case 'nout'
   varargout{1}=4;
   return
  case 'argrange'
   varargout{1}=struct('u',1:6);
   return
  case 'argsize'
   varargout{1}=struct('u',6);
   return
  case 'vector'
   varargout{1}=struct('u',1);
   return
  case 'extension'
   varargout{1}='rhs';
   return
  case 'maxorder'
   varargout{1}=2;
   return
end
nout=4;
order=varargin{1};
f=str2func(sprintf('sym_bc_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});




function [out1,out2,out3,out4] = sym_bc_rhs_0(u1,u2,u3,u4,u5,u6,u1_dev,u2_dev,u3_dev,u4_dev,u5_dev,u6_dev)
%SYM_BC_RHS_0
%    [OUT1,OUT2,OUT3,OUT4] = SYM_BC_RHS_0(U1,U2,U3,U4,U5,U6,U1_DEV,U2_DEV,U3_DEV,U4_DEV,U5_DEV,U6_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-Apr-2020 22:15:38

out1 = -u1+u3;
if nargout > 1
    out2 = -u2+u4;
end
if nargout > 2
    out3 = u5;
end
if nargout > 3
    out4 = u6-pi.*2.0;
end


function [out1,out2,out3,out4] = sym_bc_rhs_1(u1,u2,u3,u4,u5,u6,u1_dev,u2_dev,u3_dev,u4_dev,u5_dev,u6_dev)
%SYM_BC_RHS_1
%    [OUT1,OUT2,OUT3,OUT4] = SYM_BC_RHS_1(U1,U2,U3,U4,U5,U6,U1_DEV,U2_DEV,U3_DEV,U4_DEV,U5_DEV,U6_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-Apr-2020 22:15:39

out1 = -u1_dev+u3_dev;
if nargout > 1
    out2 = -u2_dev+u4_dev;
end
if nargout > 2
    out3 = u5_dev;
end
if nargout > 3
    out4 = u6_dev;
end


function [out1,out2,out3,out4] = sym_bc_rhs_2(u1,u2,u3,u4,u5,u6,u1_dev,u2_dev,u3_dev,u4_dev,u5_dev,u6_dev)
%SYM_BC_RHS_2
%    [OUT1,OUT2,OUT3,OUT4] = SYM_BC_RHS_2(U1,U2,U3,U4,U5,U6,U1_DEV,U2_DEV,U3_DEV,U4_DEV,U5_DEV,U6_DEV)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-Apr-2020 22:15:39

out1 = 0.0;
if nargout > 1
    out2 = 0.0;
end
if nargout > 2
    out3 = 0.0;
end
if nargout > 3
    out4 = 0.0;
end

