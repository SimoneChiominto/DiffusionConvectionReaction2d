function [convergence,approx_problems] = testConvergence(problem,refElement,RefiningOptions,options)
%TESTCONVERGENCE Summary of this function goes here
%   Detailed explanation goes here
arguments
    problem DiffusionConvectionReactionProblem2D
    refElement 
    RefiningOptions
    options.massLumping         logical = true
    options.correctConvection   logical = true
    options.plot                logical = false      
    options.fast                logical = false
end
fig=waitbar(0);
AreaValue=[0.02 0.01 0.008 0.007 0.005 0.003 0.002 0.001];% 0.0005 0.0002];
max_k=length(AreaValue);
if options.fast
    max_k=max_k-3;
end

for k=1:max_k
    waitbar(k/max_k,fig,sprintf("%d di %d",k,max_k))
    RefiningOptions.AreaValue  = AreaValue(k);
    mesh=Mesh(problem.domain,problem.BC,RefiningOptions);
    approx_problems(k)=ApproxDiffusionConvectionReactionProblem2D(problem,mesh,refElement);
    tic
    approx_problems(k).generateLinearSystem("correctConvection",options.correctConvection,"massLumping",options.massLumping);
    toc
    approx_problems(k).solve()

    Linf_error(k)=approx_problems(k).getLInfError();
    L2_error(k)=approx_problems(k).getL2Error();
    H0_error(k)=approx_problems(k).getH0Error();
    h(k)=approx_problems(k).mesh.getDiamMax();
    Area(k)=approx_problems(k).mesh.getAreaMax();
    n_dof(k)=approx_problems(k).mesh.getDoF();
end

convergence.h.inf = polyfit(log(h),log(Linf_error),1);
convergence.h.L2 = polyfit(log(h),log(L2_error),1);
convergence.h.H0 = polyfit(log(h),log(H0_error),1);

convergence.Area.inf = polyfit(log(Area),log(Linf_error),1);
convergence.Area.L2 = polyfit(log(Area),log(L2_error),1);
convergence.Area.H0 = polyfit(log(Area),log(H0_error),1);

convergence.nDof.inf = polyfit(log(n_dof),log(Linf_error),1);
convergence.nDof.L2 = polyfit(log(n_dof),log(L2_error),1);
convergence.nDof.H0 = polyfit(log(n_dof),log(H0_error),1);

if options.plot
    figure;loglog(h,Linf_error);
    figure;loglog(h,L2_error);
    figure;loglog(h,H0_error);
    figure;loglog(Area,Linf_error);
    figure;loglog(Area,L2_error);
    figure;loglog(Area,H0_error);
    figure;loglog(n_dof,Linf_error);
    figure;loglog(n_dof,L2_error);
    figure;loglog(n_dof,H0_error);
end

end

