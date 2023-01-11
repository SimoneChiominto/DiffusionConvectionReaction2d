% 
% %%
 addpath("Classes")
 addpath("utils")
% 
% problem=problem1
 AreaValue=0.005;
 RefiningOptions.CheckArea  = 'Y';
 RefiningOptions.CheckAngle = 'Y';
 RefiningOptions.AreaValue  = AreaValue; %l'ho modificato
 RefiningOptions.AngleValue = 20;
 RefiningOptions.Subregions = [];
 load("refElements.mat")

 %RefiningOptions= %see Triangler Documentation
RefiningOptions.AreaValue  = 0.1; %max Area of the triangles
RefiningOptions.AngleValue = 20; %min angle of the triangles

%load presaved P1 and P2 reference elements
load("refElements.mat")

%create mesh using the Triangler
mesh=Mesh(problem.domain,problem.BC,RefiningOptions);

%create object describing the approximate problem
approx_problem=ApproxParabolicProblem...
                (problem,mesh,courantEl,20);
%compute stiffness matrix and forcing vector
approx_problem.computeMassMatrix()
approx_problem.computeStiffnessMatrix()

approx_problem.solve()
%plot solution
approx_problem.video()
%%
% load("test/tests.mat")
% [convergence1,approx_problems1]=testConvergence(problem1,courantEl,RefiningOptions,"fast",true);
% save("test/resultsP1.mat","approx_problems1")
% clear approx_problems1
% 
% %%
%  load("test/tests.mat")
%  [convergence1,approx_problems1]=testConvergence(problem1,courantEl,RefiningOptions,"fast",false);
%  [convergence2,approx_problems2]=testConvergence(problem2,courantEl,RefiningOptions,"fast",false);
%  [convergence3,approx_problems3]=testConvergence(problem3,courantEl,RefiningOptions,"fast",false);
%  [convergence4_nocorrection,approx_problems4_nocorrection]=testConvergence(problem4,courantEl,RefiningOptions,"fast",false,"correctConvection",false);
%  % % save("test/resultsP1.mat","approx_problems1","approx_problems2","approx_problems3","approx_problems4_nocorrection")
%  [convergence4,approx_problems4]=testConvergence(problem4,courantEl,RefiningOptions,"fast",false);
%  [convergence5_nomasslumping,approx_problems5_nomasslumping]=testConvergence(problem5,courantEl,RefiningOptions,"fast",false,"massLumping",false);
%  [convergence5,approx_problems5]=testConvergence(problem5,courantEl,RefiningOptions,"fast",false);
%  [convergence6_nocorrection,approx_problems6_nocorrection]=testConvergence(problem6,courantEl,RefiningOptions,"fast",false,"correctConvection",false);
%  [convergence6,approx_problems6]=testConvergence(problem6,courantEl,RefiningOptions,"fast",false);
%  save("test/resultsP1.mat","approx_problems1","approx_problems2","approx_problems3","approx_problems4_nocorrection","approx_problems4","approx_problems5_nomasslumping","approx_problems5","approx_problems6_nocorrection","approx_problems6")
%  clear approx_problems6 approx_problems5 approx_problems4 approx_problems6_nocorrection approx_problems5_nomasslumping approx_problems4_nocorrection approx_problems3 approx_problems2 approx_problems1
% % %%
% [P2_convergence1,P2_approx_problems1]=testConvergence(problem1,P2El,RefiningOptions,"fast",false);
% [P2_convergence2,P2_approx_problems2]=testConvergence(problem2,P2El,RefiningOptions,"fast",false);
% [P2_convergence3,P2_approx_problems3]=testConvergence(problem3,P2El,RefiningOptions,"fast",false);
% [P2_convergence4_nocorrection,P2_approx_problems4_nocorrection]=testConvergence(problem4,P2El,RefiningOptions,"fast",false,"correctConvection",false);
% [P2_convergence4,P2_approx_problems4]=testConvergence(problem4,P2El,RefiningOptions,"fast",false);
% [P2_convergence5_nomasslumping,P2_approx_problems5_nomasslumping]=testConvergence(problem5,P2El,RefiningOptions,"fast",false,"massLumping",false);
% [P2_convergence5,P2_approx_problems5]=testConvergence(problem5,P2El,RefiningOptions,"fast",false);
% [P2_convergence6_nocorrection,P2_approx_problems6_nocorrection]=testConvergence(problem6,P2El,RefiningOptions,"fast",false,"correctConvection",false);
% [P2_convergence6,P2_approx_problems6]=testConvergence(problem6,P2El,RefiningOptions,"fast",false);
% %%
% save("test/resultsP2.mat","P2_approx_problems1","P2_approx_problems2","P2_approx_problems3","P2_approx_problems4_nocorrection","P2_approx_problems4","P2_approx_problems5_nomasslumping","P2_approx_problems5","P2_approx_problems6_nocorrection","P2_approx_problems6")
% clear P2_approx_problems6 P2_approx_problems5 P2_approx_problems4 P2_approx_problems6_nocorrection P2_approx_problems5_nomasslumping P2_approx_problems4_nocorrection P2_approx_problems3 P2_approx_problems2 P2_approx_problems1
% %
% %
% load("test/test_neumann.mat")
% [convergence1,approx_problems_neumann1]=testConvergence(problem_neumann1,courantEl,RefiningOptions,"fast",false);
% [convergence2,approx_problems_neumann2]=testConvergence(problem_neumann2,courantEl,RefiningOptions,"fast",false);
% [convergence4_nocorrection,approx_problems_nocorrection_neumann2]=testConvergence(problem_neumann2,courantEl,RefiningOptions,"fast",false,"correctConvection",false);
% 
% save("test/results_neumannP1.mat","approx_problems_neumann1","approx_problems_neumann2","approx_problems_nocorrection_neumann2");
% clear approx_problems_nocorrection_neumann2 approx_problems_neumann2 approx_problems_neumann1
% %
% [P2_convergence1,P2_approx_problems_neumann1]=testConvergence(problem_neumann1,P2El,RefiningOptions,"fast",false);
% [P2_convergence2,P2_approx_problems_neumann2]=testConvergence(problem_neumann2,P2El,RefiningOptions,"fast",false);
% [P2_convergence4_nocorrection,P2_approx_problems_nocorrection_neumann2]=testConvergence(problem_neumann2,P2El,RefiningOptions,"fast",false,"correctConvection",false);
% %%
% save("test/results_neumannP2.mat","P2_approx_problems_neumann1","P2_approx_problems_neumann2","P2_approx_problems_nocorrection_neumann2");
% clear P2_approx_problems_nocorrection_neumann2 P2_approx_problems_neumann2 P2_approx_problems_neumann1
% %%
% %
% load("test/test_parabolic.mat")
% [convergence_t,approx_problems_t]=testConvergence(problem_t,courantEl,RefiningOptions,"fast",false,"time",true);
% 
% [convergence_x,approx_problems_x]=testConvergence(problem_x,courantEl,RefiningOptions,"fast",false);
% save("test/results_parabolicP1.mat","approx_problems_t","approx_problems_x")
% clear approx_problems_x approx_problems_t
%%
%[P2_convergence_t,P2_approx_problems_t]=testConvergence(problem_t,P2El,RefiningOptions,"fast",false,"time",true);
%[P2_convergence_x,P2_approx_problems_x]=testConvergence(problem_x,P2El,RefiningOptions,"fast",false);
%save("test/results_parabolicP2.mat","P2_approx_problems_t","P2_approx_problems_x")
%clear P2_approx_problems_x P2_approx_problems_t
%
%
%conv=convResultsT(P2_approx_problems_t)
%close all
% %%
 problems={approx_problems1,approx_problems2,approx_problems3,approx_problems4_nocorrection,approx_problems5_nomasslumping,approx_problems6_nocorrection}
 figure
for p=1:length(problems)
    for k=1:length(problems{p})
        Linf_error(k)=problems{p}(k).getLInfError();
        L2_error(k)=problems{p}(k).getL2Error();
        H0_error(k)=problems{p}(k).getH0Error();
        h(k)=problems{p}(k).mesh.getDiamMax();
        Area(k)=problems{p}(k).mesh.getAreaMax();
        n_dof(k)=problems{p}(k).mesh.getDoF();
        c(k)=condest(problems{p}(k).stiffnessMatrix);
    end


%convergence.h.inf = polyfit(log(h),log(Linf_error),1);
%convergence.h.L2 = polyfit(log(h),log(L2_error),1);
%convergence.h.H0 = polyfit(log(h),log(H0_error),1);

%convergence.Area.inf = polyfit(log(Area),log(Linf_error),1);
%convergence.Area.L2 = polyfit(log(Area),log(L2_error),1);
%convergence.Area.H0 = polyfit(log(Area),log(H0_error),1);

%convergence.nDof.inf = polyfit(log(n_dof),log(Linf_error),1);
%convergence.nDof.L2 = polyfit(log(n_dof),log(L2_error),1);
%convergence.nDof.H0 = polyfit(log(n_dof),log(H0_error),1);

 
%loglog(h,Linf_error);
%figure;
%loglog(h,L2_error);
%loglog(h,H0_error);
%figure;loglog(Area,Linf_error);
%figure;
%loglog(Area,L2_error);
%figure;loglog(Area,H0_error);
%figure;loglog(n_dof,Linf_error);
%figure;loglog(n_dof,L2_error);
%figure;loglog(n_dof,H0_error);
polyfit(log(h),log(c),1)
loglog(h,c);
hold on
end
legend("approx_problems1","approx_problems2","approx_problems3","approx_problems4","approx_problems5","approx_problems6","Location","northwest",'FontSize',12)
xlabel("h",'FontSize',15)
ylabel("||u-u_h||_{L^2(\Omega)}",'FontSize',15)
% 
% 
% 
% 
% 
% 
function convergence= convResults(approx_problems)
for k=1:length(approx_problems)
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

function convergence= convResultsT(approx_problems)
for k=1:length(approx_problems)-2
    Linf_error(k)=approx_problems(k+2).getLInfError();
    L2_error(k)=approx_problems(k+2).getL2Error();
    H0_error(k)=approx_problems(k+2).getH0Error();
    dt(k)=(approx_problems(k+2).timeDiscretization(2)-approx_problems(k).timeDiscretization(1))
end


convergence.dt.inf = polyfit(log(dt),log(Linf_error),1);
convergence.dt.L2 = polyfit(log(dt),log(L2_error),1);
convergence.dt.H0 = polyfit(log(dt),log(H0_error),1);

%figure;loglog(dt,Linf_error);
figure;loglog(dt,L2_error);hold on
loglog(dt,H0_error);
legend("||u(T)-u_h(T)||_{L^2(\Omega)}","|u(T)-u_h(T)|_{H^0(\Omega)}","FontSize",12,"Location","northwest")
 xlabel("dt",'FontSize',15)
 ylabel("Error",'FontSize',15)
%figure;loglog(dt,H0_error);

end



