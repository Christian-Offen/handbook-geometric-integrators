# true motion of pendulum equation θ''+c*sin⁡θ=θ=0
# assuming vanishing initial velocity
# based on Abramowitz & Stegun, section 16.4, p571

# author: Christian Offen
# accompanying source code for the handbook article 
# "Solving ODEs for Nonlinear Dynamics with Symplectic and Geometric Integration" within "Handbook on Nonlinear Dynamics. Volume 2 Numerical Methods", World Scientific, editor: Vincent Acary


import Elliptic.Jacobi # based on Abramowitz & Stegun, section 16.4, p571

# solution theta
function Theta(c,theta0,t)
    m = sin(theta0/2)^2
    cd = Jacobi.cd(sqrt(c)*t,m)
    return 2*asin(sqrt(m)*cd)
end

# derivative of solution theta
function Theta_deriv(c,theta0,t)
    m = sin(theta0/2)^2
    md_n=-1+m
    tau=sqrt(c)*t
    cd=Jacobi.cd(tau,m)
    nd=Jacobi.nd(tau,m)
    sd=Jacobi.sd(tau,m)
    nominator = 2*sqrt(c)*md_n*sqrt(m)*nd*sd
    denominator = sqrt(1-m*cd^2)
    return nominator/denominator
end
