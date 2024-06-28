function get_ExInpinHorizon(N,t,ExInp)
    ds=size(ExInp,1);
    ExInpinHorizon=zeros(ds,N);
    period=size(ExInp,2);
    for i=0:N-1
        ind=mod(t+i,period)+1;
        ExInpinHorizon[:,i+1]=ExInp[:,ind];
    end
    return ExInpinHorizon;
end

function calc_cost(states,inps,costParams)
    N=size(inps,2);
    n=size(states,1);
    resPresB=kron(ones(1,N),costParams.reservoirPressures);
    heads=costParams.sysC*states[:,1:N]+costParams.sysD*inps;
    pow=sum((heads-resPresB).*inps,dims=1);
    cost=sum(pow.*get_ExInpinHorizon(N,costParams.t,costParams.elecPrice));
    return cost;
end

function sigmoid(z)
    return 1.0 ./ (1.0 .+ exp.(-z))
end

 mutable struct WDNCostParams
    sysC
    sysD
    elecPrice
    reservoirPressures
    t
end
# costParam.t=0;
# a=calc_cost2(s,v,costParam)
# b=calc_cost2(s2,v2,costParam)