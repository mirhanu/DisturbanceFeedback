 #Nominal MPC implementation

 mutable struct mpcController
    #Struct properties
    A::Matrix{Float64}
    B1::Matrix{Float64}
    B2::Matrix{Float64}
    E::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
    b::Matrix{Float64}
    Y::Matrix{Float64}
    z::Matrix{Float64}
    Ab::Matrix{Float64}
    B1b::Matrix{Float64}
    B2b::Matrix{Float64}
    Cb::Matrix{Float64}
    F::Matrix{Float64}
    T::Matrix{Float64}
    csmall::Matrix{Float64}
    model::Model
    n::Int
    m::Int
    s::Int
    r::Int
    l::Int
    N::Int
    feasibilityTrack::Vector{Int}
    lastv::Vector{Float64}
end
    # function definitions
     function mpcController(A, B1, B2, E, C, D, b, Y, z, N)
        n=size(A,1);
        m=size(B1,2);
        s=size(C,1);
        r=size(Y,1);
        l=size(E,1);
        model = Model(OSQP.Optimizer);
        obj = mpcController(A, B1, B2, E, C, D, b, Y, z,[;;], [;;], [;;], [;;], [;;], [;;], [;;], model, n ,m, s, r, l, N,[],[]);
        set_Matrices(obj);
        return obj
    end
    
     function gen_bigA(self::mpcController)
        n=self.n;
        N=self.N;
        A=self.A;
        Ab=zeros((N+1)*n, n);
        ai=Matrix{Float64}(I, n, n);
        for i in 1:n:n*(N+1)-1
            Ab[i:i+n-1,:]=ai;
            ai=A*ai;
        end
        self.Ab=Ab;
    end
    
     function gen_bigE(self::mpcController)
        n=self.n;
        N=self.N;
        Ab=self.Ab;
        Su=zeros((N+1)*n,n*N);
        columns=Ab[1:n*N,:];
        for i=1:n:n*N-n+1
            Su[i+n:n*(N+1),i:i+n-1]=columns[1:n*N-i+1,:];
        end
        return Su;
    end
     function set_Matrices(self::mpcController)
        n=self.n;
        N=self.N;
        m=self.m;
        gen_bigA(self);
        Eb=gen_bigE(self);
        self.Cb=[kron(Matrix{Float64}(I, N, N),self.C) zeros(N*size(self.C,1),n);zeros(size(Y,1),n*N) self.Y];
        Db=[kron(Matrix{Float64}(I, N, N),self.D); zeros(size(self.Y,1),m*N)];
        self.csmall=[kron(ones(N,1),self.b);self.z];
        self.B1b=Eb*kron(Matrix{Float64}(I, N, N),self.B1);
        self.B2b=Eb*kron(Matrix{Float64}(I, N, N),self.B2);
        self.F=self.Cb*self.B1b+Db;
        self.T=-self.Cb*self.Ab;
    end

     function construct_softVariable(self::mpcController,soft,liftedSoftInds,terminalSoftInds)
        N=self.N;
        s=self.s;
        r=self.r;
        softLifted=zeros(AffExpr,N*s+r,1);
        softCount=1;
        for i=1:length(liftedSoftInds)
            softLifted[Int.(liftedSoftInds[i])]=soft[Int.(softCount)];
            softCount=softCount+1;
        end
        for i=1:length(terminalSoftInds)
            softLifted[N*s+Int.(terminalSoftInds[i])]=soft[Int.(softCount)];
            softCount=softCount+1;
          end
        return softLifted;
    end

    

     function calc_cost(self::mpcController,states,Inps,costParams)
        N=self.N;
        m=self.m;
        n=self.n;
        # Q=10*Matrix{Float64}(I, N*n, N*n);
        # R=Matrix{Float64}(I, N*m, N*m);
        # cost=Inps'*R*Inps+states'*Q*states;
        sysCb=kron(Matrix{Float64}(I, N, N),costParams.sysC);
        sysDb=kron(Matrix{Float64}(I, N, N),costParams.sysD);
        resPresB=kron(ones(N,1),costParams.reservoirPressures);
        heads=sysCb*states[1:(N)*n]+sysDb*Inps;
        W=kron(Diagonal(vec(get_ExInpinHorizon(N,costParams.t,costParams.elecPrice))),Matrix{Float64}(I, n, n));
        # cost=transpose(vec(Inps))*W*(vec(heads)-vec(resPresB));
        cost=Inps'*W*(vec(heads)-vec(resPresB));
        # cost=Inps'*W*(sysCb*self.Ab[1:(N)*n,:]*x0+(sysCb*self.B1b[1:(N)*n,:]+sysDb)*Inps+sysCb*self.B2b[1:(N)*n,:]*vec(ExInp)-vec(resPresB));
        return cost;
    end
    
     function calc_Inp(self::mpcController,x0,ExInp,costParams,isSoft=false,softInds=[],terminalSoftInds=[];)
        n=self.n;
        m=self.m;
        N=self.N;
        s=self.s;
        self.model=Model(Ipopt.Optimizer);
        set_silent(self.model)
        @variable(self.model, v[1:m*N]);

        states=self.Ab*x0+self.B1b*v+self.B2b*ExInp';
        cost=calc_cost(self,states,v,costParams);
        if isSoft
            sC=10000000000;
            liftedSoftInds=kron(ones(1,N),(softInds)')+kron((0:s:(N-1)*s)',ones(1,length(softInds)));
            @variable(self.model,soft[1:length(liftedSoftInds)+length(terminalSoftInds)]);
            softLifted=construct_softVariable(self,soft,liftedSoftInds,terminalSoftInds);
            cost=cost+sC*sum(soft);
            @constraint(self.model, c1, self.F*v.<=  self.csmall+ self.T*x0- self.Cb* self.B2b*ExInp'+softLifted);
            inds=kron(ones(1,N),(5:self.s)')+kron((0:self.s:(N-1)*self.s)',ones(1,4));
            @constraint(self.model, c3, soft.>=0)
        else
            @constraint(self.model, c1, self.F*v.<=  self.csmall+ self.T*x0- self.Cb* self.B2b*ExInp');      
        end
        @objective(self.model, Min, cost)

        
        optimize!(self.model)
        if is_solved_and_feasible(self.model)
            push!(self.feasibilityTrack,1);
            self.lastv=JuMP.value.(v);
            return JuMP.value.(v), JuMP.value.(states[n+1:(N+1)*n]);
        else
            nSteps=length(self.feasibilityTrack);
            distInd=1;
            for i=1:nSteps
                if self.feasibilityTrack[nSteps-i+1]==1
                    distInd=i;
                    break;
                end
            end
            push!(self.feasibilityTrack,0);
            Inp=self.lastv[(distInd)*self.m+1:(distInd+1)*self.m];
            return Inp, [self.A*x0+self.B1*Inp+self.B2*ExInp[:,1]];
        end
        
    end
