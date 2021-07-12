using JuMP, CPLEX, LinearAlgebra, Ipopt

#example problem

#A = [ -3 2 4 ; 3 1 -1; 1 -2 9; -1 0 0; 0 -1 0; 0 0 -1]
#b = [9 ; 8 ; 10; 0; 0; 0]

#c = [1 ; -2; 3]
n,m = size(A)
b=zeros(n)
for i = 1:n
    b[i] = -1+sum(A[i,:].>0)
end
alpha = 0.01
b = b .+ (1 - alpha)



c = zeros(m)

function Vbounds(A,c)
    n,m = size(A)

    minV=-Inf*ones(m)
    maxV=Inf*ones(m)

    for row in eachrow(A)
        if sum(row) !=0
            rowSum = row/sum(row)
            for index in eachindex(rowSum)
                if rowSum[index] < 0
                    if -rowSum[index] < maxV[index]
                        maxV[index]=-rowSum[index]
                    end
                elseif rowSum[index] > 0
                    if -rowSum[index] > minV[index]
                        minV[index]=-rowSum[index]
                    end
                end
            end
        end
    end

    if sum(c) !=0
        rowSum = c/sum(c)
        for index in eachindex(rowSum)
            if rowSum[index] < 0
                if -rowSum[index] < maxV[index]
                    maxV[index]=-rowSum[index]
                end
            elseif rowSum[index] > 0
                if -rowSum[index] > minV[index]
                    minV[index]=-rowSum[index]
                end
            end
        end
    end


    return minV, maxV
end






# maximise c'x
# while Ax<=b

function orientationMIP(A,b,c)
    n,m = size(A)
    orientationModel = Model(CPLEX.Optimizer)

    @variable(orientationModel, 0 <= x[1:m] <= 1)

    @variable(orientationModel, u[1:m], Bin)

    @objective(orientationModel, Max, c' * x)

    @constraint(orientationModel, mainConstraint1, A * (x+u.-(1/2)) .<= b)

    @constraint(orientationModel, mainConstraint2, A * (x-u.+(1/2)) .<= b)

    TT = stdout # save original STDOUT stream
    redirect_stdout()

    optimize!(orientationModel)

    redirect_stdout(TT)

    return value.(x), 2*(value.(u) .-1/2)
end



#u=[1,1,1]
function myMethodUV(A,b,c,orientation,minV,maxV)
    n,m = size(A)

    #signCombs = unique(sign.(A), dims =2)
    #NosignCombs = size(signCombs,1)
    myMethodModel = Model(Ipopt.Optimizer)
    @variable(myMethodModel, 0<= x[1:m] <= 1)

    @variable(myMethodModel, u[1:m] >=0)

    @variable(myMethodModel, v[1:m])

    @variable(myMethodModel, q[1:m] <=0)

    @variable(myMethodModel, r[1:m] >=0)

    @variable(myMethodModel, Q )

    @variable(myMethodModel, R)

    @variable(myMethodModel, V)

    @variable(myMethodModel, y[1:m])

    @variable(myMethodModel, Y[1:m])




    #@objective(myMethodModel, Max,
    # c'*(x-(1/2)*(sign.(c).*u.+sum(u.*v.*sign.(c))))
    #)

    @objective(myMethodModel, Max,
        c'*(x-(1/2)*(sign.(c).*u.+sum(y.*sign.(c))))
    )
    #@constraint(myMethodModel, MainCon[i=1:n],
    #    A[i,:]'*(x+(1/2)*(sign.(A[i,:]).*u.+sum(u.*v.*sign.(A[i,:])))) <= b[i]
    #)

    @constraint(myMethodModel, MainCon[i=1:n],
        A[i,:]'*(x+(1/2)*(sign.(A[i,:]).*u.+sum(y.*sign.(A[i,:])))) <= b[i]
    )


    @NLconstraint(myMethodModel, y1[j=1:m],
        y[j]-u[j]*v[j]==0
    )

    @NLconstraint(myMethodModel, y2[j=1:m],
        Y[j]-u[j]*V==0
    )

    @constraint(myMethodModel, orientation11[j=1:m],
        Q-q[j]+orientation[j]*(v[j]-V)>=-Y[j]
    )

    @constraint(myMethodModel, orientation12[j=1:m],
        R-r[j]+orientation[j]*(v[j]-V)<=Y[j]
    )

    @constraint(myMethodModel, orientation21[j=1:m],
        Q-q[j]>=-Y[j]
    )

    @constraint(myMethodModel, orientation22[j=1:m],
        R-r[j]<=Y[j]
    )

    @constraint(myMethodModel, vBoundsUpper[j=1:m],
        v[j]<=0
    )

    for j=1:m
        set_normalized_rhs(vBoundsUpper[j], maxV[j])
    end

    @constraint(myMethodModel, vBoundsLower[j=1:m],
        -v[j]<=0
    )

    for j=1:m
        set_normalized_rhs(vBoundsLower[j], -minV[j])
    end

    @constraint(myMethodModel, qBound,
        q .<= orientation.*v
    )

    @constraint(myMethodModel, rBound,
        r .>= orientation.*v
    )

    @constraint(myMethodModel, qSum,
        Q == sum(q)
    )

    @constraint(myMethodModel, rSum,
        R == sum(r)
    )

    @constraint(myMethodModel, vSum,
        V == sum(v)+1
    )

    for j = 1:m
        set_start_value(u[j],1)
        set_start_value(v[j],0)
    end
    #TT = stdout # save original STDOUT stream
    #redirect_stdout()
    #set_optimizer_attribute(myMethodModel, "CPX_PARAM_EPRHS", 1e-4)

    optimize!(myMethodModel)

    #println(termination_status(myMethodModel))
    #redirect_stdout(TT)

    return objective_value(myMethodModel), value.(x), value.(u), value.(v)

    #if termination_status(myMethodModel) == MOI.OPTIMAL
    #    return objective_value(myMethodModel), value.(x), value.(v)
    #else
    #    return NaN, ones(m)*NaN, ones(m)*NaN
    #end
end



function intFinder(c,initX,u,v)

    m=length(initX)
    #normMatrix = inv((u .* v)' .+I(m).*u)


    intFinderModel = Model(CPLEX.Optimizer)

    @variable(intFinderModel, x[1:m], Bin)

    @variable(intFinderModel, s[1:m] >=0)

    @variable(intFinderModel, S)

    @objective(intFinderModel, Max, c'*x)

    @constraint(intFinderModel, sMaxCon[j=1:m], s[j]<=u[j])

    @constraint(intFinderModel, sSumCon, S==sum(s.*v))

    @constraint(intFinderModel, mainCon[j=1:m],
    x[j]+s[j]+S==initX[j] - 1/2*(u[j]+sum(u.*v))
    )
    #@constraint(intFinderModel, LowerCon,normMatrix*x.>=-ones(m).*1/2 + normMatrix*initX)

    #@constraint(intFinderModel, UpperCon,normMatrix*x.<=ones(m).*1/2 + normMatrix*initX)


    optimize!(intFinderModel)

    if termination_status(intFinderModel) == MOI.OPTIMAL
        return objective_value(intFinderModel), value.(x)
    else
        return NaN, ones(m)*NaN
    end
end

#trueObjVal = zeros(size(objValCompact,1))
#trueXVal = zeros(size(objValCompact,1),m)


#@showprogress for i =1:size(objValCompact,1)
#    trueObjVal[i], trueXVal[i,:]= intFinder(c,xValCompact[i,:],uCombCompact[i,:],vValCompact[i,:],m)

#end




function ActualSol(A,b,c)
    actualModel = Model(CPLEX.Optimizer)
    @variable(actualModel, x[1:size(A,2)], Bin)
    @objective(actualModel, Max, c'*x)
    @constraint(actualModel,generic, A*x.<=b)
    optimize!(actualModel)
    return value.(x), objective_value(actualModel)
end

#my method
minV, maxV = Vbounds(A,c)

flotSol, orientation = orientationMIP(A,b,c)

objVal, xVal, uVal, vVal =myMethodUV(A,b,c,orientation,minV,maxV)

trueObjVal, trueXVal = intFinder(c,xVal,uVal,vVal)


#actual optimal
actualXVal, actualObjVal = ActualSol(A,b,c)
