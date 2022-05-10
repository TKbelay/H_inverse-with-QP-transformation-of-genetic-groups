#run as julia GEUPG.jl (#for julia version 1.6.1)

#A) Constructed inverse pedigree based relationship matrix

using DelimitedFiles, LinearAlgebra, SparseArrays, Statistics, Serialization, ProgressMeter, DataFrames

function Amax()
    genot=trunc.(Int64,readdlm("genotyped_renum.dat"))[:];#Read a list of genotyped animal ids (151437)
    ngenot=trunc.(Int64,readdlm("non-genotyped_renum.dat"))[:]; #Read a list of non-genotyped animal ids (4790606)
    ped=trunc.(Int64,readdlm("renumped_21_nobd.ped")); 
    animals=length(ped[:,1]);
    IT=sparse(1.0I, animals, animals); #inverse of lower triangular matrix (T) in A-matrix: (inv(A)=(inv(T))'*inv(D)*inv(T)
    ID=2*ones(animals);#calculate inverse of diagonal matrix(D), which is reciprocal of diagonal elements of A-matrix

    for i=1:animals
    far=ped[i,2];
    mor=ped[i,3];
    if (far*mor)==0
    ID[i]=4/3;
    end
    if (far+mor)==0
    ID[i]=1;
    end
    #elements in inverse D are finshed;
    if far!=0
    IT[i,far]=-0.5;
    end
    if mor!=0
    IT[i,mor]=-0.5;
    end
    #elements in inverse T are finshed;
    end

    Ainv=transpose(IT)*sparse(Diagonal(ID))*IT;#can also be IT'spdiagm(ID)IT
    serialize("Ainv.ser", Ainv)
    #Ainv = deserialize("Ainv.ser");
    A11=Ainv[ngenot,ngenot];#inverse of A for non-genotyped animals
    A12=Ainv[ngenot,genot];
    A21=Ainv[genot,ngenot];
    A22=Ainv[genot,genot];
    J1=-A11\(A12*ones(151437));
    J2=ones(151437);
    J=[J1
      J2];
    Ids=[ngenot
        genot];
    J_id=[Ids J];
    writedlm("J_valus.dat",J_id,"\t")
    
    J1=0;
    J2=0;
    Ids=0;
    J=0;
    Ainv=0;
end
#calculate inverse of pedigree relationship between genotyped animals (A22inv) as A22 - A21*inv(A11)*A12 (blockwise)

function A22max(; nblk = 100)
    @info "Reading the data and genotyped ID"
    Ai = deserialize("Ainv.ser");
    jx = begin
    tmp = Int[]
    for i in eachline("genotyped_renum.dat") #renumbered genotyped ids
    push!(tmp, parse(Int, i))
    end
    tmp
    end
    N = size(Ai, 1)
    n2 = length(jx)
    ix = sort(collect(setdiff(Set(1:N), Set(jx))))
    A11 = Ai[ix, ix]
    A12 = Ai[ix, jx]
    A22 = Matrix(Ai[jx, jx])
    nr = Int(floor(n2/nblk))
    stops = collect(nr:nr:n2)
    stops[end] = n2
    start = 1
    @info "Calculating by blocks"
    @showprogress for stop in stops
    println(stop)
    rhs = Matrix(A12[:, start:stop])
    sol = A11 \ rhs
    blk = A12'sol
    k = 1
    for i in start:stop
    A22[:, i] -= blk[:, k]
    k += 1
    end
    start = stop + 1
    end
    serialize("A22inv.ser", A22)
end

function AinvQ()

    #B)correct Ainv for genetic groups (Q) contributions (QP transformation of Q) and compile the sub-matrixes into large matrix (aka Ainv_mod)

    ggrc=readdlm("kgmilk305d_ggrc_bothid_sort.dat"); # 4942043x119 Read genetic group contribution for all animals in pedigree sorted based renum id
    Q=ggrc[:,3:end]; # gentic group contributions that start at 3rd column for all animals in pedigree (4942043x117)
    #construct the sub-matrixes
    Ainv=deserialize("Ainv.ser");
    aq=sparse(-Ainv*Q);
    qa=sparse(-Q'*Ainv);
    qaq=sparse(Q'*Ainv*Q) + sparse(1.0I, 119, 119);#Q fitted as random considering identity relationships between groups
    #compile the sub-matrixes
    top=[Ainv aq];
    bot=[qa qaq];
    Ainv_mod=[top
              bot];
    serialize("Ainv_mod.ser", Ainv_mod) 
    ped=0;
    ggrc=0;
    Q=0;
    Ainv=0;
    aq=0;
    qa=0;
    qaq=0;
    top=0;
    bot=0;
end

function Gmax()
    #C) Construct G matrix without correcting it for genetic group contributions(Q)
    X=readdlm("genotype50k_renum.dat"); #read genotype file for genotyped animals(30729 genotyped animals) 
    P=0.5*mean(X[:,2:end], dims=1); # calculate allele frequencey per column
    X_cent=X[:,2:end].-(2*P); # center genotype data 
    rho=2*sum(P.*(1 .- P)); # calculate sum of heterozygozity (i.e. 2 * sum product of p*q)
    S=sparse(1.0I, 151437, 151437)*0.01;
    G=(X_cent*X_cent')/rho + S; #sparse(1.0I, 30729, 30729)*0.01; #scale and add 0.01 to diagonal of G-matrix 
    Ginv=sparse(inv(G));
    X=0;
    X_cent=0;
end

function Zmax()
    #D)construct index matrix Z of size mxn (toal number of animals in pedigree including genetic group levels by number of genotyped animals),this matrix help to correctly located genotyped position in MME

    Is=genot;
    Js=collect(1:151437);
    Vs=ones(151437);
    m=4942043;
    n=151437;
    Z=sparse(Is,Js,Vs,m,n);
    Is=0;
    Js=0;
    Vs=0;
    m=0;
    n=0;
end
    
function Hmax(Ginv,Z)    
    A22inv = deserialize("A22inv.ser"); #already constructed above using function A22max
    D=sparse(Ginv-A22inv);
    Ginv=0;
    A22inv=0;
    B=Z*D*Z';
    Z=0;
    D=0;

    Ainv = deserialize("Ainv.ser");
    Hinv=Ainv + B;
    Ai=0;
    B=0;
end
    #write out the Hinv (original ids not used) with suitable format(Id1, Id2,relationship)(only non-zero values written out)

function prtsparse(Hinv)
    df = begin
    x, y, z = findnz(Hinv)
    tmp = DataFrame(x=x, y=y, z=z)
    sort!(tmp, [:x, :y])
    tmp
    end
    open("Hinv.dat", "w") do io
    for (i, j, v) in eachrow(df)
    (i < j) && continue
    println(io, "$i $j $v")
    end
    end
end

function Allmax()
   Amax()
   A22max(nblk = 1000)
   AinvQ(Ainv)
   Gmax()
   Zmax()
   Hmax(Ginv,Z)
   prtsparse(Hinv)
end
