def cat(L) {
    while ( L != [] ) {
        print(car(L),0);
        L = cdr(L);
    }
    print(""); 
}

/*
    if Poly is not Dpoly then convert to Dpoly.
*/
def chkPoly(Poly, Vars){
    if (9 != type(Poly)){
        if(2 == type(Poly)){
            Poly = dp_ptod(Poly, Vars); // 
        }
        else{
            print(Poly, 0); error(" is not a Polynomial");
        }
    }
    return Poly;
}

/*
    make new order defined as [Vector, [Matrix]]
*/
def new_order(Vector, Matrix){
    List = [vtol(Vector)];
    N = size(Matrix)[0];
    for(I = 0; I < N; I++){
        List = cons(vtol(Matrix[I]), List);
    }
    List = reverse(List);
    Result = newmat(length(List), length(Vector), List);

    return Result;
}

/*
    return innerproduct of V1 and V2
*/
def innerProduct(V1, V2){
    if ( V1 == 0 || V2 == 0 ) return 0;
    N = length(V1);
    if(N != length(V2)){
        error("V1 and V2 have different lengths");
    }
    Sum = 0;
    for(I = 0; I < N; I++){
        Sum += V1[I]*V2[I];
    }
    return Sum;
}

/*
    return w-degree of Poly
*/
def wdeg(Poly, Weight){
    Result = 0;
    for( ; Poly; Poly = dp_rest(Poly)){
        Deg = innerProduct(Weight, dp_etov(Poly));
        if(Deg > Result) Result = Deg;
    }

    return Result;
}

/*
    return w-initial form of Poly
*/
def initialForm(Poly, Weight){
    InitialForm = HT = dp_hm(Poly);
    Deg = wdeg(HT, Weight);
    Poly = dp_rest(Poly);
    for(; Poly; Poly = dp_rest(Poly)){
        HT = dp_hm(Poly);
        if(Deg != wdeg(HT, Weight)) break;
        InitialForm += HT;
        //if(Deg == wdeg(HT, Weight)) InitialForm += HT;
    }

    return InitialForm;
}

/*
    return abolute value of N
*/
def abs_r(N){
    if(N < 0) return -N;
    return N;
}

/*
    make set list of L.
    same with remove duplicate of L and sort.
*/
def set_of(L){
    if(5 == type(L)) L = vtol(L);
    if([] == L) error("set_of : L is empty");
    L = qsort(L);
    Result = [L[0]];
    L = cdr(L);
    while(L != []){
        if(L[0] != Result[0]) Result = cons(L[0], Result);
        L = cdr(L);
    }
    return reverse(Result);
}

/*
    return LCM of list V.
    every elements of V must be non-negative.
*/
def lcm_of(V){
    V = reverse(set_of(V));
    N = length(V);
    if (0 == N) error("lcm_of: V is empty.");
    LCM = V[0];
    for(I = 1; I < N; I++){
        if(1 - (V[I] && (1 - V[I]))) return LCM;
        LCM = ilcm(LCM, V[I]);
    }
    return LCM;
}

/*
    return GCD of list V.
    every elements of V must be non-negative.
*/
def gcd_of(V){
    V = set_of(V);
    N = length(V);
    if (0 == N) error("gcd_of: V is empty.");
    if(1 - (V[0] && (V[0] - 1))) return 1;
    GCD = V[0];
    for(I = 1; I < N; I++){
        GCD = igcd(GCD, V[I]);
        if(1 == GCD) return 1;
    }
    return GCD;
}

/**/
/*
    convert rational vector V to integer vector every coordinates are relatively prime
*/
def conv2int(V){
    if(4 == type(V)){ // if V is list then change to vector
        V = ltov(V);
    }
    if(5 != type(V)){ 
        error(" V is not a vector.");
    }
    Denom_V = qsort(map(dn, V));
    LCM = lcm_of(Denom_V);
    V *= LCM;
    Numer_V = qsort(map(nm, V));
    Numer_V = map(abs_r, Numer_V);
    GCD = gcd_of(Numer_V);
    V /= GCD;
    return V;
}

/*
    make a list of monomials of a polynomial V
*/
def terms(V){
    Result = [];
    while(V != 0){
        Result = cons(dp_hm(V), Result);
        V = dp_rest(V);
    }
    return reverse(Result);
}





/*
    subroutine of nextcone(compute_last_w).
    Max is degree vector of largest term of Poly.
    return list of [Max - "degree vector of other terms of Poly"]
*/
def compute_df(Poly){
    HT = dp_ht(Poly);
    Poly -= HT;
    HT = dp_etov(HT);
    for(Result = []; Poly; Poly = dp_rest(Poly)){
        Result = cons(HT - dp_etov(Poly), Result);
    }
    return Result;
}



/*
    compute last_w for grwalk.
*/
def compute_last_t(G, W, W_t){
    DF = [];
    N = length(G);
    for(I = 0; I < N; I++){
        DF = append(compute_df(G[I]), DF);
        DF = set_of(DF);
    }
    N = length(DF);
    U_j = U_last = 1;
    for(I = 0; I < N; I ++){
        WtB = innerProduct(W_t, DF[I]);
        if(WtB < 0){
            WB = innerProduct(W, DF[I]);
            U_j = WB/(WB - WtB);
            if(U_j < U_last){
                U_last = U_j;
            }
        }
        
    }
    return U_last;
}

/*
    G is list of polynomials.
    P is a polynomial in ideal of G.
    make list of quotient divide P by G. 
    for list L, L[I] is coefficient polynomial of G[I] with representation of G.
*/
def make_quotient(P, G){
    N = length(G);
    Quotients = newvect(N);
    while(P){
        for(I = 0; I < N; I++){
            if(dp_redble(P, G[I])){
                Q = (dp_hc(P)/dp_hc(G[I]))*dp_subd(P, G[I]);
                Quotients[I] += Q;
                P = P - Q*G[I];
                if(P == 0) break;
            }
        }
    }
    return vtol(Quotients);
}

/*
    for In-representation of InG, replace In[I] to G[I]
*/
def lifting(In, InG , G, Order){
    N = length(In);
    M = length(InG);
    H = newvect(M);
    Order_new = dp_ord();
    dp_ord(Order);
    In = map(dp_sort, In);
    InG = map(dp_sort, InG);
    for(I = 0; I < M; I++){
        Quotient = make_quotient(InG[I], In);
        for(J = 0; J < N; J++){
            H[I] += Quotient[J]*G[J];
        }
    }
    dp_ord(Order_new);
    H = map(dp_sort, H);
    return H;
}

/*
    with G, compute remainder of P 
    need to fix
*/
def reducing(P, G){
    M = length(G);
    T = terms(P);
    N = length(T);
    for(I = 0; I < N; I++){
        for(J = 0; J < M; J++){
            if((T[I] != 0) && (dp_redble(T[I], G[J]))){
                P -= (dp_hc(T[I])/dp_hc(G[J]))*dp_subd(T[I], G[J])*G[J];
                break;
            }
        }
        if(P == 0) return 0;
    }
    return P;
}

/*
    for non-reduce GB G. return reduce GB respect to Order
*/
def interreduce2(G, Order){
    dp_ord(Order);
    G = map(dp_sort, G);
    G = vtol(G);
    N = length(G);
    Result = [];
    for(I = 0; I < N; I++){
        P = car(G);
        G = cdr(G);
        Q = reducing(P, append(Result, G));
        //if(Q){
            Result = cons(Q, Result);
        //}
    }
    return Result;
}

/*
    for non-reduce GB G. return reduce GB respect to Order
*/
def interreduce(G, Order){  //Slow
    Result = [];
    dp_ord(Order);
    G = map(dp_sort, G);
    N = length(G);

    for(I = 0; I < N; I++){
        T = terms(G[I]);
        M = length(T);
        for(J = 0; J < N; J++){//루프 순서 고칠것. G 정렬 함수도 만들것.
            if(J != I){
                for(K = 0; K < M; K++){
                    if(dp_redble(T[K], G[J])){
                        G[I] -= (dp_hc(T[K])/dp_hc(G[J]))*dp_subd(T[K], G[J])*G[J];
                    }
                }
            if(G[I] == 0) break;
            } 
        }
        
        if(G[I] != 0){
            Result = cons(G[I], Result);
        }
    }
    return reverse(Result);
}

end$