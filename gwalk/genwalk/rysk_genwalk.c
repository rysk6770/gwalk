load("rysk.c")$
load("noro_gw.rr")$


extern Count;
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
    make list of degree vectors of F
*/
def supp(F){
    V = [];
    for(; F; F = dp_rest(F)){
        V = cons(dp_etov(F), V);
    }
    return reverse(V);
}

/*
    U is maximum vector respect to order of supp(F)
    return list of del_<F := {U - U'| U' \in supp(F)-{U}}
*/
def round_poly(F, U){
    D = [];
    F = supp(F);
    for(; F != []; F = cdr(F)){
        T = car(F);
        if(T != U){
            D = cons(U - T, D);
        }
    }
    return D;
}

/*
    G is list of polynomials
    GM is list of makred term of G
    return union of round_poly for each(G[I],GM[I])
*/
def round_set(G, GM){
    D = [];
    for(; G!=[]; G = cdr(G), GM = cdr(GM)){
        D = append(D, round_poly(car(G), dp_etov(car(GM))));
    }
    D = set_of(D);
    return D;
}

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
    compare vectors V1 and V2 by weight(or matrix) M.
    if V1 < V2 then return 1
    if V1 > V2 then return -1
    if V1 == V2 then return 0
*/
def compare(V1, V2, M){
    // when M is vector.
    if(5 == type(M)){
        T = innerProduct(V1 - V2, M);
        if(T < 0) return 1;
        else if(T > 0) return -1;
        else return 0;
        //T1 = innerProduct(V1, M);
        //T2 = innerProduct(V2, M);
        //if(T1 < T2) return 1;
        //else if(T1 > T2) return -1;
        //else return 0;
    }
    // when M is matrix.
    if(6 == type(M)){
        N = size(M)[0];
        for(I = 0; I < N; I++){
            T = innerProduct(V1 - V2, M[I]);
            if(T < 0) return 1;
            else if(T > 0) return -1;
            //T1 = innerProduct(V1, M[I]);
            //T2 = innerProduct(V2, M[I]);
            //if(T1 < T2) return 1;
            //else if(T1 > T2) return -1;
        }
        return 0;
    }
}

/*
    U, V is vector
    M1, M2 if matrix
    return 1 when U < V as facet preorder defined as M1, M2
    else return 0
*/
def facet_preorder(U, V, M1, M2){
    if ( U == 0 ) return 1;
	if ( U == 1 ) return 0;
	N = size(V)[0];
	for ( I = 0; I < N; I++ ) {
		Left = innerProduct(M2[I],U)*V;
		Right = innerProduct(M2[I],V)*U;
		R = compare(Left,Right,M1);
		if (R != 0) return (R > 0) ? 1 : 0;
	}
	return 1;
}
/*
    check V is in Cone_<1,<2
    M1 = O(<1), M2 = O(<2)
*/
def in_C12(V, M1, M2){
    T1 = compare(0, V, M1);
    T2 = compare(0, V, M2);
    return (T1 == 1 && T2 == -1) ? 1 : 0;
}

/*
    G is list of polynomials
    GM is marked terms of G
    W is previous W
    M1, M2 is matrix s.t. define order
*/
def compute_last_w(G, GM, W, M1, M2){
    
    R = round_set(G, GM);
    V = [];
    for(; R != []; R = cdr(R)){
        if(in_C12(car(R), M1, M2) && facet_preorder(W, car(R), M1, M2)) V = cons(car(R), V);
    }
    if(V == []) return 1;

    Min = car(V);
    V = cdr(V);
    for(; V!=[]; V=cdr(V)){
        if(facet_preorder(car(V), Min, M1, M2)){
            Min = car(V);
        }
    }
    
    return Min;
}

/*
    GM is marked term of G
    M1 and M2 define facet preorder
*/
def initial_W(P, PM, W, M1, M2){
    P -= PM;
    U = dp_etov(PM);
    
    for(; P; P = dp_rest(P)){
        HM = dp_hm(P);
        T = U - dp_etov(HM);
        if( facet_preorder(T, W, M1, M2) && facet_preorder(W, T, M1, M2)){
            PM += HM;
        }
    }
    
    return PM;
}

def reducing_mark(P, G, GM){
    M = length(G);
    while(P){
        for(T=P; T; T= dp_rest(T)){
            for(J = 0; J < M; J++){
                if((dp_redble(dp_hm(T), GM[J]))){
                    P -= (dp_hc(T)/dp_hc(GM[J]))*dp_subd(T, GM[J])*G[J];
                    
                    break;
                }
            }
            if(J != M) break;
        }
        
        if(!T) return P;
    }
    return P;
}

def removeKth(L, K){
    R = [];
    N = length(L);
    for(I=0;I<N;I++){
        if(I!=K) R = cons(L[I], R);
    }
    return reverse(R);
}

def interreduce_mark(G, GM){
    N = length(G);
    Result = ResultM = [];
    for(I = 0; I < N; I++){
        P = G[I];
        Q = reducing_mark(P, removeKth(G, I), removeKth(GM, I));
        if(Q){
            Result = cons(Q, Result);
            ResultM = cons(GM[I], ResultM);
        }
    }
    return [Result, ResultM];
}

def lifting_mark(H, G, GM){
    N = length(H);
    Result = newvect(N);
    for(I=0; I<N; I++){
        Result[I] = H[I] - reducing_mark(H[I], G, GM);
    }
    return Result;
}



def generic_walk_r(G, StartOrder, TargetOrder, Vars){
    dp_ord(StartOrder); 
    G = map(chkPoly, G, Vars); 
    G = map(dp_sort, G);
    GM = map(dp_hm, G); //Marking initial as Order1
    dp_ord(TargetOrder);
    G = map(dp_sort,G); 
    Tw= Tw2= Tini = Tlift = Tred = Step = 0;
    while(1){
T0 = time()[0];
        W = compute_last_w(G,GM, W, StartOrder, TargetOrder);
Tw += time()[0]-T0;
//print(W);

//print("W:", 0); print(W,0); //print(L[0]);
        if(1 == W){
            print("");
            print(["w",Tw,"ini",Tini,"lift",Tlift,"red",Tred]);
            break;
        }
print("step ",0); print(++Step, 2); print(["w",Tw,"ini",Tini,"lift",Tlift,"red",Tred]);
        

T0 = time()[0];
//print("Ini:", 0);
        N = length(G);
        IniG = newvect(N);
        for(I=0; I<N; I++){
            IniG[I] = initial_W(G[I], GM[I], W, StartOrder, TargetOrder);
        }
Tini += time()[0]-T0;

 T0 = time()[0];   
//print("nd_gr:", 0);
        H = nd_gr(vtol(IniG), Vars, 0, TargetOrder|dp=1);
//print("done",0);
        H = map(dp_ptozp, H);
        Hmark = map(dp_hm, H);
//print("lifting:", 0);
        Hd = lifting_mark(H, G, GM);
//print("done",0);
Tlift += time()[0]-T0;

T0 = time()[0];
//print("reducing:", 0);
    [G, GM] = interreduce_mark(Hd,Hmark);
//print("done");
Tred += time()[0]-T0;
    }

print("genwalk : Complete");
    return G;
}

print("genwalk loaded")$ 

S = generic_walk_r(G_s, Order_s, Order_t, Vars)$

G_tp = qsort(G_t)$//map(dp_dtop, G_t, Vars)$
S_p = map(dp_dtop, S, Vars)$
S_p = map(ptozp, S_p)$
S_p = qsort(S_p)$
cat(["S == G_t? ", gb_comp(G_tp, S_p)? "True":"False"])$

end$
