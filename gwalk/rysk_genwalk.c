load("rysk.c")$
load("noro_gw.rr")$

def vec_to_mat(V){
    N = length(V);
    Mat = newmat(N, 1);
    for(I = 0; I<N; I++){
        Mat[I][0] = V[I];
    }
    return Mat;
}

def mat_transpose(Mat){
    [N, M] = size(Mat);
    NewMat = newmat(M,N);
    for(I=0; I<N; I++){
        for(J=0; J<M; J++){
            NewMat[J][I] = Mat[I][J];
        }
    }
    return NewMat;
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
    compare vectors V1 and V2 by weight(or matrix) M.
    if V1 < V2 then return 1
    if V1 > V2 then return -1
    if V1 == V2 then return 0
*/
def compare(V1, V2, M){
    // when M is vector.
    if(5 == type(M)) V = innerProduct(M, V1-V2);
    // when M is matrix.
    if(6 == type(M)){
        N = size(M)[0];
        V = M*(V1 - V2);
        for(I = 0; I < N; I++){
            if(V[I] != 0){
                V = V[I];
                break;
            }
        }
    }
    if(V > 0) return -1;
    else if(V < 0) return 1;
    else return 0;
}
/*
    return degree vectors of polynomial F
*/
def supp(Poly){
    Result = [];
    while(Poly){
        Result = cons(dp_etov(Poly), Result);
        Poly = dp_rest(Poly);
    }
    return reverse(Result);
}
/*
    if V is in C12 return 1
    else return 0
*/
def in_C12(V, Order1, Order2){
    T1 = compare(0, V, Order1);
    T2 = compare(0, V, Order2);
    if((T1 == 1) && (T2 == -1)) return 1;
    else return 0;
}

/*
    for facet preorder < defined by matrice M1 and M2 = O(T)
    if u < v return 1
    else return 0
*/
def facet_preorder(U, V, M1, M2){
    if ( U == 0 || V == 1) return 1;
	if ( U == 1 || V == 0) return 0;
	N = size(V)[0];
    
	for ( I = 0; I < N; I++ ) {
		Left = innerProduct(M2[I],U)*V;
		Right = innerProduct(M2[I],V)*U;
		R = compare(Left,Right,M1);
		if ( R < 0 ) return 1;
		else if ( R > 0 ) return 0;
	}
	return 1;
}

def compute_last_w(G, W, Order1, Order2){
    DG = []; N = length(G);
    //if(W == 0) W=newvect(size(Order1)[0]);
    for(I = 0; I < N; I++){
        DG = append(compute_df(G[I]), DG);
    }
    DG = set_of(DG); V = []; N = length(DG);
    for(I = 0; I < N; I++){
        if( in_C12(DG[I], Order1, Order2) && facet_preorder(W, DG[I], Order1, Order2)){
            V = cons(DG[I], V);
        }
    }
    
    if(V == []) return 1;
    for(T = V; T!=[]; T=cdr(T)){
        S = car(T);
        for(U = V; U!=[]; U=cdr(U)){
            if(!facet_preorder(S, car(U), Order1, Order2)){
                break;
            }
        }
        if(U == []){
            return S;
        }
    }
    print("c");
    error("compute_last_w");
}

/*
    GM is marked term of G
    M1 and M2 define facet preorder
*/
def initial_W(G, W, M1, M2){
    GM = dp_hm(G);
    G = dp_rest(G);
    U = dp_etov(GM);
    S = [];
    for(; G; G = dp_rest(G)){
        V = dp_etov(G);
        T = U - V;
        if( facet_preorder(T, W, M1, M2) && facet_preorder(W, T, M1, M2)){
            GM += dp_hm(G);
        }
    }
    return GM;
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

def interreduce_mark(G, GM){
    
    G = vtol(G);
    N = length(G);
    Result = [];
    GM_L = [];
    for(I = 0; I < N; I++){
        P = car(G); PM = car(GM);
        G = cdr(G); GM = cdr(GM);
        Q = reducing_mark(P, append(Result, G),append(GM_L, GM));
        Result = cons(Q, Result);   
        GM_L = cons(PM, GM_L);     
    }
    return Result;
}


def generic_walk(G, StartOrder, TargetOrder, Vars){
    dp_ord(StartOrder); 
    G = map(chkPoly, G, Vars); 
    G = map(dp_sort, G);
    G_mark = map(dp_hm, G); //Marking initial as Order1
   // GM = map(dp_hm, G);
    dp_ord(TargetOrder);
    T1 = T0 = Tw = Step = 0;
    while(1){
    print("step : ",0); print(++Step);
        W = compute_last_w(G, W, StartOrder, TargetOrder);
    print(["W : ", W]);
        if(1 == W) break;
    
        N = length(G);
        IniG = newvect(N);
        for(I=0; I<N; I++){
            IniG[I] = initial_W(G[I], W, StartOrder, TargetOrder);
        }
    //print(IniG);
    //print(N);
        H = map(chkPoly, nd_gr(vtol(IniG), Vars, 0, TargetOrder), Vars);
        H_mark = map(dp_hm, H);
    //print(H); //Marking initial term as Order2
        //G = map(dp_sort, G);
        H_r = map(reducing_mark, H, G, G_mark);
        H_d = [];
        N = length(H);
        for(I = 0; I<N; I++){
            F = H[I] - H_r[I];
            //print(["I", I, H[I] == H_r[I] ]);
            //print(F);
            if(F == 0){
                H_d = cons(H[I], H_d);
            }else{
                H_d = cons(F, H_d);
            }
        }
        //print(H_d);
        H_d = reverse(H_d);
        //H_d = 
        //print(H_d);
    print(N);
        dp_ord(TargetOrder);
        H_d = map(dp_sort, H_d);
        G = interreduce_mark(H_d,H_mark);
        G_mark = H_mark;
        //print(G);
        G = map(dp_sort, G);
    }
    print("genwalk : Complete");
    return G;
}




N = 5$
Cyclic = 1$
if(Cyclic) B = cyclic(N)$
else B = katsura(N)$
Vars = vars(B)$
N = length(Vars)$
Order_s = create_order_mat(N, 0)$
Order_t = create_order_mat(N, 2)$

if(1){
    G_s = nd_gr(B, Vars, 0, 0 | dp = 1)$
    G_t = noro_gw.generic_walk(map(dp_dtop,G_s, Vars), Vars, Order_s, Order_t)$
}else{
    G_s = dp_f4_main(B, Vars, Order_s)$
    G_t = noro_gw.generic_walk(G_s, Vars, Order_s, Order_t)$
}


print("start genwalk")$ cputime(1)$
//S= generic_walk(G_s, Order_s, Order_t, Vars)$
cputime(0)$
G_tp = qsort(G_t)$//map(dp_dtop, G_t, Vars)$
S_p = map(dp_dtop, S, Vars)$
S_p = map(ptozp, S_p)$
S_p = qsort(S_p)$
cat(["S == G_t? ", gb_comp(G_tp, S_p)? "True":"False"])$
end$
