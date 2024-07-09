load("rysk.c")$
load("noro_gw.rr")$
/*
    G is groebner basis respect to StartOrder.
    G must be Dpoly(분산).
*/
def grwalk(G, StartOrder, TargetOrder, Vars){
    dp_ord(StartOrder); // Set Order as StartOrder.
    G = map(chkPoly, G, Vars); // if G is not Dpoly then convert to Dpoly.

    Mat_old = StartOrder; 
    G_old = G;
    Weight_new = StartOrder[0];
    Weight_target = TargetOrder[0];
    Mat_new = new_order(Weight_new, TargetOrder); 

    T1 = T0 = Tw = Step = 0;
    while(1){
    print("step : ",0); print(++Step);
        Ini_new = map(initialForm, G_old, Mat_new[0]); // Initial form of Weight_new

        G_Ini = map(chkPoly, nd_gr(vtol(Ini_new), Vars, 0, Mat_new), Vars);
    print("nd_gr: done    ");

    print("lifting     : ", 0); T0 = time()[0]$
        G_new = lifting(Ini_new, G_Ini, G_old, Mat_old);
    print("done. ", 0); print("Tlift : ", 0); print(time()[0] - T0);
    
    print("interreduce : ", 0); T0 = time()[0]$
        G_new = interreduce2(G_new, Mat_new);
    print("done. ", 0); print(" Tred : ", 0); print(time()[0] - T0);
        T = compute_last_t(G_new, Weight_new, Weight_target);
        if(Weight_new == Weight_target){
            break;
        }else{
            Mat_old = Mat_new;
            G_old = G_new;
            Weight_new += T*(Weight_target - Weight_new);
            Weight_new = conv2int(Weight_new); 

            Mat_new = new_order(Weight_new, TargetOrder);

        }

    }
    print("grwalk : Complete");
    return G_new;
}


N = 6$
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
    G_s = nd_f4_main(B, Vars, Order_s)$
    G_t = noro_gw.generic_walk(G_s, Vars, Order_s, Order_t)$
}
// 
//G_t = nd_gr(G_s, Vars, 0, 2 | dp = 1)$

print("start grwalk")$ cputime(1)$
S= grwalk(G_s, Order_s, Order_t, Vars)$
cputime(0)$
G_tp = qsort(G_t)$//map(dp_dtop, G_t, Vars)$
S_p = map(dp_dtop, S, Vars)$
S_p = map(ptozp, S_p)$
S_p = qsort(S_p)$
cat(["S == G_t? ", gb_comp(G_tp, S_p)? "True":"False"])$
end$
