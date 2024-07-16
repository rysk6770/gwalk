load("rysk.c")$
load("noro_gw.rr")$

N = 7$
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

end$