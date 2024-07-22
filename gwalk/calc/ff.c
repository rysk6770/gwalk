#include "calc.h"

void print_ff(Coef a){
    printf("%lld", a.f);
}

int zero_ff(Coef t){
    return t.f == 0;
}

Coef one_ff(){
    Coef c;
    c.f = 1;
    return c;
}

Coef add_ff(Coef a, Coef b){
    LONG c,p;
    Coef r;

    p = CurrentRing->chr;
    c = (a.f + b.f)%p;
    if(c < 0) c += p;
    r.f = c;
    return r;
}

Coef neg_ff(Coef a){
    LONG c,p;
    Coef r;

    p = CurrentRing->chr;
    c = (p-a.f)%p;
    if(c < 0) c += p;
    r.f = c;
    return r;
}

Coef sub_ff(Coef a, Coef b){
    LONG c,p;
    Coef r;

    p = CurrentRing->chr;
    c = (a.f - b.f)%p;
    if(c < 0) c += p;
    r.f = c;
    return r;
}

Coef mul_ff(Coef a, Coef b){
    LONG c,p;
    Coef r;

    if(a.f == 0 || b.f == 0){
        c = 0;
    }else{
        p = CurrentRing->chr;
        c = (a.f * b.f)%p;
        if(c < 0) c += p;
    }
    r.f = c;
    return r;
}

// c-a*b
Coef mulsub_ff(Coef c,Coef a,Coef b)
{
  LONG d,p;
  Coef r;

  if ( a.f == 0 || b.f == 0 ) return c;
  p = CurrentRing->chr;
  d = (c.f-a.f*b.f) % p;
  if ( d < 0 ) d += p;
  r.f = d;
  return r;
}

Coef inv_ff(Coef a){
    LONG f1,f2,a1,a2,q,r,p;
    Coef u;
    p=CurrentRing ->chr;
    f1=a.f; f2=p; a1=1; a2=0;
    while(1){
        q=f1/f2; r=f1-f2*q;
        f1=f2; f2=r;
        if(f2 == 0) 
            break;
        r=(a2*q)%p; r=a1-r;
        if(r < 0) r+=p;
        a1=a2; a2=r;
    } 
    u.f = a2;
    return u;
}

Coef div_ff(Coef a, Coef b){
    LONG c,p;
    Coef inv, r;
    if(b.f == 0)
        error("division by 0");
    if(a.f == 0) return a;
    inv = inv_ff(b);
    p = CurrentRing->chr;
    c = (a.f * inv.f)%p;
    if(c < 0) c += p;
    r.f = c;
    return r;
}

Coef ntoc_ff(char *n){
    LONG a,p;
    Coef r;
    p = CurrentRing->chr;
    a = strtol(n,0,10)%p;
    if(a<0) a+=p;
    r.f=a;
    return r;
}

