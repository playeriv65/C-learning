#include<bits/stdc++.h>
#include<iostream>
#include<cstdio>
#include<cmath>
using namespace std;
using ll = long double;
using ld = long double;
const int N = 1e5 + 5;
const ll inf = 1e10;
const ld eps = 1e-10;
int num, m;
ll suma, sumb;
ld pre[N], sub[N], ans;
ll lenp[N], lens[N];
struct node{
    ll l, r, len;
    ld mid;
    node(ll a, ll b){l = a,r = b, len = r-l,mid = (l+r)/(ld)2.0;}
}a[N], b[N];

int main(){
    // ios::sync_with_stdio(0),cin.tie(0);
    scanf("%d%d",&num,&m);
    for(int i = 1; i <= num; i ++){
        scanf("%Lf%Lf",&a[i].l,&a[i].r);
        a[i].len = a[i].r - a[i].l;
        suma += a[i].len;
    }
    if(!suma){
        for(int i = 1; i <= num; i ++){
            a[i].len = 1;
        }
        suma = num;
    }
    for(int i = 1; i <= m; i ++){
        scanf("%Lf%Lf",&b[i].l,&b[i].r);
        b[i].len = b[i].r - b[i].l;
        sumb += b[i].len;
    }
    if(!sumb){
        for(int i = 1; i <= m; i ++){
            b[i].len = 1;
        }
        sumb = m;
    }
    sort(a + 1, a + 1 + num, [](node x, node y){
        return x.l < y.l;
    });
    sort(b + 1, b + 1 + m, [](node x, node y){
        return x.l < y.l;
    });
    for(int i = 1; i <= m; i ++){
        lenp[i] = lenp[i-1] + b[i].len;
        pre[i] = pre[i-1] + (ld)(b[i].l + b[i].r) / 2 * b[i].len;
    }
    for(int i = m; i; i --){
        lens[i] = lens[i+1] + b[i].len;
        sub[i] = sub[i+1] + (ld)(b[i].l + b[i].r) / 2 * b[i].len; 
    }
    b[0].l = b[0].r = -inf;
    b[m+1].l = b[m+1].r = inf;
    for(int i = 1, l, r, mid, L, R; i <= num; i ++){
        // l = 0, r = m + 1;
        // while(l < r){
        //     mid = (l + r) / 2;
        //     if(b[mid].r < a[i].l) l = mid + 1;
        //     else r = mid;
        // }
        // L = l;

        while(b[L+1].r <= a[i].l){
            L++;
        }

        // l = 0, r = m+1;
        // while(l < r){
        //     mid = (l + r + 1) / 2;
        //     if(b[mid].l > a[i].r) r = mid - 1;
        //     else l = mid;
        // }
        // R = l;

        while(b[R].l < a[i].r){
            R++;
        }
        
        // cout << "TEST!! " << a[i].l << " " << a[i].r << " " << L << " " << R << '\n';
        if(L == m + 1){
            ans += (ld)(a[i].l + a[i].r) / 2 * a[i].len * lenp[L-1] - pre[L-1];
        }
        else if(R == 0){
            ans += sub[R+1] - (ld)(a[i].l + a[i].r) / 2 * a[i].len * lens[R+1];
        }
        else if(L > R){
            ans += (ld)(a[i].l + a[i].r) / 2 * a[i].len * lenp[L-1] - pre[L-1];
            ans += sub[R+1] - (ld)(a[i].l + a[i].r) / 2 * a[i].len * lens[R+1];
        }
        else{
            ans += (ld)(a[i].l + a[i].r) / 2 * a[i].len * lenp[L-1] - pre[L-1];
            ans += sub[R+1] - (ld)(a[i].l + a[i].r) / 2 * a[i].len * lens[R+1];
            for(int j = L; j <= R; j ++){
                ld l = max(a[i].l, b[j].l);
                ld r = min(a[i].r, b[j].r);
                if(a[i].l == a[i].r && a[i].len == 1){
                    ans += (ld)max((ld)(b[j].r - b[j].l) / 2, abs(a[i].l - (ld)(b[j].l + b[j].r) / 2));
                }
                else if(b[j].l == b[j].r && b[j].len == 1){
                    ans += (ld)max((ld)(a[i].r - a[i].l) / 2, abs(b[j].l - (ld)(a[i].l + a[i].r) / 2));
                }
                else{
                    ans += (ld)(r - l) * (r - l) * (r - l) / 3;
                    ll qa[4] = {a[i].l, l, r, a[i].r};
                    ll qb[4] = {b[j].l, l, r, b[j].r};
                    for(int k = 0; k < 3; k ++){
                        for(int p = 0; p < 3; p ++){
                            // cout << qa[k] << ' ' << qa[k + 1] << " " << qb[p] << " " << qb[p+1] << '\n';
                            if(qa[k] == qb[p] && qa[k+1]==qb[k+1])
                                continue;
                            ans += (ld)(qa[k+1]-qa[k]) * (qb[p+1]-qb[p]) * abs((qa[k+1]+qa[k])/2 - (qb[p+1]+qb[p])/2);
                        }
                    }
                }
            }
        }
    }
    //cout << suma << " " << sumb << '\n';
    ans /= suma * sumb;
    printf("%0.11Lf\n",ans);
    return 0;
}
/*
1 1
0 2
1 3
*/