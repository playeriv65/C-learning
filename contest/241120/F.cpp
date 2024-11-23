#include<bits/stdc++.h>
using namespace std;
#define ll long long
const ll N = 1e6+10;

ll T,n,m,sn,sm,ans[N];
bool tr;

ll sum(ll a){return a*(a+1)/2;}

int main(){
    cin>>T;
    while(T--){
        tr = 0;
        scanf("%lld%lld",&n,&m);
        sn = sum(n), sm = sum(m);
        if(sn%2 && sm%2){
            printf("No\n");
            continue;
        }
        printf("Yes\n");

        if(sn%2){
            tr = 1;
            swap(n,m);
            swap(sn,sm);
        }

        ll len = 0,tot = 0;
        ans[0] = 0;
        for(int i=1;i<=n;i++){
            if(tot+sum(len+1)+(n-i) > sn/2){
                tot += sum(len);
                len = 1;
                ans[i] = ans[i-1]^1;
            }
            else{
                len ++;
                ans[i] = ans[i-1];
            }

            // cout<<tot<<" "<<len<<" "<<endl;
        }

        // cout<<n<<sn<<endl;
        if(!tr){
            for(int i=1;i<=n;i++){
                for(int j=1;j<=m;j++){
                    printf("%d ", ans[i]);
                }
                printf("\n");
            }
        }
        else{
            for(int i=1;i<=m;i++){
                for(int j=1;j<=n;j++){
                    printf("%d ", ans[j]);
                }
                printf("\n");
            }
        }
    }
    return 0;
}