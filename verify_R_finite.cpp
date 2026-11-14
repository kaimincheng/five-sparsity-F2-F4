#include <cassert>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <map>
#include <utility>
using namespace std;
static const int L=326, LIM=1630;
uint8_t coeff[L+1];
int g[3][5], ng[3];
int toggle(int y){coeff[y]^=1;return y;}
int cubic(int t,int a,int b){
 memset(coeff,0,sizeof(coeff));
 ng[0]=ng[1]=ng[2]=0;
 int E[5]={0,1,b,a,t};
 for(int e:E)g[e%3][ng[e%3]++]=e/3;
 for(int r=0;r<3;r++){
  for(int i=0;i<ng[r];i++)for(int j=0;j<ng[r];j++) toggle(r+g[r][i]+2*g[r][j]);
 }
 for(int i=0;i<ng[0];i++)for(int j=0;j<ng[1];j++)for(int k=0;k<ng[2];k++)toggle(1+g[0][i]+g[1][j]+g[2][k]);
 int n=0;
 for(int i=0;i<=t;i++)n+=coeff[i];
 return n;
}
int f16mul(int a,int b){int z=0;for(int i=0;i<4;i++)for(int j=0;j<4;j++)if((a>>i&1)&&(b>>j&1))z^=(1<<(i+j));for(int i=6;i>=4;--i)if(z>>i&1)z^=0b10011<<(i-4);return z;}
int f16pow(int a,int n){int ans=1;while(n){if(n&1)ans=f16mul(ans,a); a=f16mul(a,a);n>>=1;}return ans;}
int zeta[5],cnt[LIM];
int fifth(int t,int a,int b){
 int E[5]={0,1,b,a,t};
 memset(cnt,0,sizeof(cnt));
 for(int u=0;u<5;u++)for(int v=0;v<5;v++)for(int w=0;w<5;w++)for(int x=0;x<5;x++)for(int y=0;y<5;y++){
  int ee=E[u]+E[v]+E[w]+E[x]+E[y];
  int e=(E[v]+2*E[w]+3*E[x]+4*E[y])%5;
  cnt[ee]^=zeta[e];
 }
 int nondiv=0,nonbinary=0,wt=0;
 for(int i=0;i<=5*t;i++){
  if(i%5 && cnt[i])nondiv++;
  if(i%5==0){if(cnt[i]!=0 && cnt[i]!=1)nonbinary++; if(cnt[i])wt++;}
 }
 if(nondiv||nonbinary){cerr<<"ERROR F16 descent "<<nondiv<<" "<<nonbinary<<endl;exit(1);}
 return wt;
}
int main(){
 int z=f16pow(2,3);
 for(int i=0;i<5;i++)zeta[i]=f16pow(z,i);
 assert(f16pow(z,5)==1 && z!=1);
 long long total=0, small=0, both=0;
 int min5=9999;
 map<int,int> hist;
 for(int t=15;t<=325;t+=2){
  for(int a=3;a<t;a++)for(int b=2;b<a;b++){
   total++;
   int wt=cubic(t,a,b);
   if(wt<=5){
    small++;
    int w5=fifth(t,a,b);
    min5=min(min5,w5);hist[w5]++;
    if(w5<=5){cerr<<"COUNTEREXAMPLE "<<t<<" "<<a<<" "<<b<<endl;both++;}
   }
  }
 }
 cout<<"Total candidates: "<<total<<"\nCubic survivors: "<<small<<"\nBoth small: "<<both<<"\nMin fifth norm weight: "<<min5<<"\n";
 cout<<"Fifth weight histogram:";
 for(auto [w,c]:hist)cout<<" "<<w<<":"<<c;
 cout<<"\n";
 if(total!=2821078 || small!=272 || both!=0 || min5!=7){cerr<<"C1 FAILED"<<endl; return 1;}
 cout<<"C1: PASS"<<endl;
}
