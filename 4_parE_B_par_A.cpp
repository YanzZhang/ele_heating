#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#define err_rel pow(10,-12)
#define err_tol pow(10,-12)
using namespace std;




int main()
{
  int flag,totalnum,gap;;
  cout<<"restar (1) or not (0)"<<endl;
  cin>>flag;
  ifstream infile;
  double del,k,k_A,a0,z0,t0,Cy=0;
  if(flag==0){

  cout<< "please input the time step you want to store the data"<<endl;
  cin>>gap;
  del=1;k=0.3,k_A=0.5;a0=1;z0=0;t0=0.0;}
  else{
  infile.open("restart.txt");
  infile>>a0;
  infile>>k;
  infile>>k_A;
  infile>>t0;
  infile>>del;
  infile>>z0;
  infile>>gap;
  infile.close();
  }
  cout<<"please input your expected points in mapping"<<endl;
  cin>>totalnum;

double h=0.0001,kesi=t0,E0;
int num=1,N=totalnum;
double Kessi[N],energy[N],Ax,Ab,Uz,dAb,dUz,d2Ab,d2Uz;
vector<double> z;
vector<double> delta;
vector<double> allenergy;
vector<double> Xi;
Kessi[0]=t0;
z.push_back(z0);
delta.push_back(del);
Xi.push_back(kesi);
Ax=a0*sin(kesi);
Ab=k_A*z0*z0/2.0;
Uz=k*z0*z0/2.0;
energy[0]=0.5*(1+(Cy+Ab)*(Cy+Ab)+Ax*Ax)/del+del/2.0+Uz;
allenergy.push_back(energy[0]);
double b[2]={0},x_next[2]={0},dfdx[2][2]={0},f[2]={0},diff[2]={0};
double z_n,del_n,b1,b2;
ofstream outfile;

int gapp=0;

outfile.open("parameters.txt");
outfile<<"a0 "<<a0<<endl;
outfile<<"k "<<k<<endl;
outfile<<"h "<<h<<endl;
outfile<<"gap "<<gap<<endl;
outfile<<"k_A "<<k_A<<endl;
outfile<<"Ci "<<Cy<<endl;

outfile.close();

while(num<N){
  z_n=z0;
  del_n=del;
  Ax=a0*sin(kesi);
  Ab=k_A*z0*z0/2.0;
  Uz=k*z0*z0/2.0;
  dAb=k_A*z0;
  dUz=k*z0;
  b[0]=del+h/2.0*(dUz+(Cy+Ab)*dAb/del);
  b[1]=z0-h/2.0+h/4.0/(del*del)*(1.0+(Cy+Ab)*(Cy+Ab)+Ax*Ax);
  while(1){
    Ax=a0*sin(kesi);
    Ab=k_A*z0*z0/2.0;
    Uz=k*z0*z0/2.0;
    dAb=k_A*z0;
    dUz=k*z0;
    d2Ab=k_A;
    d2Uz=k;
    f[0]=del-h/2.0*(dUz+(Cy+Ab)*dAb/del);
    f[1]=z0-h/4.0/(del*del)*(1.0+(Cy+Ab)*(Cy+Ab)+Ax*Ax);
    dfdx[0][0]=-h/2.0*(d2Uz+(Cy+Ab)*d2Ab/del+dAb*dAb/del);
    dfdx[0][1]=1.0+h/2.0*(Cy+Ab)*dAb/(del*del);
    dfdx[1][0]=1.0-h/2.0/(del*del)*(Cy+Ab)*dAb;
    dfdx[1][1]=h/2.0/pow(del,3)*(1+(Cy+Ab)*(Cy+Ab)+Ax*Ax);
    b1=b[0]-f[0];
    b2=b[1]-f[1];
    diff[1]=(b1-b2*dfdx[0][0]/dfdx[1][0])/(dfdx[0][1]-dfdx[1][1]*dfdx[0][0]/dfdx[1][0]);
    diff[0]=(b2-dfdx[1][1]*diff[1])/dfdx[1][0];

    x_next[0]=z0+diff[0];
    x_next[1]=del+diff[1];
    if(fabs(diff[0])<err_rel&&fabs(diff[1])<err_rel&&fabs(f[0]-b[0])<err_tol&&fabs(f[1]-b[1])<err_tol) break;
    z0=x_next[0];
    del=x_next[1];
  }
  /*gapp+=1;
  if(gapp==gap){
  z.push_back(z0);
  delta.push_back(del);
    Xi.push_back(kesi);
    Ax=a0*sin(kesi);
    Ab=k_A*z0*z0/2.0;
    Uz=k*z0*z0/2.0;
  E0=0.5*(1+(Cy+Ab)*(Cy+Ab)+Ax*Ax)/del+del/2.0+Uz;
  allenergy.push_back(E0);
  gapp=0;
}*/

   if(z0*z_n<0){
     Ax=a0*sin(kesi);
     Ab=k_A*z0*z0/2.0;
     Uz=k*z0*z0/2.0;
     energy[num]=0.5*(1+(Cy+Ab)*(Cy+Ab)+Ax*Ax)/del+del/2.0+Uz;

     Kessi[num]=kesi;

  cout << "The current number of points is "<<num<<endl;
    num+=1;
   }
  kesi=kesi+h;
}

if(flag==0){
outfile.open("Xi.txt");
for(int i=0;i<N;i++){
outfile<<Kessi[i]<<endl;}
outfile.close();
//cout<<"test"<<endl;
outfile.open("energy.txt");

for(int i=0;i<N;i++){
outfile<<energy[i]<<endl;}
outfile.close();
/*outfile.open("z.txt");
for(int i=0;i<z.size();i++){
  outfile<<z[i]<<endl;}
outfile.close();
outfile.open("delta.txt");
for(int i=0;i<z.size();i++){
  outfile<<delta[i]<<endl;}
outfile.close();
outfile.open("allenergy.txt");
for(int j=0;j<z.size();j++){
  outfile<<allenergy[j]<<endl;}
outfile.close();
outfile.open("AllXi.txt");
for(int j=0;j<z.size();j++){
  outfile<<Xi[j]<<endl;}
outfile.close();*/
}
else{
  outfile.open("Xi.txt",ios::app);
  for(int i=0;i<N;i++){
  outfile<<Kessi[i]<<endl;}
  outfile.close();
  //cout<<"test"<<endl;
  outfile.open("energy.txt",ios::app);

  for(int i=0;i<N;i++){
  outfile<<energy[i]<<endl;}
  outfile.close();
  outfile.open("z.txt",ios::app);
  for(int i=0;i<z.size()/gap-1;i++){
    outfile<<z[i*gap]<<endl;}
  outfile.close();
  outfile.open("delta.txt",ios::app);
  for(int i=0;i<z.size()/gap-1;i++){
    outfile<<delta[i*gap]<<endl;}
  outfile.close();
  outfile.open("allenergy.txt",ios::app);
  for(int j=0;j<z.size()/gap-1;j++){
    outfile<<allenergy[j*gap]<<endl;}
  outfile.close();
  outfile.open("AllXi.txt",ios::app);
  for(int j=0;j<z.size()/gap-1;j++){
    outfile<<Xi[j*gap]<<endl;}
  outfile.close();

}
outfile.open("restart.txt");
outfile<<a0<<endl;
outfile<<k<<endl;
outfile<<k_A<<endl;
outfile<<kesi<<endl;
outfile<<del<<endl;
outfile<<z0<<endl;
outfile<<gap<<endl;

outfile.close();



}
