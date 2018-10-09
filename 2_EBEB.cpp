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
double py,k,k_A,a0,y0,t0,Ci;
if(flag==0){
  cout<< "please input the time step you want to store the data"<<endl;
  cin>>gap;
py=5;k=0.3;k_A=0.5;a0=1.0;y0=0;t0=0.0;}
else{
infile.open("restart.txt");
infile>>a0;
infile>>k;
infile>>k_A;
infile>>t0;
infile>>py;
infile>>y0;
infile>>Ci;
infile>>gap;
infile.close();
}
cout<<"please input your expected points in mapping"<<endl;
cin>>totalnum;
double h=0.0001,kesi=t0,E0;
int num=1,N=totalnum;
double Kessi[N],energy[N],Ax,w_posi,w_neg,dw_posi,dw_neg,dw2_posi,dw2_neg;
vector<double> y;
vector<double> pyta;
vector<double> allenergy;
vector<double> Xi;
Kessi[0]=t0;
y.push_back(y0);
pyta.push_back(py);
Xi.push_back(kesi);


Ax=a0*sin(kesi);
w_posi=k*y0*y0/2.0+k_A*y0*y0/2.0;
w_neg=k*y0*y0/2.0-k_A*y0*y0/2.0;
//Ci=sqrt(1+py*py)+w_posi;
Ci=8;
energy[0]=0.5*((1+py*py+Ax*Ax)/(Ci-w_posi)+w_neg);
allenergy.push_back(energy[0]);
double b[2]={0},x_next[2]={0},dfdx[2][2]={0},f[2]={0},diff[2]={0};
double y_n,py_n,b1,b2;
ofstream outfile;
outfile.open("parameters.txt");
outfile<<"a0 "<<a0<<endl;
outfile<<"k "<<k<<endl;
outfile<<"h "<<h<<endl;
outfile<<"gap "<<gap<<endl;
outfile<<"k_A "<<k_A<<endl;
outfile.close();
int gapp=0;
while(num<N){
  y_n=y0;
  py_n=py;
  Ax=a0*sin(kesi);
  w_posi=k*y0*y0/2.0+k_A*y0*y0/2.0;
  w_neg=k*y0*y0/2.0-k_A*y0*y0/2.0;
  dw_posi=k*y0+k_A*y0;
  dw_neg=k*y0-k_A*y0;
  b[0]=py-h/4*((1+py*py+Ax*Ax)/((Ci-w_posi)*(Ci-w_posi))*dw_posi+dw_neg);
  b[1]=y0+h/2.0*py/(Ci-w_posi);
  while(1){
    w_posi=k*y0*y0/2.0+k_A*y0*y0/2.0;
    w_neg=k*y0*y0/2.0-k_A*y0*y0/2.0;
    dw_posi=k*y0+k_A*y0;
    dw_neg=k*y0-k_A*y0;
    dw2_posi=k+k_A;
    dw2_neg=k-k_A;
    f[0]=py+h/4*((1+py*py+Ax*Ax)/((Ci-w_posi)*(Ci-w_posi))*dw_posi+dw_neg);
    f[1]=y0-h/2.0*py/(Ci-w_posi);
    dfdx[0][0]=h/4*(1+py*py+Ax*Ax)/((Ci-w_posi)*(Ci-w_posi))*dw2_posi+h/4*2.0*(1+py*py+Ax*Ax)/pow((Ci-w_posi),3)*dw_posi*dw_posi
               +h/4.0*dw2_neg;
    dfdx[0][1]=1.0+h/2*py/((Ci-w_posi)*(Ci-w_posi))*dw_posi;
    dfdx[1][0]=1.0-h/2*py/((Ci-w_posi)*(Ci-w_posi))*dw_posi;
    dfdx[1][1]=-h/2.0/(Ci-w_posi);
    b1=b[0]-f[0];
    b2=b[1]-f[1];
    diff[1]=(b1-b2*dfdx[0][0]/dfdx[1][0])/(dfdx[0][1]-dfdx[1][1]*dfdx[0][0]/dfdx[1][0]);
    diff[0]=(b2-dfdx[1][1]*diff[1])/dfdx[1][0];

    x_next[0]=y0+diff[0];
    x_next[1]=py+diff[1];
    if(fabs(diff[0])<err_rel&&fabs(diff[1])<err_rel&&fabs(f[0]-b[0])<err_tol&&fabs(f[1]-b[1])<err_tol) break;
    y0=x_next[0];
    py=x_next[1];
  }
  gapp+=1;

  if(gapp==gap){
 // y.push_back(y0);
 // pyta.push_back(py);
//  Xi.push_back(kesi);
//cout<<"test"<<endl;
Ax=a0*sin(kesi);
w_posi=k*y0*y0/2.0+k_A*y0*y0/2.0;
w_neg=k*y0*y0/2.0-k_A*y0*y0/2.0;
E0=0.5*((1+py*py+Ax*Ax)/(Ci-w_posi)+w_neg);
allenergy.push_back(E0);
  gapp=0;}

   if(py*py_n<0){
     Ax=a0*sin(kesi);
     w_posi=k*y0*y0/2.0+k_A*y0*y0/2.0;
     w_neg=k*y0*y0/2.0-k_A*y0*y0/2.0;
     energy[num]=0.5*((1+py*py+Ax*Ax)/(Ci-w_posi)+w_neg);

     Kessi[num]=kesi;

  cout << "The current number of points is "<<num<<endl;
    num+=1;
   }
  kesi=kesi+h;
//cout<<kesi<<endl;
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
//outfile.open("z.txt");
/*for(int i=0;i<y.size();i++){
  outfile<<y[i]<<endl;}
outfile.close();
outfile.open("delta.txt");
for(int i=0;i<y.size();i++){
  outfile<<pyta[i]<<endl;}
outfile.close();
outfile.open("allenergy.txt");
for(int j=0;j<y.size();j++){
  outfile<<allenergy[j]<<endl;}
outfile.close();
outfile.open("AllXi.txt");
for(int j=0;j<y.size();j++){
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
/*  outfile.open("z.txt",ios::app);
  for(int i=0;i<y.size()/gap-1;i++){
    outfile<<y[i*gap]<<endl;}
  outfile.close();
  outfile.open("delta.txt",ios::app);
  for(int i=0;i<y.size()/gap-1;i++){
    outfile<<pyta[i*gap]<<endl;}
  outfile.close();
  outfile.open("allenergy.txt",ios::app);
  for(int j=0;j<y.size()/gap-1;j++){
    outfile<<allenergy[j*gap]<<endl;}
  outfile.close();
  outfile.open("AllXi.txt",ios::app);
  for(int j=0;j<y.size()/gap-1;j++){
    outfile<<Xi[j*gap]<<endl;}
  outfile.close();*/

}

outfile.open("restart.txt");
outfile<<a0<<endl;
outfile<<k<<endl;
outfile<<k_A<<endl;
outfile<<kesi<<endl;
outfile<<py<<endl;
outfile<<y0<<endl;
outfile<<Ci<<endl;
outfile<<gap<<endl;
outfile.close();


}
