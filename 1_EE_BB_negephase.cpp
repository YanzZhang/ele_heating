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
double py,k,k_A,a1,y0,t0,Ci,alpha,px_tilt,a2;
if(flag==0){
  //cout<< "please input the time step you want to store the data"<<endl;
//  cin>>gap;
//gap=5000;
//py=8;k=0.3;k_A=0.5;a1=1.0,a2=1.0;y0=0.0;t0=0.0;alpha=1.1;px_tilt=0;
//Ci=8.0;//Ci=gamma-alpha*Pz for initial y0=0 which is less than sqrt(1+P_perp^2)
infile.open("InitialBC12.txt");
infile>>a1;
infile>>a2;
infile>>alpha;
infile>>k;
infile>>k_A;
infile>>t0;
infile>>py;
infile>>y0;
infile>>Ci;
infile>>px_tilt;
infile.close();
}
else{
infile.open("restart.txt");
infile>>a1;
infile>>a2;
infile>>k;
infile>>k_A;
infile>>t0;
infile>>py;
infile>>y0;
infile>>Ci;
infile>>gap;
infile>>alpha;
infile>>px_tilt;
infile.close();
}
gap=5000;
cout<<"please input your expected points in mapping"<<endl;
cin>>totalnum;
double h=0.0001,kesi=t0,E0;
double AB,U,dAB,dU,d2AB,d2U,w_AB,w_U,dw_AB,dw_U,dw2_AB,dw2_U,p_perp;
int num=1,N=totalnum;
double Kessi[N],energy[N],Ay,Ax;
vector<double> y;
vector<double> pyta;
vector<double> allenergy;
vector<double> Xi;
Kessi[0]=t0;
y.push_back(y0);
pyta.push_back(py);
Xi.push_back(kesi);
Ax=a2*sin(kesi);
Ay=a1*sin(kesi);
///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
U=k*y0*y0/2.0;
AB=k_A*y0*y0/2.0;
///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
w_AB=U+alpha*AB-Ci;
w_U=alpha*U+AB-Ci;
p_perp=1+(py+Ay)*(py+Ay)+(px_tilt+Ax)*(px_tilt+Ax);
energy[0]=alpha/(alpha*alpha-1)*(w_U+sqrt(w_AB*w_AB+(alpha*alpha-1)*p_perp));
allenergy.push_back(energy[0]);
double b[2]={0},x_next[2]={0},dfdx[2][2]={0},f[2]={0},diff[2]={0};
double y_n,py_n,b1,b2;
ofstream outfile;
outfile.open("parameters.txt");
outfile<<"a1 "<<a1<<endl;
outfile<<"a2 "<<a2<<endl;
outfile<<"k "<<k<<endl;
outfile<<"h "<<h<<endl;
outfile<<"gap "<<gap<<endl;
outfile<<"k_A "<<k_A<<endl;
outfile<<"Ci "<<Ci<<endl;
outfile<<"alpha "<<alpha<<endl;
outfile.close();

int gapp=0;
while(num<N){
  y_n=y0;
  py_n=py;
  Ay=a1*sin(kesi);
  Ax=a2*sin(kesi);

  ///////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  U=k*y0*y0/2.0;
  AB=k_A*y0*y0/2.0;
  dU=k*y0;
  dAB=k_A*y0;
  ///////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  w_AB=U+alpha*AB-Ci;
  w_U=alpha*U+AB-Ci;
  dw_AB=dU+alpha*dAB;
  dw_U=alpha*dU+dAB;
  p_perp=1+(py+Ay)*(py+Ay)+(px_tilt+Ax)*(px_tilt+Ax);
  b[0]=py-h/2*alpha/(alpha*alpha-1)*(dw_U+w_AB/sqrt(w_AB*w_AB+(alpha*alpha-1)*p_perp)*dw_AB);
  b[1]=y0+h/2*alpha/sqrt(w_AB*w_AB+(alpha*alpha-1)*p_perp)*(py+Ay);
  while(1){
  /*  if(w_AB>0) { cout << "w_AB "<<w_AB<<endl;
      return 1;
    }*/
    ///////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    p_perp=1+(py+Ay)*(py+Ay)+(px_tilt+Ax)*(px_tilt+Ax);
    U=k*y0*y0/2.0;
    AB=k_A*y0*y0/2.0;
    dU=k*y0;
    dAB=k_A*y0;
    d2U=k;
    d2AB=k_A;
    ///////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    w_AB=U+alpha*AB-Ci;
    dw_AB=dU+alpha*dAB;
    dw2_AB=d2U+alpha*d2AB;
    w_U=alpha*U+AB-Ci;
    dw_U=alpha*dU+dAB;
    dw2_U=alpha*d2U+d2AB;
    f[0]=py+h/2*alpha/(alpha*alpha-1)*(dw_U+w_AB/sqrt(w_AB*w_AB+(alpha*alpha-1)*p_perp)*dw_AB);
    f[1]=y0-h/2*alpha/sqrt(w_AB*w_AB+(alpha*alpha-1)*p_perp)*(py+Ay);

    dfdx[0][0]=h/2*alpha/(alpha*alpha-1)*(dw2_U+w_AB/sqrt(w_AB*w_AB+(alpha*alpha-1)*p_perp)*dw2_AB
               +(alpha*alpha-1)*p_perp*dw_AB*dw_AB/pow(w_AB*w_AB+(alpha*alpha-1)*p_perp,1.5));

    dfdx[0][1]=1.0-h/2*alpha/pow(w_AB*w_AB+(alpha*alpha-1)*p_perp,1.5)*w_AB*dw_AB*(py+Ay);
    dfdx[1][0]=1.0+h/2*alpha/pow(w_AB*w_AB+(alpha*alpha-1)*p_perp,1.5)*w_AB*dw_AB*(py+Ay);

    dfdx[1][1]=-h/2.0*alpha*(w_AB*w_AB+(alpha*alpha-1)*((px_tilt+Ax)*(px_tilt+Ax)+1))/pow(w_AB*w_AB+(alpha*alpha-1)*p_perp,1.5);
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
 // gapp+=1;

 /* if(gapp=gap){
  y.push_back(y0);
  pyta.push_back(py);
  Xi.push_back(kesi);
//cout<<"test"<<endl;
Ay=a1*sin(kesi);
Ax=a2*sin(kesi);

///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
U=k*y0*y0/2.0;
AB=k_A*y0*y0/2.0;
///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
w_AB=U+alpha*AB-Ci;
w_U=alpha*U+AB-Ci;
p_perp=1+(py+Ay)*(py+Ay)+(px_tilt+Ax)*(px_tilt+Ax);
E0=alpha/(alpha*alpha-1)*(w_U+sqrt(w_AB*w_AB+(alpha*alpha-1)*p_perp));
allenergy.push_back(E0);
  gapp=0;*/

   if(py*py_n<0){
     Ay=a1*sin(kesi);
     Ax=a2*sin(kesi);

     ///////////////////////////////////////////////////////////////////////////////////////////
     //////////////////////////////////////////////////////////////////////////////////////////
     U=k*y0*y0/2.0;
     AB=k_A*y0*y0/2.0;
     ///////////////////////////////////////////////////////////////////////////////////////////
     //////////////////////////////////////////////////////////////////////////////////////////
     w_AB=U+alpha*AB-Ci;
     w_U=alpha*U+AB-Ci;
     p_perp=1+(py+Ay)*(py+Ay)+(px_tilt+Ax)*(px_tilt+Ax);
     energy[num]=alpha/(alpha*alpha-1)*(w_U+sqrt(w_AB*w_AB+(alpha*alpha-1)*p_perp));
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
/*outfile.open("z.txt");
for(int i=0;i<y.size()/gap-1;i++){
  outfile<<y[i*gap]<<endl;}
outfile.close();
outfile.open("delta.txt");
for(int i=0;i<y.size()/gap-1;i++){
  outfile<<pyta[i*gap]<<endl;}
outfile.close();
outfile.open("allenergy.txt");
for(int j=0;j<y.size()/gap-1;j++){
  outfile<<allenergy[j*gap]<<endl;}
outfile.close();
outfile.open("AllXi.txt");
for(int j=0;j<y.size()/gap-1;j++){
  outfile<<Xi[j*gap]<<endl;}
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
 /* outfile.open("z.txt",ios::app);
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
outfile<<a1<<endl;
outfile<<a2<<endl;
outfile<<k<<endl;
outfile<<k_A<<endl;
outfile<<kesi<<endl;
outfile<<py<<endl;
outfile<<y0<<endl;
outfile<<Ci<<endl;
outfile<<gap<<endl;
outfile<<alpha<<endl;
outfile<<px_tilt<<endl;
outfile.close();


}
