//程序使用备注
/*
此程序通过C/C++编译实现单像空间后方交会的计算机程序，程序编译环境为Win10+vs2010控制台，
参考题目为摄影测量学第二版（张剑清等著）P44页第27题，
程序输入采用txt文件输入的形式，命名为data.txt，存放在D盘中，文件中的数据格式为：
  比例尺m，内方位元素x0,y0,f，控制点的个数
  各控制点的影像坐标和地面坐标(x,y,X,Y,Z)
此例中的数据文件内容为：
50000,0,0,153.24,4
-86.15,-68.99,36589.41,25273.32,2195.17
-53.40,82.21,37631.08,31324.51,728.69
-14.78,-76.63,39100.97,24934.98,2386.50
10.46,64.43,40426.54,30319.81,757.31
程序输出采用txt文件输出的形式，命名为Space_Resection.Output.txt，存放在D盘中
*/

#include "stdafx.h"
#include<iostream>
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <windows.h>
using namespace std;

//定义点对的影像坐标和地面坐标的数据结构
struct PtoP {
	double x;
	double y;
	double X;
	double Y;
	double Z;
};

//确定未知数的初值
void set_initial_unknownvalue(PtoP *pp,int n,double *m,double *f,double *Xs,double *Ys,double *Zs,double *phi,double *omega,double *kappa){
	double sum_Xs=0,sum_Ys=0;
	for(int i=0;i<n;i++){
		sum_Xs+= pp[i].X;
		sum_Ys+= pp[i].Y;
	}
	*Xs=sum_Xs/n;
	*Ys=sum_Ys/n;
	*Zs=(*m)*(*f);
	*phi=0;
	*omega=0;
	*kappa=0;
}

//计算旋转矩阵
void RotationMatrix(double phi_temp,double omega_temp,double kappa_temp,double *matrix){
	matrix[0]=cos(phi_temp)*cos(kappa_temp)-sin(phi_temp)*sin(omega_temp)*sin(kappa_temp);
	matrix[1]=-cos(phi_temp)*sin(kappa_temp)-sin(phi_temp)*sin(omega_temp)*cos(kappa_temp);
	matrix[2]=-sin(phi_temp)*cos(omega_temp);

	matrix[3]=cos(omega_temp)*sin(kappa_temp);
	matrix[4]=cos(omega_temp)*cos(kappa_temp);
	matrix[5]=-sin(omega_temp);

	matrix[6]=sin(phi_temp)*cos(kappa_temp) + cos(phi_temp)*sin(omega_temp)*sin(kappa_temp);
	matrix[7]=-sin(phi_temp)*sin(kappa_temp) + cos(phi_temp)*sin(omega_temp)*cos(kappa_temp);
	matrix[8]=cos(phi_temp)*cos(omega_temp);
}

//计算一对点的误差方程V=AX-L的参数和常数项
void ErrorEquation(PtoP *pp,double R[],double Ai[12],double Li[2],int t,double Xs,double Ys,double Zs,double phi,double omega,double kappa,double f){
	//辅助参数X、Y、Z（上横杠）
	double X_up,Y_up,Z_up;
	X_up= R[0]*(pp[t].X-Xs)+R[3]*(pp[t].Y-Ys)+R[6]*(pp[t].Z-Zs);
	Y_up= R[1]*(pp[t].X-Xs)+R[4]*(pp[t].Y-Ys)+R[7]*(pp[t].Z-Zs);
	Z_up= R[2]*(pp[t].X-Xs)+R[5]*(pp[t].Y-Ys)+R[8]*(pp[t].Z-Zs);

	//计算x部分的偏导数a11至a16:
	Ai[0+0]=(R[0]*f+R[2]*pp[t].x)/Z_up;
	Ai[0+1]=(R[3]*f+R[5]*pp[t].x)/Z_up;
	Ai[0+2]=(R[6]*f+R[8]*pp[t].x)/Z_up;
	Ai[0+3]=pp[t].y*sin(omega)-(pp[t].x/f*(pp[t].x*cos(kappa)-pp[t].y*sin(kappa))+f*cos(kappa))*cos(omega);
	Ai[0+4]=-f*sin(kappa)-pp[t].x/f*(pp[t].x*sin(kappa)+pp[t].y*cos(kappa));
	Ai[0+5]=pp[t].y;

	//计算y部分的偏导数a21至a26:
	Ai[6+0]=(R[1]*f+R[2]*pp[t].y)/Z_up;
	Ai[6+1]=(R[4]*f+R[5]*pp[t].y)/Z_up;
	Ai[6+2]=(R[7]*f+R[8]*pp[t].y)/Z_up;
	Ai[6+3]=-pp[t].x*sin(omega)-(pp[t].y/f*(pp[t].x*cos(kappa)-pp[t].y*sin(kappa))-f*sin(kappa))*cos(omega);
	Ai[6+4]=-f*cos(kappa)-pp[t].y/f*(pp[t].x*sin(kappa)+pp[t].y*cos(kappa));
	Ai[6+5]=-pp[t].x;

	//计算常数项
	Li[0]=pp[t].x+f*X_up/Z_up;
	Li[1]=pp[t].y+f*Y_up/Z_up;
}

//矩阵转置
double* matrixTransport(double *a, int m, int n) {
	double *newMatrix = new double[m*n];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			newMatrix[j*m+i] = a[i*n+j];
		}
	}
	return newMatrix;
}

//矩阵相乘：C=A x B，A为m,t阶矩阵,B为t,n阶矩阵,C为m,n阶矩阵
double* matrix_multiplication(double *A,double *B,int m,int t,int n){
	double *C = new double[m*n];
	//对C矩阵进行初始化
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			C[i*n+j]=0;
		}
	}

	for(int i=0;i<m;i++){//前面矩阵的行数
		for(int j=0;j<n;j++){//后面矩阵的列数
			for(int k=0;k<t;k++){//前面矩阵的列数
				C[i*n+j] += A[i*t+k]*B[k*n+j];
			}
		}
	}
	return C;
}

//矩阵求逆
void swap(double &a,double &b){
	double c=a;a=b;b=c;
}
int matrix_Inverse(double *A,int n){
	int i,j,k;
	double d;
	int *JS=new int[n];
	int *IS=new int[n];
	for (k=0;k<n;k++){
		d=0;
		for (i=k;i<n;i++){
			for (j=k;j<n;j++){
				if (fabs(A[i*n+j])>d){
					d=fabs(A[i*n+j]);
					IS[k]=i;
					JS[k]=j;
				}
			}
		}
		if (d+1.0==1.0) return 0;
		if (IS[k]!=k){
			for (j=0;j<n;j++){swap(A[k*n+j],A[IS[k]*n+j]);}
		}
		if (JS[k]!=k){
			for (i=0;i<n;i++){swap(A[i*n+k],A[i*n+JS[k]]);}
		}
		A[k*n+k]=1/A[k*n+k];
		for (j=0;j<n;j++){if (j!=k) A[k*n+j]=A[k*n+j]*A[k*n+k];}
		for (i=0;i<n;i++){
			if (i!=k){
				for (j=0;j<n;j++){if (j!=k) A[i*n+j]=A[i*n+j]-A[i*n+k]*A[k*n+j];}
			}
		}
		for (i=0;i<n;i++){if (i!=k) A[i*n+k]=-A[i*n+k]*A[k*n+k];}
	}
	for (k=n-1;k>=0;k--){
		for (j=0;j<n;j++){if (JS[k]!=k) swap(A[k*n+j],A[JS[k]*n+j]);}
		for (i=0;i<n;i++){if (IS[k]!=k) swap(A[i*n+k],A[i*n+IS[k]]);}
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	//计算算法时间效率的参数
	time_t start_time,end_time;
	//解算程序起始时间点
	start_time=clock();

	//已知数据：影像比例尺m，内方位元素x_0,y_0,f，未知数的个数point_number
	double m,x_0,y_0,f;
	int point_number;

	//未知数
	double Xs, Ys, Zs, H;
	double phi, omega, kappa;

	//从输入文件中读取数据并对未知数初始化
	FILE *fp; 
    fp = fopen("d:\\data.txt", "r"); 
    if  (fp == NULL) {
		cout<<"error:can`t open the file"<<endl;
    } 
	fscanf(fp,"%lf,%lf,%lf,%lf,%d", &m,&x_0,&y_0,&f,&point_number);
	PtoP *p = new PtoP[point_number];
	for (int i = 0;i < point_number;i++) {
		fscanf(fp,"%lf,%lf,%lf,%lf,%lf", &p[i].x,&p[i].y,&p[i].X,&p[i].Y,&p[i].Z);
	}
	//将f和x，y的mm单位统一为m
	f=f*0.001;
	for (int i = 0;i < point_number;i++) {
		p[i].x=p[i].x*0.001;
		p[i].y=p[i].y*0.001;
	}
	set_initial_unknownvalue(p,point_number,&m,&f,&Xs,&Ys,&Zs,&phi,&omega,&kappa);
	
	//定义最大迭代次数为15次
	int iteration_times=15;
	//定义法方程第i个未知数的权倒数的矩阵Qii[6]
	double Qii[6]={0,0,0,0,0,0};
	double *V = new double[point_number*2];
	double vv=0;
	double MEPUW=0;//单位权中误差Median error per unit weight
	double RMSE_Xs=0,RMSE_Ys=0,RMSE_Zs=0,RMSE_phi=0,RMSE_omega=0,RMSE_kappa=0;//6个外方位元素的中误差root mean square error

	do{
		vv=0;
		MEPUW=0;
		//定义旋转矩阵R
		double R[9];
		RotationMatrix(phi,omega,kappa,R);
	
		//输出R矩阵检查结果
		/*cout<<R[0]<<"  "<<R[1]<<"  "<<R[2]<<endl;
		cout<<R[3]<<"  "<<R[4]<<"  "<<R[5]<<endl;
		cout<<R[6]<<"  "<<R[7]<<"  "<<R[8]<<endl;*/

		double *A = new double[point_number*2*6];
		double *L = new double[point_number*2];
	
		//逐个点建立总的误差方程矩阵形式
		for(int i = 0;i < point_number;i++){
			double Ai[12];
			double Li[2];
			ErrorEquation(p,R,Ai,Li,i,Xs,Ys,Zs,phi,omega,kappa,f);
			for(int j=0;j<2*6;j++){
				A[i*12+j]=Ai[j];
			}
			for(int k=0;k<2;k++){
				L[i*2+k]=Li[k];
			}
		}

		//输出A矩阵检验结果
		/*for(int i=0;i<point_number*2;i++){
			for(int j=0;j<6;j++){
				cout<<A[i*6+j]<<"  ";
			}
			cout<<endl;
		}
		cout<<endl;*/

		double *At = new double[point_number*2*6];
		At = matrixTransport(A,2*point_number,6);
	
		//输出At检验结果
		/*cout<<"At"<<endl;
		for(int i=0;i<6;i++){
			for(int j=0;j<point_number*2;j++){
				cout<<At[i*point_number*2+j]<<"  ";
			}
			cout<<endl;
		}
		cout<<endl;*/

		double *AtA = new double[6*6];
		AtA=matrix_multiplication(At,A,6,point_number*2,6);

		//输出AtA矩阵检验结果
		/*cout<<"AtA"<<endl;
		for(int i=0;i<6;i++){
			for(int j=0;j<6;j++){
				cout<<AtA[i*6+j]<<"  ";
			}
			cout<<endl;
		}
		cout<<endl;*/

		double *AtA_Inverse=new double[6*6];
		for(int i=0;i<6;i++){
			for(int j=0;j<6;j++){
				AtA_Inverse[i*6+j]=AtA[i*6+j];
			}
		}
		matrix_Inverse(AtA_Inverse,6);
		//输出AtA的逆的矩阵检验结果
		/*cout<<"AtA_Inverse"<<endl;
		for(int i=0;i<6;i++){
			for(int j=0;j<6;j++){
				cout<<AtA_Inverse[i*6+j]<<"  ";
			}
			cout<<endl;
		}
		cout<<endl;*/
		
		//取出AtA_Inverse的主对角线元素作为法方程第i个未知数的权倒数的矩阵Qii
		for(int i=0;i<6;i++){Qii[i] = AtA_Inverse[7*i];}
		//输出Qii矩阵检验结果
		/*for(int i=0;i<6;i++){cout<<Qii[i]<<"    ";}
		cout<<endl;*/

		double *AtL = new double[6*1];
		AtL=matrix_multiplication(At,L,6,point_number*2,1);
		
		double *x = new double[6*1];
		x=matrix_multiplication(AtA_Inverse,AtL,6,6,1);
		//输出x矩阵检验结果
		/*for(int i=0;i<6;i++){
			cout<<x[i]<<"    ";
		}*/

		double *AX = new double[point_number*2*1];
		AX=matrix_multiplication(A,x,point_number*2,6,1);

		for (int i=0;i<point_number*2;i++){
			V[i]=AX[i]-L[i];
			vv+=pow(V[i],2);
		}

		MEPUW=sqrt(vv/(2*point_number-6));
		RMSE_Xs=MEPUW*sqrt(Qii[0]);
		RMSE_Ys=MEPUW*sqrt(Qii[1]);
		RMSE_Zs=MEPUW*sqrt(Qii[2]);
		RMSE_phi=MEPUW*sqrt(Qii[3]);
		RMSE_omega=MEPUW*sqrt(Qii[4]);
		RMSE_kappa=MEPUW*sqrt(Qii[5]);

		//外方位元素近似值与改正数求和得到新的近似值
		Xs+=x[0];
		Ys+=x[1];
		Zs+=x[2];
		phi+=x[3];
		omega+=x[4];
		kappa+=x[5];

		//输出每次迭代的结果作为检验
		/*cout<<"time"<<15-iteration_times<<endl;
		cout<<"dertaXs= "<<x[0]<<endl;
		cout<<"dertaYs= "<<x[1]<<endl;
		cout<<"dertaZs= "<<x[2]<<endl;
		cout<<"dertaphi= "<<x[3]<<endl;
		cout<<"dertaomega= "<<x[4]<<endl;
		cout<<"dertakappa= "<<x[5]<<endl;*/

		//判断外方位元素的改正数与限差的关系

		if((abs(x[0])<1e-3)&&(abs(x[1])<1e-3)&&(abs(x[2])<1e-3)&&(abs(x[3])<1e-6)&&(abs(x[4])<1e-6)&&(abs(x[5])<1e-6)){
			break;
		}

		delete[] AtA_Inverse,AtA,AtL,p,A,L,At,x,AX;
		iteration_times--;

		if(iteration_times==0){cout<<"结果是不收敛的";}

	}while(iteration_times>0);

	delete[] V;
	
	//解算程序终止时间点
	end_time=clock();

	if(iteration_times>0){
		//程序完成，开始输出结果及其数据文件
		//建立一个输出文件的操作指针
		FILE *fp2=fopen("d:\\Space_Resection.Output.txt","w");
		time_t timep;
		struct tm *ptime;
		time(&timep);
		ptime = gmtime(&timep);

		cout<<"**********这里是单像空间后方交会程序的结果显示界面**********"<<endl;
		cout<<"生成此报告的时间是："<<(1900 + ptime->tm_year)<<"年"<<(1 + ptime->tm_mon)<<"月"<<ptime->tm_mday<<"日 "<<ptime->tm_hour+8<<":"<<ptime->tm_min<<":"<<ptime->tm_sec<<endl;
		cout<<"单像空间后方交会程序运行共花费时间为："<<difftime(end_time,start_time)<<" ms"<<endl<<endl; 
		
		fprintf(fp2,"**********这里是单像空间后方交会程序的结果显示界面**********\n");
		fprintf(fp2,"生成此报告的时间是：%d年%d月%d日 %d:%d:%d\n",(1900 + ptime->tm_year),(1 + ptime->tm_mon),ptime->tm_mday,ptime->tm_hour+8,ptime->tm_min,ptime->tm_sec);
		fprintf(fp2,"单像空间后方交会程序运行共花费时间为：%d ms\n\n",difftime(end_time,start_time));

		cout<<"使用到的像点坐标和地面坐标的信息如下所示:"<<endl;
		for(int i=0;i<point_number;i++){
			cout<<"第"<<i+1<<"个坐标对为：x="<<p[i].x<<" y="<<p[i].y<<" X="<<p[i].X<<" Y="<<p[i].Y<<" Z="<<p[i].Z<<endl;
			fprintf(fp2,"第%d个坐标对为：x=%lf y=%lf X=%lf Y=%lf Z=%lf\n",i+1,p[i].x,p[i].y,p[i].X,p[i].Y,p[i].Z);
		}
		
		cout<<endl<<"在解算中一共经过了"<<(15+1)-iteration_times<<"次迭代计算,解算得到的六个外方位元素的结果为："<<endl;
		cout<<"Xs="<<Xs<<"m,Ys="<<Ys<<"m,Zs="<<Zs<<"m"<<endl<<"phi="<<phi<<",omega="<<omega<<",kappa="<<kappa<<endl;
		fprintf(fp2,"\n在解算中一共经过了%d次迭代计算，解算得到的六个外方位元素的结果为：\n", (15+1)-iteration_times);
		fprintf(fp2, "Xs=%lfm,Ys=%lfm,Zs=%lfm\nPhi=%lf rad,Omega=%lf rad,Kappa=%lf rad\n",Xs,Ys,Zs,phi,omega,kappa);

		double R_final[9];
		RotationMatrix(phi,omega,kappa,R_final);
		cout<<endl<<"由解算得到的外方位元素计算得到的R矩阵为:"<<endl;
		cout<<R_final[0]<<"  "<<R_final[1]<<"  "<<R_final[2]<<endl;
		cout<<R_final[3]<<"  "<<R_final[4]<<"  "<<R_final[5]<<endl;
		cout<<R_final[6]<<"  "<<R_final[7]<<"  "<<R_final[8]<<endl;
		fprintf(fp2, "\n由解算得到的外方位元素计算得到的R矩阵为:\n %lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n\n",R_final[0],R_final[1],R_final[2],R_final[3],R_final[4],R_final[5],R_final[6],R_final[7],R_final[8]);

		cout<<endl<<"单位权中误差为："<<MEPUW<<endl;
		cout<<"六个外方位元素的解算精度为："<<endl;
		cout<<"Xs的中误差RMSE_Xs(mXs)= "<<RMSE_Xs<<endl;
		cout<<"Ys的中误差RMSE_Ys(mYs)= "<<RMSE_Ys<<endl;
		cout<<"Zs的中误差RMSE_Zs(mZs)= "<<RMSE_Zs<<endl;
		cout<<"phi的中误差RMSE_phi(mphi)= "<<RMSE_phi<<endl;
		cout<<"omega的中误差RMSE_omega(momega)= "<<RMSE_omega<<endl;
		cout<<"kappa的中误差RMSE_kappa(mkappa)= "<<RMSE_kappa<<endl;

		fprintf(fp2, "单位权中误差为：%lf\n",MEPUW);
		fprintf(fp2, "六个外方位元素的解算精度为：\nXs的中误差RMSE_Xs(mXs)= %lf\nYs的中误差RMSE_Ys(mYs)= %lf\nZs的中误差RMSE_Zs(mZs)= %lf\nphi的中误差RMSE_phi(mphi)= %lf\nomega的中误差RMSE_omega(momega)= %lf\nkappa的中误差RMSE_kappa(mkappa)= %lf\n",RMSE_Xs,RMSE_Ys,RMSE_Zs,RMSE_phi,RMSE_omega,RMSE_kappa);

		cout<<endl<<"**********结果显示完毕，感谢使用**********"<<endl;
		fprintf(fp2,"\n**********结果显示完毕，感谢使用**********");

		fclose(fp2);
	}
	fclose(fp);
	system("pause");
	return 0;
}

