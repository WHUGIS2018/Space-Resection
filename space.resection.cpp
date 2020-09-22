//����ʹ�ñ�ע
/*
�˳���ͨ��C/C++����ʵ�ֵ���ռ�󷽽���ļ�������򣬳�����뻷��ΪWin10+vs2010����̨��
�ο���ĿΪ��Ӱ����ѧ�ڶ��棨�Ž��������P44ҳ��27�⣬
�����������txt�ļ��������ʽ������Ϊdata.txt�������D���У��ļ��е����ݸ�ʽΪ��
  ������m���ڷ�λԪ��x0,y0,f�����Ƶ�ĸ���
  �����Ƶ��Ӱ������͵�������(x,y,X,Y,Z)
�����е������ļ�����Ϊ��
50000,0,0,153.24,4
-86.15,-68.99,36589.41,25273.32,2195.17
-53.40,82.21,37631.08,31324.51,728.69
-14.78,-76.63,39100.97,24934.98,2386.50
10.46,64.43,40426.54,30319.81,757.31
�����������txt�ļ��������ʽ������ΪSpace_Resection.Output.txt�������D����
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

//�����Ե�Ӱ������͵�����������ݽṹ
struct PtoP {
	double x;
	double y;
	double X;
	double Y;
	double Z;
};

//ȷ��δ֪���ĳ�ֵ
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

//������ת����
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

//����һ�Ե������V=AX-L�Ĳ����ͳ�����
void ErrorEquation(PtoP *pp,double R[],double Ai[12],double Li[2],int t,double Xs,double Ys,double Zs,double phi,double omega,double kappa,double f){
	//��������X��Y��Z���Ϻ�ܣ�
	double X_up,Y_up,Z_up;
	X_up= R[0]*(pp[t].X-Xs)+R[3]*(pp[t].Y-Ys)+R[6]*(pp[t].Z-Zs);
	Y_up= R[1]*(pp[t].X-Xs)+R[4]*(pp[t].Y-Ys)+R[7]*(pp[t].Z-Zs);
	Z_up= R[2]*(pp[t].X-Xs)+R[5]*(pp[t].Y-Ys)+R[8]*(pp[t].Z-Zs);

	//����x���ֵ�ƫ����a11��a16:
	Ai[0+0]=(R[0]*f+R[2]*pp[t].x)/Z_up;
	Ai[0+1]=(R[3]*f+R[5]*pp[t].x)/Z_up;
	Ai[0+2]=(R[6]*f+R[8]*pp[t].x)/Z_up;
	Ai[0+3]=pp[t].y*sin(omega)-(pp[t].x/f*(pp[t].x*cos(kappa)-pp[t].y*sin(kappa))+f*cos(kappa))*cos(omega);
	Ai[0+4]=-f*sin(kappa)-pp[t].x/f*(pp[t].x*sin(kappa)+pp[t].y*cos(kappa));
	Ai[0+5]=pp[t].y;

	//����y���ֵ�ƫ����a21��a26:
	Ai[6+0]=(R[1]*f+R[2]*pp[t].y)/Z_up;
	Ai[6+1]=(R[4]*f+R[5]*pp[t].y)/Z_up;
	Ai[6+2]=(R[7]*f+R[8]*pp[t].y)/Z_up;
	Ai[6+3]=-pp[t].x*sin(omega)-(pp[t].y/f*(pp[t].x*cos(kappa)-pp[t].y*sin(kappa))-f*sin(kappa))*cos(omega);
	Ai[6+4]=-f*cos(kappa)-pp[t].y/f*(pp[t].x*sin(kappa)+pp[t].y*cos(kappa));
	Ai[6+5]=-pp[t].x;

	//���㳣����
	Li[0]=pp[t].x+f*X_up/Z_up;
	Li[1]=pp[t].y+f*Y_up/Z_up;
}

//����ת��
double* matrixTransport(double *a, int m, int n) {
	double *newMatrix = new double[m*n];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			newMatrix[j*m+i] = a[i*n+j];
		}
	}
	return newMatrix;
}

//������ˣ�C=A x B��AΪm,t�׾���,BΪt,n�׾���,CΪm,n�׾���
double* matrix_multiplication(double *A,double *B,int m,int t,int n){
	double *C = new double[m*n];
	//��C������г�ʼ��
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			C[i*n+j]=0;
		}
	}

	for(int i=0;i<m;i++){//ǰ����������
		for(int j=0;j<n;j++){//������������
			for(int k=0;k<t;k++){//ǰ����������
				C[i*n+j] += A[i*t+k]*B[k*n+j];
			}
		}
	}
	return C;
}

//��������
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
	//�����㷨ʱ��Ч�ʵĲ���
	time_t start_time,end_time;
	//���������ʼʱ���
	start_time=clock();

	//��֪���ݣ�Ӱ�������m���ڷ�λԪ��x_0,y_0,f��δ֪���ĸ���point_number
	double m,x_0,y_0,f;
	int point_number;

	//δ֪��
	double Xs, Ys, Zs, H;
	double phi, omega, kappa;

	//�������ļ��ж�ȡ���ݲ���δ֪����ʼ��
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
	//��f��x��y��mm��λͳһΪm
	f=f*0.001;
	for (int i = 0;i < point_number;i++) {
		p[i].x=p[i].x*0.001;
		p[i].y=p[i].y*0.001;
	}
	set_initial_unknownvalue(p,point_number,&m,&f,&Xs,&Ys,&Zs,&phi,&omega,&kappa);
	
	//��������������Ϊ15��
	int iteration_times=15;
	//���巨���̵�i��δ֪����Ȩ�����ľ���Qii[6]
	double Qii[6]={0,0,0,0,0,0};
	double *V = new double[point_number*2];
	double vv=0;
	double MEPUW=0;//��λȨ�����Median error per unit weight
	double RMSE_Xs=0,RMSE_Ys=0,RMSE_Zs=0,RMSE_phi=0,RMSE_omega=0,RMSE_kappa=0;//6���ⷽλԪ�ص������root mean square error

	do{
		vv=0;
		MEPUW=0;
		//������ת����R
		double R[9];
		RotationMatrix(phi,omega,kappa,R);
	
		//���R��������
		/*cout<<R[0]<<"  "<<R[1]<<"  "<<R[2]<<endl;
		cout<<R[3]<<"  "<<R[4]<<"  "<<R[5]<<endl;
		cout<<R[6]<<"  "<<R[7]<<"  "<<R[8]<<endl;*/

		double *A = new double[point_number*2*6];
		double *L = new double[point_number*2];
	
		//����㽨���ܵ����̾�����ʽ
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

		//���A���������
		/*for(int i=0;i<point_number*2;i++){
			for(int j=0;j<6;j++){
				cout<<A[i*6+j]<<"  ";
			}
			cout<<endl;
		}
		cout<<endl;*/

		double *At = new double[point_number*2*6];
		At = matrixTransport(A,2*point_number,6);
	
		//���At������
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

		//���AtA���������
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
		//���AtA����ľ��������
		/*cout<<"AtA_Inverse"<<endl;
		for(int i=0;i<6;i++){
			for(int j=0;j<6;j++){
				cout<<AtA_Inverse[i*6+j]<<"  ";
			}
			cout<<endl;
		}
		cout<<endl;*/
		
		//ȡ��AtA_Inverse�����Խ���Ԫ����Ϊ�����̵�i��δ֪����Ȩ�����ľ���Qii
		for(int i=0;i<6;i++){Qii[i] = AtA_Inverse[7*i];}
		//���Qii���������
		/*for(int i=0;i<6;i++){cout<<Qii[i]<<"    ";}
		cout<<endl;*/

		double *AtL = new double[6*1];
		AtL=matrix_multiplication(At,L,6,point_number*2,1);
		
		double *x = new double[6*1];
		x=matrix_multiplication(AtA_Inverse,AtL,6,6,1);
		//���x���������
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

		//�ⷽλԪ�ؽ���ֵ���������͵õ��µĽ���ֵ
		Xs+=x[0];
		Ys+=x[1];
		Zs+=x[2];
		phi+=x[3];
		omega+=x[4];
		kappa+=x[5];

		//���ÿ�ε����Ľ����Ϊ����
		/*cout<<"time"<<15-iteration_times<<endl;
		cout<<"dertaXs= "<<x[0]<<endl;
		cout<<"dertaYs= "<<x[1]<<endl;
		cout<<"dertaZs= "<<x[2]<<endl;
		cout<<"dertaphi= "<<x[3]<<endl;
		cout<<"dertaomega= "<<x[4]<<endl;
		cout<<"dertakappa= "<<x[5]<<endl;*/

		//�ж��ⷽλԪ�صĸ��������޲�Ĺ�ϵ

		if((abs(x[0])<1e-3)&&(abs(x[1])<1e-3)&&(abs(x[2])<1e-3)&&(abs(x[3])<1e-6)&&(abs(x[4])<1e-6)&&(abs(x[5])<1e-6)){
			break;
		}

		delete[] AtA_Inverse,AtA,AtL,p,A,L,At,x,AX;
		iteration_times--;

		if(iteration_times==0){cout<<"����ǲ�������";}

	}while(iteration_times>0);

	delete[] V;
	
	//���������ֹʱ���
	end_time=clock();

	if(iteration_times>0){
		//������ɣ���ʼ���������������ļ�
		//����һ������ļ��Ĳ���ָ��
		FILE *fp2=fopen("d:\\Space_Resection.Output.txt","w");
		time_t timep;
		struct tm *ptime;
		time(&timep);
		ptime = gmtime(&timep);

		cout<<"**********�����ǵ���ռ�󷽽������Ľ����ʾ����**********"<<endl;
		cout<<"���ɴ˱����ʱ���ǣ�"<<(1900 + ptime->tm_year)<<"��"<<(1 + ptime->tm_mon)<<"��"<<ptime->tm_mday<<"�� "<<ptime->tm_hour+8<<":"<<ptime->tm_min<<":"<<ptime->tm_sec<<endl;
		cout<<"����ռ�󷽽���������й�����ʱ��Ϊ��"<<difftime(end_time,start_time)<<" ms"<<endl<<endl; 
		
		fprintf(fp2,"**********�����ǵ���ռ�󷽽������Ľ����ʾ����**********\n");
		fprintf(fp2,"���ɴ˱����ʱ���ǣ�%d��%d��%d�� %d:%d:%d\n",(1900 + ptime->tm_year),(1 + ptime->tm_mon),ptime->tm_mday,ptime->tm_hour+8,ptime->tm_min,ptime->tm_sec);
		fprintf(fp2,"����ռ�󷽽���������й�����ʱ��Ϊ��%d ms\n\n",difftime(end_time,start_time));

		cout<<"ʹ�õ����������͵����������Ϣ������ʾ:"<<endl;
		for(int i=0;i<point_number;i++){
			cout<<"��"<<i+1<<"�������Ϊ��x="<<p[i].x<<" y="<<p[i].y<<" X="<<p[i].X<<" Y="<<p[i].Y<<" Z="<<p[i].Z<<endl;
			fprintf(fp2,"��%d�������Ϊ��x=%lf y=%lf X=%lf Y=%lf Z=%lf\n",i+1,p[i].x,p[i].y,p[i].X,p[i].Y,p[i].Z);
		}
		
		cout<<endl<<"�ڽ�����һ��������"<<(15+1)-iteration_times<<"�ε�������,����õ��������ⷽλԪ�صĽ��Ϊ��"<<endl;
		cout<<"Xs="<<Xs<<"m,Ys="<<Ys<<"m,Zs="<<Zs<<"m"<<endl<<"phi="<<phi<<",omega="<<omega<<",kappa="<<kappa<<endl;
		fprintf(fp2,"\n�ڽ�����һ��������%d�ε������㣬����õ��������ⷽλԪ�صĽ��Ϊ��\n", (15+1)-iteration_times);
		fprintf(fp2, "Xs=%lfm,Ys=%lfm,Zs=%lfm\nPhi=%lf rad,Omega=%lf rad,Kappa=%lf rad\n",Xs,Ys,Zs,phi,omega,kappa);

		double R_final[9];
		RotationMatrix(phi,omega,kappa,R_final);
		cout<<endl<<"�ɽ���õ����ⷽλԪ�ؼ���õ���R����Ϊ:"<<endl;
		cout<<R_final[0]<<"  "<<R_final[1]<<"  "<<R_final[2]<<endl;
		cout<<R_final[3]<<"  "<<R_final[4]<<"  "<<R_final[5]<<endl;
		cout<<R_final[6]<<"  "<<R_final[7]<<"  "<<R_final[8]<<endl;
		fprintf(fp2, "\n�ɽ���õ����ⷽλԪ�ؼ���õ���R����Ϊ:\n %lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n\n",R_final[0],R_final[1],R_final[2],R_final[3],R_final[4],R_final[5],R_final[6],R_final[7],R_final[8]);

		cout<<endl<<"��λȨ�����Ϊ��"<<MEPUW<<endl;
		cout<<"�����ⷽλԪ�صĽ��㾫��Ϊ��"<<endl;
		cout<<"Xs�������RMSE_Xs(mXs)= "<<RMSE_Xs<<endl;
		cout<<"Ys�������RMSE_Ys(mYs)= "<<RMSE_Ys<<endl;
		cout<<"Zs�������RMSE_Zs(mZs)= "<<RMSE_Zs<<endl;
		cout<<"phi�������RMSE_phi(mphi)= "<<RMSE_phi<<endl;
		cout<<"omega�������RMSE_omega(momega)= "<<RMSE_omega<<endl;
		cout<<"kappa�������RMSE_kappa(mkappa)= "<<RMSE_kappa<<endl;

		fprintf(fp2, "��λȨ�����Ϊ��%lf\n",MEPUW);
		fprintf(fp2, "�����ⷽλԪ�صĽ��㾫��Ϊ��\nXs�������RMSE_Xs(mXs)= %lf\nYs�������RMSE_Ys(mYs)= %lf\nZs�������RMSE_Zs(mZs)= %lf\nphi�������RMSE_phi(mphi)= %lf\nomega�������RMSE_omega(momega)= %lf\nkappa�������RMSE_kappa(mkappa)= %lf\n",RMSE_Xs,RMSE_Ys,RMSE_Zs,RMSE_phi,RMSE_omega,RMSE_kappa);

		cout<<endl<<"**********�����ʾ��ϣ���лʹ��**********"<<endl;
		fprintf(fp2,"\n**********�����ʾ��ϣ���лʹ��**********");

		fclose(fp2);
	}
	fclose(fp);
	system("pause");
	return 0;
}

