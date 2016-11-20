// Copyright <duan> 2016.11.15
#include<iostream>
#include<cmath>
using namespace std;
#include<stdio.h>

int times = 0;
double m1 = 0.5;
double m2 = 0.5;
double l1 = 0.15;
double l2 = 0.15;
double lg1 = 0.5*l1;
double lg2 = 0.5*l2;
double I1 = (1.0/3)*m1*l1*l1;
double I2 = (1.0/12)*m2*l2*l2;
double g = 9.8;
double derutaT = 0.001;
double theta1ref = 90.0*M_PI/180;
double theta2ref = 0.0*M_PI/180;
int period = 500;
int T = 2000;
double P_1 = 400;//350;//200;
double I_1 = 0.02;//0.03;
double D_1 = 20;//50;//17;
double P_2 = 400;//350;//400;
double I_2 = 0.1;//0.05;
double D_2 = 20;//50;//12;
double Error[2] = {0, 0};
double sumError[2] = {0, 0};
double trac[2];
double docTrac[2];
double doc2Trac[2];
int N = 2;
bool Gauss(double A[][2], double B[][2], int n) {
    int i, j, k;
    double max, temp;
    double t[N][N];              //临时矩阵
    // 将A矩阵存放在临时矩阵t[n][n]中
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        t[i][j] = A[i][j];
      }
    }
    // 初始化B矩阵为单位阵
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        B[i][j] = (i == j) ? (float)1 : 0;
      }
    }
    for (i = 0; i < n; i++) {
        //寻找主元
        max = t[i][i];
        k = i;
        for (j = i + 1; j < n; j++) {
          if (fabs(t[j][i]) > fabs(max)) {
            max = t[j][i];
            k = j;
          }
        }
        // 如果主元所在行不是第i行，进行行交换
        if (k != i) {
          for (j = 0; j < n; j++) {
            temp = t[i][j];
            t[i][j] = t[k][j];
            t[k][j] = temp;
            // B伴随交换
            temp = B[i][j];
            B[i][j] = B[k][j];
            B[k][j] = temp;
          }
        }
        // 判断主元是否为0, 若是, 则矩阵A不是满秩矩阵,不存在逆矩阵
        if (t[i][i] == 0) {
            cout << "There is no inverse matrix!";
            return false;
        }
        // 消去A的第i列除去i行以外的各行元素
        temp = t[i][i];
        for (j = 0; j < n; j++) {
            t[i][j] = t[i][j] / temp;        // 主对角线上的元素变为1
            B[i][j] = B[i][j] / temp;        // 伴随计算
        }
        for (j = 0; j < n; j++) {
          if (j != i) {
                temp = t[j][i];
                for (k = 0; k < n; k++) {
                    t[j][k] = t[j][k] - t[i][k] * temp;
                    B[j][k] = B[j][k] - B[i][k] * temp;
                }
            }
        }
    }
    return true;
}


int main() {
  
  FILE *PIDdata1;
  if ((PIDdata1 = fopen("PIDdata1.txt", "w+")) == NULL) {
    cout<< "打开文件失败!";
  }
  FILE *PIDdata2;
  if ((PIDdata2 = fopen("PIDdata2.txt", "w+")) == NULL) {
    cout<< "打开文件失败!";
  }
  FILE *PIDdata3;
  if ((PIDdata3 = fopen("PIDdata3.txt", "w+")) == NULL) {
    cout<< "打开文件失败!";
  }
  FILE *PIDdata4;
  if ((PIDdata4 = fopen("PIDdata4.txt", "w+")) == NULL) {
    cout<< "打开文件失败!";
  }

  FILE *PIDdata5;
  if ((PIDdata5 = fopen("TrainSymble4.txt", "w+")) == NULL) {
    cout<< "打开文件失败!";
  }

  FILE *PIDdata6;
  if ((PIDdata6 = fopen("TrainResult4.txt", "w+")) == NULL) {
    cout<< "打开文件失败!";
  }

  double thetad1 = 30;
  double thetad2 = 30;
  double omega1 = 0;
  double omega2 = 0;
  double theta1 = 1.0*thetad1*M_PI/180;
  double theta2 = 1.0*thetad2*M_PI/180;
  double u1 = 0;
  double u2 = 0;
  double fx = l1*cos(theta1)+l2*cos(theta1+theta2);
  double fy = l1*sin(theta1)+l2*cos(theta1+theta2);
  double X[4][period];
  X[0][0] = theta1;
  X[1][0] = theta2;
  X[2][0] = omega1;
  X[3][0] = omega2;
  double M[2][2];
  double _M[2][2];
  double h[2];
  double G[2];
  // Get u1, u2
  for (int k = 0; k < period; k++) {

    // Error[0] = theta1ref - theta1;
    // Error[1] = theta2ref - theta2;
    trac[0] = M_PI/6+4*M_PI*pow((k*derutaT), 2)-16*M_PI/3*pow((k*derutaT), 3);
    trac[1] = M_PI/6-2*M_PI*pow((k*derutaT), 2)+8*M_PI/3*pow((k*derutaT), 3);
    docTrac[0] = 8*M_PI*(k*derutaT)-16*M_PI*pow((k*derutaT), 2);
    docTrac[1] = -4*M_PI*(k*derutaT)+8*M_PI*pow((k*derutaT), 2);
    doc2Trac[0] = 0;//M_PI/2-M_PI/2*(k*derutaT);
    doc2Trac[1] = 0;//-M_PI/4+M_PI/4*pow((k*derutaT), 2);
    Error[0] = trac[0] - theta1;
    Error[1] = trac[1] - theta2;
    sumError[0] += Error[0];
    sumError[1] += Error[1];
    u1 = P_1*Error[0]+D_1*(docTrac[0]-omega1)+I_1*sumError[0];
    u2 = P_2*Error[1]+D_2*(docTrac[1]-omega2)+I_2*sumError[1];
    M[0][0] = m1*pow(lg1, 2)+m2*pow(l1, 2)+m2*pow(lg2, 2)
          +I1+I2+2*m2*l1*lg2*cos(theta2);
    M[0][1] = m2*pow(lg2, 2)+I2+m2*l1*lg2*cos(theta2);
    M[1][0] = M[0][1];
    M[1][1] = m2*pow(lg2, 2)+I2;
    bool Is = Gauss(M, _M, 2);
    h[0] = -m2*l1*lg2*(2*omega1+omega2)*omega2*sin(theta2);
    h[1] = m2*l1*lg2*pow(omega1, 2)*sin(theta2);
    G[0] = (g*l1*m2+g*lg1*m1)*cos(theta1)+m2*g*lg2*cos(theta1+theta2);
    G[1] = m2*g*lg2*cos(theta1+theta2);
    double taw[2];
    taw[0] = M[0][0]*u1+M[0][1]*u2;
    taw[1] = M[1][0]*u1+M[1][1]*u2;
    // taw[0] = h[0]+G[0]-M[0][0]*u1-M[0][1]*u2+M[0][0]*doc2Trac[0]+M[0][1]*doc2Trac[1];
    // taw[1] = h[1]+G[1]-M[1][0]*u1-M[1][1]*u2+M[1][0]*doc2Trac[0]+M[1][1]*doc2Trac[1];
    double f[2];
    /*   f[0] = _M[0][0]*(-h[0]-G[0])+_M[0][1]*(-h[1]-G[1])
          +_M[0][0]*taw[0]+_M[0][1]*taw[1];
    f[1] = _M[1][0]*(-h[1]-G[1])+_M[1][1]*(-h[1]-G[1])
    +_M[1][0]*taw[0]+_M[1][1]*taw[1];*/
    f[0] = _M[0][0]*(-h[0]-G[0])+_M[0][1]*(-h[1]-G[1])+u1;
    f[1] = _M[1][0]*(-h[0]-G[0])+_M[1][1]*(-h[1]-G[1])+u2;
    X[0][k+1] = X[0][k]+derutaT*omega1;
    X[1][k+1] = X[1][k]+derutaT*omega2;
    X[2][k+1] = X[2][k]+derutaT*f[0];
    X[3][k+1] = X[3][k]+derutaT*f[1];
    
    theta1 = X[0][k+1];
    theta2 = X[1][k+1];
    omega1 = X[2][k+1];
    omega2 = X[3][k+1];
    fprintf(PIDdata3, "%f\n", taw[0]);
    fprintf(PIDdata4, "%f\n", taw[1]);
    if (k%5 == 0) {
      fprintf(PIDdata6, "%f %f\n", taw[0], taw[1]);
    }
  }

  for (int i = 0; i < period+1; i++) {
    theta1 = X[0][i];
    theta2 = X[1][i];
    omega1 = X[2][i];
    omega2 = X[3][i];
    double fx = l1*cos(theta1)+l2*cos(theta1+theta2);
    double fy = l1*sin(theta1)+l2*cos(theta1+theta2);
    fprintf(PIDdata1, "%f\n", theta1);  // *180/M_PI);
    fprintf(PIDdata2, "%f\n", theta2);  // *180/M_PI);
    if (i%5 == 0) {
      fprintf(PIDdata5, "%f %f %f %f\n", theta1, theta2, omega1, omega2);
    }
  }

  cout<< "Finial!"<< times<< endl;
  return 1;
}
