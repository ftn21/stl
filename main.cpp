#include <QCoreApplication>
#include <iostream>
#include "linalg.h"
#include <cstring>
#include <math.h>

#include <QVector>
#include <QDebug>

#define pi 3.141592653589793
using namespace std;
using namespace MyLinearAlgebra;

/////////////////////////////////////////////////////////////////////////////////////////////////

const double C = 299792458; // скорость света [м/с]

//точка
TVector P (3);

//спутники
TVector S08 (3);
TVector S11 (3);
TVector S20 (3);
TVector S21 (3);
TVector S23 (3);
TVector S27 (3);

TVector C1 (6);     //псевдодальность [м]
TVector dt (6);     //оценка времени полета сигнала
TVector dT (6);     //сдвиг часов [мкс]
TVector D1(6);      //доплер [Гц]

//sattelite s08;

int main(int argc, char *argv[]) {
    QCoreApplication a(argc, argv);

    //координаты

    P[0] = 0;
    P[1] = 0;
    P[2] = 0;

    S08[0] = 16678843.241;
    S08[1] = -2634637.336;
    S08[2] = 20579559.047;

    S11[0] = 13020628.571;
    S11[1] = -12527567.852;
    S11[2] = 19099436.863;

    S20[0] = -2555568.528;
    S20[1] = 16762578.066;
    S20[2] = 20311833.019;

    S21[0] = 16181563.666;
    S21[1] = -957262.138;
    S21[2] = 21758956.397;

    S23[0] = -1405082.371;
    S23[1] = 15628413.976;
    S23[2] = 21394716.744;

    S27[0] = 20637508.416;
    S27[1] = 8906797.175;
    S27[2] = 14267852.574;

    C1[0] = 21293946.348;
    C1[1] = 22744995.395;
    C1[2] = 21565628.372;
    C1[3] = 21513278.753;
    C1[4] = 21488700.402;
    C1[5] = 21246116.669;

    dt[0] = C1[0] / C;
    dt[1] = C1[1] / C;
    dt[2] = C1[2] / C;
    dt[3] = C1[3] / C;
    dt[4] = C1[4] / C;
    dt[5] = C1[5] / C;

    D1[0] = 2677.134;
    D1[1] = 5331.631;
    D1[2] = -484.137;
    D1[3] = 3400.669;
    D1[4] = -347.145;
    D1[5] = 143.621;

    //последовательные приближения
    for (quint8 it = 0; it < 8; it++) {

        TMatrix H (6,4);

        TVector r08 = S08 - P;
        H(0, 0) = -r08[0] / r08.length();
        H(0, 1) = -r08[1] / r08.length();
        H(0, 2) = -r08[2] / r08.length();
        H(0, 3) = 1.0;

        TVector r11 = S11 - P;
        H(1, 0) = -r11[0] / r11.length();
        H(1, 1) = -r11[1] / r11.length();
        H(1, 2) = -r11[2] / r11.length();
        H(1, 3) = 1; //поправка в метрах и считать сразу с*deltaT

        TVector r20 = S20 - P;
        H(2, 0) = -r20[0]/ r20.length();
        H(2, 1) = -r20[1] / r20.length();
        H(2, 2) = -r20[2] / r20.length();
        H(2, 3) = 1;

        TVector r21 = S21 - P;
        H(3, 0) = -r21[0] / r21.length();
        H(3, 1) = -r21[1] / r21.length();
        H(3, 2) = -r21[2] / r21.length();
        H(3, 3) = 1;

        TVector r23 = S23 - P;
        H(4, 0) = -r23[0] / r23.length();
        H(4, 1) = -r23[1] / r23.length();
        H(4, 2) = -r23[2] / r23.length();
        H(4, 3) = 1;

        TVector r27 = S27 - P;
        H(5, 0) = -r27[0] / r27.length();
        H(5, 1) = -r27[1] / r27.length();
        H(5, 2) = -r27[2] / r27.length();
        H(5, 3) = 1;

        //матрица невязок
        TMatrix V(6, 1);

        //вращение приемника вместе с землёй
        TVector Ve(3);
        Ve[0] = 0;
        Ve[1] = 0;
        Ve[2] = 7.2921151467 / 100000;
        Ve = Ve * P;

        //рассточние приемник-спутник
        TVector R(6);
        R = S08 - P;

        //нормаль на спутник
        TVector e = R.norm();

        //double Vr =

        //           C1                        dT                       D
        V(0, 0) = C1[0] - r08.length() - 50.014071*(C/1000000)  - (D1[0]*0.1902)*dt[0];
        V(1, 0) = C1[1] - r11.length() - 148.368549*(C/1000000) - (D1[1]*0.1902)*dt[1];
        V(2, 0) = C1[2] - r20.length() + 526.543696*(C/1000000) - (D1[2]*0.1902)*dt[2];
        V(3, 0) = C1[3] - r21.length() + 53.373334*(C/1000000)  - (D1[3]*0.1902)*dt[3];
        V(4, 0) = C1[4] - r23.length() + 78.206308*(C/1000000)  - (D1[4]*0.1902)*dt[4];
        V(5, 0) = C1[5] - r27.length() - 409.011384*(C/1000000) - (D1[5]*0.1902)*dt[5];

        /*cout << "V[0,0] = " << V(0,0) << endl;
        cout << "V[1,0] = " << V(1,0) << endl;
        cout << "V[2,0] = " << V(2,0) << endl;
        cout << "V[3,0] = " << V(3,0) << endl;
        cout << "V[4,0] = " << V(4,0) << endl;
        cout << "V[5,0] = " << V(5,0) << endl;
        */

        //MHK: dV = (H^t * H)^(-1) * H^t * V
        TMatrix Ht = H.t();
        TMatrix HtH = Ht*H;

        TMatrix Q = !HtH;
        TMatrix dV = Q*Ht*V;

        //вычисление поправки
        double dx = dV(0, 0);
        double dy = dV(1, 0);
        double dz = dV(2, 0);

        //новое приближение
        P[0] += dx;
        P[1] += dy;
        P[2] += dz;

        if ( sqrt(dx*dx + dy*dy + dz*dz) < 0.00001 ) {
            break;
        }
        else {
            cout << "calculating. . .";
        }

    }

    TVector approx_pos (3);
    approx_pos[0] = 2849830.5060;
    approx_pos[1] = 2186822.2813;
    approx_pos[2] = 5252937.0124;

    TVector error = P - approx_pos;

    cout << endl << "ошибка = " << error.length() << endl << endl;

    return a.exec();
}
