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

const double C = 299.792458;

//точка
TVector P (3);

//спутники
TVector S08 (3);
TVector S11 (3);
TVector S20 (3);
TVector S21 (3);
TVector S23 (3);
TVector S27 (3);

int main(int argc, char *argv[]) {
    QCoreApplication a(argc, argv);

    //координаты

    P[0] = 0;
    P[1] = 0;
    P[2] = 0;

    S08[0] = 16678843;
    S08[1] = -2634637;
    S08[2] = 20579559;

    S11[0] = 13020629;
    S11[1] = -12527568;
    S11[2] = 19099437;

    S20[0] = -2555569;
    S20[1] = 16762578;
    S20[2] = 20311833;

    S21[0] = 16181564;
    S21[1] = -957262;
    S21[2] = 21758956;

    S23[0] = -1405082;
    S23[1] = 15628414;
    S23[2] = 21394717;

    S27[0] = 20637508;
    S27[1] = 8906797;
    S27[2] = 14267853;

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
        V(0, 0) = 21293946 - r08.length() - 50.014071*C;
        V(1, 0) = 22744995 - r11.length() - 148.368549*C;
        V(2, 0) = 21565628 - r20.length() + 526.543696*C;
        V(3, 0) = 21513278 - r21.length() + 53.373334*C;
        V(4, 0) = 21488700 - r23.length() + 78.206308*C;
        V(5, 0) = 21246117 - r27.length() - 409.011384*C;

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

        if (sqrt(dx*dx + dy*dy + dz*dz) < 0.0001) {
            break;
        }
        else {
            cout << "calculating. . .";
        }

    }

    TVector approx_pos (3);
    approx_pos[0] = 2849831;
    approx_pos[1] = 2186822;
    approx_pos[2] = 5252937;

    TVector error = P - approx_pos;
    int err = error.length();

    cout << "ошибка = " << err << endl << endl;

    return a.exec();
}
