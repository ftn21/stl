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

struct sattelite {
public:
    string prn;
    TVector xyz; //координаты
    sattelite() : xyz(3) {};
    double C1,  //псевдодальность [м]
    dt,         //оценка времени полета сигнала
    dT,         //сдвиг часов [мкс]
    D1,         //доплер [Гц]
    P1;         //псевдодальность с коррекцией часов
};

class vec3 {
public:
    TVector xyz;
    vec3() : xyz(3) {};
};

const double C = 299792458; // скорость света [м/с]

//точка
TVector P (3);

//спутники
sattelite sats[6];

int main(int argc, char *argv[]) {
    QCoreApplication a(argc, argv);

    //prn
    sats[0].prn = "S08";
    sats[1].prn = "S11";
    sats[2].prn = "S20";
    sats[3].prn = "S21";
    sats[4].prn = "S23";
    sats[5].prn = "S27";

    //координаты
    P[0] = 0;
    P[1] = 0;
    P[2] = 0;

    sats[0].xyz[0] = 16678843.241;
    sats[0].xyz[1] = -2634637.336;
    sats[0].xyz[2] = 20579559.047;

    sats[1].xyz[0] = 13020628.571;
    sats[1].xyz[1] = -12527567.852;
    sats[1].xyz[2] = 19099436.863;

    sats[2].xyz[0] = -2555568.528;
    sats[2].xyz[1] = 16762578.066;
    sats[2].xyz[2] = 20311833.019;

    sats[3].xyz[0] = 16181563.666;
    sats[3].xyz[1] = -957262.138;
    sats[3].xyz[2] = 21758956.397;

    sats[4].xyz[0] = -1405082.371;
    sats[4].xyz[1] = 15628413.976;
    sats[4].xyz[2] = 21394716.744;

    sats[5].xyz[0] = 20637508.416;
    sats[5].xyz[1] = 8906797.175;
    sats[5].xyz[2] = 14267852.574;

    //псевдодальность [м]
    sats[0].C1 = 21293946.348;
    sats[1].C1 = 22744995.395;
    sats[2].C1 = 21565628.372;
    sats[3].C1 = 21513278.753;
    sats[4].C1 = 21488700.402;
    sats[5].C1 = 21246116.669;

    //доплер [Гц]
    sats[0].D1 = 2677.134;
    sats[1].D1 = 5331.631;
    sats[2].D1 = -484.137;
    sats[3].D1 = 3400.669;
    sats[4].D1 = -347.145;
    sats[5].D1 = 143.621;

    //сдвиг часов [мкс]
    sats[0].dT = -50.014071;
    sats[1].dT = -148.368549;
    sats[2].dT = 526.543696;
    sats[3].dT = 53.373334;
    sats[4].dT = 78.206308;
    sats[5].dT = -409.011384;

    //последовательные приближения
    for (quint8 it = 0; it < 8; it++) {

        TMatrix H (6,4);

        vec3 r[6];
        for (int i = 0; i < 6; i++) {
            r[i].xyz = sats[i].xyz - P;
            for (int j = 0; j < 3; j++) {
                H(i, j) = - r[i].xyz[j] / r[i].xyz.length();
            }
            H(i, 3) = 1;    //поправка в метрах и считать сразу с*deltaT
        }

        //матрица невязок
        TMatrix V(6, 1);

        //вращение приемника вместе с землёй

        if (it == 7) {
            cout << endl << "gg is coming" << endl;
        }

        TVector Ve(3);
        Ve[0] = 0;
        Ve[1] = 0;
        Ve[2] = 7.2921151467 / 100000;
        Ve = Ve ^ P;

        //расстояние приемник-спутник
        vec3 R[6];

        //нормаль на спутник
        vec3 e[6];

        //скорость изменения дальности
        TVector Vr(6);

        for (int i = 0; i < 6; i++) {
            sats[i].P1 = sats[i].C1 + sats[i].dT*(C / 1000000); // - Rdcv
            sats[i].dt = sats[i].P1 / C;

            R[i].xyz = sats[i].xyz - P;
            e[i].xyz = R[i].xyz.norm();

            Vr[i] = sats[i].D1*0.1902 - Ve*e[i].xyz;

            V(i, 0) = sats[i].C1 - r[i].xyz.length() + sats[i].dT*(C/1000000) - (sats[i].D1*0.1902)*sats[i].dt;
            //V(i, 0) = sats[i].P1 - R[i].xyz.length() - (sats[i].D1*0.1902)*sats[i].dt;
            //V(i, 0) = sats[i].P1 - (R[i].xyz.length() + Vr[i]*sats[i].dt);
        }

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
