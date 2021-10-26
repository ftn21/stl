#include <iostream>
#include "linalg.h"
#include <cstring>
#include <math.h>
#include <fstream>
#include <bitset>
#include <iomanip>
#include <iterator>
#include <vector>
#include <cmath>

// константы
#define pi 3.1415926535897932
#define GM 398600500000000    //GM = 3.986005*10^14 [m^3 /s^2]
#define omega_e 0.00007292115 //value of the Earth’s rotation rate
#define RTCM3PREAMB 211
#define light_speed_m 299792.458

using namespace std;
using namespace MyLinearAlgebra;

//-----------------------------------------------------------------------------------------------------------------------------------------------
// структура для хранения пакета rtcm1019
typedef struct RTCM1019_STRUCTURE
{
    /* header | 24  bits */
    int32_t preamble; //      |  8  bits
    //6 zeroes        //      |  6  bits
    int32_t length; //      | 10  bits

    /* data message | 488  bits */
    int32_t type;         // DF002 | 12  bits | u
    int32_t satellite_id; // DF009 |  6  bits | u
    int32_t week;         // DF076 | 10  bits | u
    int32_t sv_accuracy;  // DF077 |  4  bits | u
    int32_t code_on_l2;   // DF078 |  2  bits | u
    double idot;          // DF079 | 14  bits | s
    int32_t iode;         // DF071 |  8  bits | u
    double t_oc;          // DF081 | 16  bits | u
    double a_f2;          // DF082 |  8  bits | s
    double a_f1;          // DF083 | 16  bits | s
    double a_f0;          // DF084 | 22  bits | s
    int32_t idoc;         // DF085 | 10  bits | u
    double c_rs;          // DF086 | 16  bits | s
    double delta_n;       // DF087 | 16  bits | s
    double M0;            // DF088 | 32  bits | s
    double c_uc;          // DF089 | 16  bits | s
    double e;             // DF090 | 32  bits | u
    double c_us;          // DF091 | 16  bits | s
    double A_sqrt;        // DF092 | 32  bits | u
    double t_oe;          // DF093 | 16  bits | u | время действия эфемирид
    double c_ic;          // DF094 | 16  bits | s
    double OMEGA0;        // DF095 | 16  bits | s
    double c_is;          // DF096 | 16  bits | s
    double i0;            // DF097 | 32  bits | s
    double c_rc;          // DF098 | 16  bits | s
    double omega;         // DF099 | 32  bits | s
    double OMEGADOT;      // DF100 | 24  bits | s
    double t_GD;          // DF101 |  8  bits | s
    int32_t sv_health;    // DF102 |  6  bits | u
    bool flag;            // DF103 |  1  bits | u
    char fit;             // DF137 |  1  bits | u

    /* parity | 24  bits */
    int32_t parity; //      | 24  bits
} rtcm1019_pack;

// структура для хранения пакета rtcm1002
typedef struct RTCM1002_STRUCTURE
{
    /* data */
    uint8_t satellite_id;    // DF009 |  6  bits | u
    uint32_t L1_pseudorange; // DF011 | 24  bits | u
    uint32_t L1_ambiguity;   // DF014 |  8  bits | u

    double C1; // L1_pseudorange + L1_ambiguity
} rtcm1002_pack;

// структура для координат спутника + его id 
typedef struct XYZ_COORD
{
    uint8_t sat_id;
    double x;
    double y;
    double z;
} xyz_coord;

// для crc-проверки целостности данных
const uint32_t tbl_CRC24Q[] = {
    0x000000, 0x864CFB, 0x8AD50D, 0x0C99F6, 0x93E6E1, 0x15AA1A, 0x1933EC, 0x9F7F17,
    0xA18139, 0x27CDC2, 0x2B5434, 0xAD18CF, 0x3267D8, 0xB42B23, 0xB8B2D5, 0x3EFE2E,
    0xC54E89, 0x430272, 0x4F9B84, 0xC9D77F, 0x56A868, 0xD0E493, 0xDC7D65, 0x5A319E,
    0x64CFB0, 0xE2834B, 0xEE1ABD, 0x685646, 0xF72951, 0x7165AA, 0x7DFC5C, 0xFBB0A7,
    0x0CD1E9, 0x8A9D12, 0x8604E4, 0x00481F, 0x9F3708, 0x197BF3, 0x15E205, 0x93AEFE,
    0xAD50D0, 0x2B1C2B, 0x2785DD, 0xA1C926, 0x3EB631, 0xB8FACA, 0xB4633C, 0x322FC7,
    0xC99F60, 0x4FD39B, 0x434A6D, 0xC50696, 0x5A7981, 0xDC357A, 0xD0AC8C, 0x56E077,
    0x681E59, 0xEE52A2, 0xE2CB54, 0x6487AF, 0xFBF8B8, 0x7DB443, 0x712DB5, 0xF7614E,
    0x19A3D2, 0x9FEF29, 0x9376DF, 0x153A24, 0x8A4533, 0x0C09C8, 0x00903E, 0x86DCC5,
    0xB822EB, 0x3E6E10, 0x32F7E6, 0xB4BB1D, 0x2BC40A, 0xAD88F1, 0xA11107, 0x275DFC,
    0xDCED5B, 0x5AA1A0, 0x563856, 0xD074AD, 0x4F0BBA, 0xC94741, 0xC5DEB7, 0x43924C,
    0x7D6C62, 0xFB2099, 0xF7B96F, 0x71F594, 0xEE8A83, 0x68C678, 0x645F8E, 0xE21375,
    0x15723B, 0x933EC0, 0x9FA736, 0x19EBCD, 0x8694DA, 0x00D821, 0x0C41D7, 0x8A0D2C,
    0xB4F302, 0x32BFF9, 0x3E260F, 0xB86AF4, 0x2715E3, 0xA15918, 0xADC0EE, 0x2B8C15,
    0xD03CB2, 0x567049, 0x5AE9BF, 0xDCA544, 0x43DA53, 0xC596A8, 0xC90F5E, 0x4F43A5,
    0x71BD8B, 0xF7F170, 0xFB6886, 0x7D247D, 0xE25B6A, 0x641791, 0x688E67, 0xEEC29C,
    0x3347A4, 0xB50B5F, 0xB992A9, 0x3FDE52, 0xA0A145, 0x26EDBE, 0x2A7448, 0xAC38B3,
    0x92C69D, 0x148A66, 0x181390, 0x9E5F6B, 0x01207C, 0x876C87, 0x8BF571, 0x0DB98A,
    0xF6092D, 0x7045D6, 0x7CDC20, 0xFA90DB, 0x65EFCC, 0xE3A337, 0xEF3AC1, 0x69763A,
    0x578814, 0xD1C4EF, 0xDD5D19, 0x5B11E2, 0xC46EF5, 0x42220E, 0x4EBBF8, 0xC8F703,
    0x3F964D, 0xB9DAB6, 0xB54340, 0x330FBB, 0xAC70AC, 0x2A3C57, 0x26A5A1, 0xA0E95A,
    0x9E1774, 0x185B8F, 0x14C279, 0x928E82, 0x0DF195, 0x8BBD6E, 0x872498, 0x016863,
    0xFAD8C4, 0x7C943F, 0x700DC9, 0xF64132, 0x693E25, 0xEF72DE, 0xE3EB28, 0x65A7D3,
    0x5B59FD, 0xDD1506, 0xD18CF0, 0x57C00B, 0xC8BF1C, 0x4EF3E7, 0x426A11, 0xC426EA,
    0x2AE476, 0xACA88D, 0xA0317B, 0x267D80, 0xB90297, 0x3F4E6C, 0x33D79A, 0xB59B61,
    0x8B654F, 0x0D29B4, 0x01B042, 0x87FCB9, 0x1883AE, 0x9ECF55, 0x9256A3, 0x141A58,
    0xEFAAFF, 0x69E604, 0x657FF2, 0xE33309, 0x7C4C1E, 0xFA00E5, 0xF69913, 0x70D5E8,
    0x4E2BC6, 0xC8673D, 0xC4FECB, 0x42B230, 0xDDCD27, 0x5B81DC, 0x57182A, 0xD154D1,
    0x26359F, 0xA07964, 0xACE092, 0x2AAC69, 0xB5D37E, 0x339F85, 0x3F0673, 0xB94A88,
    0x87B4A6, 0x01F85D, 0x0D61AB, 0x8B2D50, 0x145247, 0x921EBC, 0x9E874A, 0x18CBB1,
    0xE37B16, 0x6537ED, 0x69AE1B, 0xEFE2E0, 0x709DF7, 0xF6D10C, 0xFA48FA, 0x7C0401,
    0x42FA2F, 0xC4B6D4, 0xC82F22, 0x4E63D9, 0xD11CCE, 0x575035, 0x5BC9C3, 0xDD8538};

uint32_t crc24q(const uint8_t *buff, int32_t len)
{
    uint32_t crc = 0;
    int i;
    for (i = 0; i < len; i++)
    {
        uint32_t crc0 = (crc << 8) & 0xFFFFFF;
        uint32_t crc1 = crc >> 16;
        crc = crc0 ^ tbl_CRC24Q[crc1 ^ buff[i]];
    }
    return crc;
}

// считывание беззнакового числа
uint32_t bit32u(const uint8_t *buff, int32_t pos, int32_t len)
{
    uint32_t bits = 0;
    for (int32_t i = pos; i < pos + len; i++) {
        bits = (bits << 1) + ((buff[i / 8] >> (7 - i % 8)) & 1u);
    }
    return bits;
}

// считывания числа со знаком
int32_t bit32s(const uint8_t *buff, int32_t pos, int32_t len)
{
    uint32_t bits = bit32u(buff, pos, len);
    if (len <= 0 || 32 <= len || !(bits & (1u << (len - 1)))) {
        return int32_t(bits);
    }
    return int(bits | (~0u << len)); /* extend sign */
}

// проверка на то, является ли пакет rtcm-пакетом
bool is_rtcm3(const uint8_t *data, int32_t &data_size, uint32_t &tp)
{

    if (*data != RTCM3PREAMB || data_size < 3 + 1 + 3)
        return false;

    int32_t len = bit32s(data, 14, 10);

    if (len < 0 || len > 5000 || data_size < 3 + len + 3)
        return false;

    uint32_t par = bit32u(data, 24 + len * 8, 24);
    uint32_t par_chk = crc24q(data, 3 + len);

    if (par_chk != par)
        return false;

    tp = bit32u(data, 24, 12);
    data_size = 3 + len + 3;

    return tp > 0;
}

// считывание rtcm-пакета
void read_rtcmPack(const uint8_t *data, rtcm1019_pack *pack);

// crc-проверка целостности данных
bool crc_chek(const uint8_t *data, rtcm1019_pack *pack)
{
    int32_t crc_pos = 24 + 8 * pack->length;
    pack->parity = bit32u(&data[0], crc_pos, 24); //crc in msg
    uint32_t crc = crc24q(&data[0], 64);          //calc crc
    if (crc = pack->parity)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------

// структура для одного спутника
struct sattelite
{
public:
    int prn;
    TVector ecef; //координаты в ECEF
    sattelite() : ecef(3){};
    double C1; //псевдодальность [м]
    double dt; //оценка времени полета сигнала
    double dT; //сдвиг часов [мкс]
    double P1; //псевдодальность с коррекцией часов
};

// вектор координат 
class vec3
{
public:
    TVector xyz;
    vec3() : xyz(3){};
};

const double C = 299792458;  // скорость света [м/с]

TVector P(3);  // вектор координат точки, которую мы аппроксимируем последовательными приближениями

sattelite *sats;  // динамический массив спутников

// функция рассчета координат в ECEF
xyz_coord algorithm(rtcm1019_pack pack19, rtcm1002_pack pack02, double tow, int st_num, sattelite *sats_arr);

int main(int argc, char *argv[])
{

    /* 1019 */
    // считвание файла с пакетом rtcm1019
    string fname = "/home/ftn21/Documents/MAI/6sem/app_sredstva/sources/16.04.2021/1019.rtcm";
    ifstream rtcm_st6;
    rtcm_st6.open(fname, ios::binary);
    int rtcm19_filelength = 0;

    // проверка на корректное открытие файла
    if (!rtcm_st6.is_open())
    {
        cout << "error opening " << fname << endl;
    }
    else
    {
        // получить размер файла:
        rtcm_st6.seekg(0, rtcm_st6.end);
        rtcm19_filelength = rtcm_st6.tellg();
        rtcm_st6.seekg(0, rtcm_st6.beg);
        cout << "fileLength of the " << fname << endl
             << " is " << rtcm19_filelength << endl
             << endl;
    }

    // считывание данных из файла
    vector<uint8_t> buf(istreambuf_iterator<char>(rtcm_st6), {});

    int32_t s = 67;
    uint32_t c = 1019;

    int rtcm19_num = 0; //кол-во пакетов rtcm19

    //поиск пакетов rtcm19
    for (int j = 0; j < (rtcm19_filelength - 66); j++)
    {
        if (is_rtcm3(&buf[j], s, c))
        {
            rtcm19_num++;
            cout << "the " << rtcm19_num << " package is rtcm. \n";
        }
    }

    rtcm1019_pack *rtcm19_msg = new rtcm1019_pack[rtcm19_num];

    //чтение пакетов rtcm19
    int rtcm19_count = 0;
    for (int j = 0; j < (rtcm19_filelength - 66); j++)
    {
        if (is_rtcm3(&buf[j], s, c))
        {
            rtcm19_count++;
            cout << "reading the " << rtcm19_count << " package. \n";

            // rtcm read
            rtcm1019_pack pack;
            read_rtcmPack(&buf[j], &rtcm19_msg[rtcm19_count - 1]);

            // parity | 24  bits
            if (crc_chek(&buf[j], &rtcm19_msg[rtcm19_count - 1]))
            {
                cout << rtcm19_count << " rtcm 1019 msg is fine. \n\n";
            }
            else
            {
                cout << rtcm19_count << " rtcm 1019 msg is damaged or crc error. \n\n";
            }
        }
    }

    /* 1002 */
    // считвание файла с пакетом rtcm1002
    string fname1002 = "/home/ftn21/Documents/MAI/6sem/app_sredstva/sources/16.04.2021/1002.rtcm";
    ifstream rtcm02_st6;
    rtcm02_st6.open(fname1002, ios::binary);
    int rtcm02_filelength = 0;

    // проверка на корректное открытие файла
    if (!rtcm02_st6.is_open())
    {
        cout << "error opening " << fname << endl;
    }
    else
    {
        // получить размер файла
        rtcm02_st6.seekg(0, rtcm02_st6.end);
        rtcm02_filelength = rtcm02_st6.tellg();
        rtcm02_st6.seekg(0, rtcm02_st6.beg);
        cout << "fileLength is " << rtcm02_filelength << endl
             << endl;
    }

    // чтение данных
    vector<uint8_t> buf02(istreambuf_iterator<char>(rtcm02_st6), {});
    int32_t s02 = 74;
    uint32_t c02 = 1002;
    int i_1002 = 0;
    for (int i = 0; i < (rtcm02_filelength); i++)
    {
        if (is_rtcm3(&buf02[i], s02, c02))
        {
            cout << "the package is rtcm 1002. " << endl
                 << endl;
            i_1002 = i;
        }
        else
        {
            //cout << "the package is NOT rtcm. " << endl
            //<< endl;
        }
    }

    /* заголовок | 24  bits*/
    int32_t preamble02 = bit32u(&buf02[i_1002], 0, 8);
    int32_t length02 = bit32u(&buf02[i_1002], 14, 10);

    /* данные | 64  bits*/
    uint32_t msg_num = bit32u(&buf02[i_1002], 24, 12);
    uint32_t sta_id = bit32u(&buf02[i_1002], 36, 12);
    uint32_t tow = bit32u(&buf02[i_1002], 48, 30) / 1000;
    uint32_t no_stl = bit32u(&buf02[i_1002], 79, 5); //кол-во спутников в 1002 пакете

    int sat_num = 0; //кол-во спутников
    sat_num = no_stl;

    rtcm1002_pack *rtcm02_msg = new rtcm1002_pack[no_stl];

    /* 6 blocks | 6*74  bits*/
    for (int i = 0; i < no_stl; i++)
    {
        rtcm02_msg[i].satellite_id = bit32u(&buf02[0], 88 + i * 74, 6);
        rtcm02_msg[i].L1_pseudorange = bit32u(&buf02[0], 95 + i * 74, 24);
        rtcm02_msg[i].L1_ambiguity = bit32u(&buf02[0], 146 + i * 74, 8);
        rtcm02_msg[i].C1 = rtcm02_msg[i].L1_pseudorange * 0.02 + rtcm02_msg[i].L1_ambiguity * light_speed_m;
    }

    sats = new sattelite[sat_num]; //спутники
    for (int i = 0; i < sat_num; i++)
    {
        sats[i].prn = rtcm02_msg[i].satellite_id;
        sats[i].C1 = rtcm02_msg[i].C1;
    }

    xyz_coord res[sat_num];

    for (int i = 0; i < rtcm19_num; i++)
    {
        for (int j = 0; j < no_stl; j++)
        {
            if (rtcm19_msg[i].satellite_id == rtcm02_msg[j].satellite_id)
            {
                res[j] = algorithm(rtcm19_msg[i], rtcm02_msg[j], tow, sat_num, sats);
            }
        }
    }

    //ecef
    for (int i = 0; i < sat_num; i++)
    {
        for (int j = 0; j < sat_num; j++)
        {
            if (sats[i].prn == res[j].sat_id)
            {
                sats[i].ecef[0] = res[j].x;
                sats[i].ecef[1] = res[j].y;
                sats[i].ecef[2] = res[j].z;
                cout << sats[i].prn << "  " << float(sats[i].ecef[0]) << "  " << float(sats[i].ecef[1]) << "  " << float(sats[i].ecef[2]) << endl ;
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------------

    // координаты точки приёмника, которые ищем методом последовательных приближений
    P[0] = 0;
    P[1] = 0;
    P[2] = 0;

    //поправка по времени
    double Rdcv = 0;

    sattelite sats_static[sat_num]; 
    for (int i=0; i<sat_num; i++) {
        sats_static[i] = sats[i];
        cout << i << "  " << "ecef: " << sats_static[i].ecef[0] << "  " << sats_static[i].ecef[1] << sats_static[i].ecef[2] << "  ";
    }
    // очистка памяти
    delete [] sats;

    delete [] rtcm19_msg;
    delete [] rtcm02_msg;

    // последовательные приближения
    for (uint8_t it = 0; it < 8; it++)
    {

        TMatrix H(sat_num, 4);

        cout << endl;
        vec3 r[sat_num];
        for (int i = 0; i < sat_num; i++)
        {
            r[i].xyz = sats_static[i].ecef - P;
            for (int j = 0; j < 3; j++)
            {
                H(i, j) = -r[i].xyz[j] / r[i].xyz.length();
            }
            H(i, 3) = 1; // поправка в метрах и считать сразу с*deltaT
        }

        // матрица невязок
        TMatrix V(sat_num, 1);

        // вращение приемника вместе с землёй
        TVector V_e(3);
        TVector Ve(3);
        V_e[0] = 0;
        V_e[1] = 0;
        V_e[2] = 7.2921151467 / 100000;
        Ve = V_e ^ P;

        // расстояние приемник-спутник
        vec3 R[sat_num];

        // нормаль на спутник
        vec3 e[sat_num];

        // скорость изменения дальности
        TVector Vr(sat_num);

        for (int i = 0; i < sat_num; i++)
        {
            sats_static[i].P1 = sats_static[i].C1 + sats_static[i].dT * (C / 1000000) - Rdcv;
            sats_static[i].dt = sats_static[i].P1 / C;   // время полета сигнала

            R[i].xyz = sats_static[i].ecef - P;
            cout << i << ":  " << R[i].xyz[0] << "  ,   " << R[i].xyz[1] << "  ,   " << R[i].xyz[2] << "  ;   " << endl;
            //e[i].xyz = R[i].xyz.norm(); //строка - убийца программы

            double b = (Ve * R[i].xyz) / R[i].xyz.length();

            Vr[i] = /*sats_static[i].D1*0.1902*/ 0 - b; // считаем доплера без доплера, не уверена, на сколько корректно

            if (i == sat_num)
                cout << endl;

            V(i, 0) = sats_static[i].P1 - (R[i].xyz.length() + Vr[i] * sats_static[i].dt);
            //cout << "V(" << i << ", 0) = " << V(i, 0) << "    " ;
        }

        //MHK: dV = (H^t * H)^(-1) * H^t * V
        TMatrix Ht = H.t();
        TMatrix HtH = Ht * H;

        TMatrix Q = !HtH;
        TMatrix dV = Q * Ht * V;

        //вычисление поправки
        double dx = dV(0, 0);
        double dy = dV(1, 0);
        double dz = dV(2, 0);
        Rdcv = dV(3, 0);

        //новое приближение
        P[0] += dx;
        P[1] += dy;
        P[2] += dz;

        if (sqrt(dx * dx + dy * dy + dz * dz) < 0.00001)  // если выполняется условие по точности
        {
            break; // прекращаем последовательные приближения
        }
        else
        {
            cout << "calc . . . x" << int(it) + 1 << endl; // номер итерации
        }

        // правильные координаты приёмника. те координаты, с которыми нужно сравнить
        // те, что мы высчитвали последовательными приближениями
        TVector approx_pos (3);
        approx_pos[0] = 2849830.5060;
        approx_pos[1] = 2186822.2813;
        approx_pos[2] = 5252937.0124;

        // как раз сравниваем и считаем ошибку на каждой итерации
        TVector error = P - approx_pos;

        cout << "ошибка = " << error.length() << endl;
        
    }

    //широта (фи)
    double lat = asin(P[2] / P.length());
    //долгота (лямбда)
    double lon = atan2(P[1], P[0]);

    cout << endl
         << "широта: " << lat << endl
         << "долгота: " << lon << endl;

    //матрица перехода из ECEF в ENU
    TMatrix ENU(3, 3);
    ENU(0, 0) = -sin(lon);
    ENU(0, 1) = cos(lon);
    ENU(0, 2) = 0;
    ENU(1, 0) = -sin(lat) * cos(lon);
    ENU(1, 1) = -sin(lat) * sin(lon);
    ENU(1, 2) = cos(lat);
    ENU(2, 0) = cos(lat) * cos(lon);
    ENU(2, 1) = cos(lat) * sin(lon);
    ENU(2, 2) = sin(lat);

    // перевод ECEF->ENU и рассчет азимута и элевации
    cout << endl;
    for (int i = 0; i < sat_num; i++)
    {
        TVector e = (sats_static[i].ecef - P).norm();
        TVector e_enu = ENU * e;
        double az = atan2(e_enu[0], e_enu[1]);
        double el = asin(e_enu[2]);

        cout << sats_static[i].prn << ":   az = " << float(az * 180 / pi) << "    el = " << float(el * 180 / pi) << endl;
    }

    // правильные координаты приёмника. те координаты, с которыми нужно сравнить
    // те, что мы высчитвали последовательными приближениями
    TVector approx_pos(3);
    approx_pos[0] = 2849830.5060;
    approx_pos[1] = 2186822.2813;
    approx_pos[2] = 5252937.0124;

    // как раз сравниваем и считаем итоговую ошибку
    TVector error = P - approx_pos;

    cout << endl
         << "ошибка по x = " << error[0] << endl; //
    cout << "ошибка по y = " << error[1] << endl;
    cout << "ошибка по z = " << error[2] << endl;
    cout << endl
         << "ошибка = " << error.length() << endl;

    //return 0;
}


// реализация функций, объявленных в начале кода
void read_rtcmPack(const uint8_t *data, rtcm1019_pack *pack)
{

    /* header | 24  bits*/
    pack->preamble = bit32u(&data[0], 0, 8);
    pack->length = bit32u(&data[0], 14, 10);

    /* data message | 488  bits*/
    pack->type = bit32u(&data[0], 24, 12);
    pack->satellite_id = bit32u(&data[0], 36, 6);
    pack->week = bit32u(&data[0], 42, 10);
    pack->sv_accuracy = bit32u(&data[0], 52, 4);
    pack->code_on_l2 = bit32u(&data[0], 56, 2);
    pack->idot = pow(2, -43) * bit32s(&data[0], 58, 14) * pi; //4.7066350816749036e-11  /pi ?
    pack->iode = bit32u(&data[0], 72, 8);
    pack->t_oc = pow(2, 4) * bit32u(&data[0], 80, 16);
    pack->a_f2 = bit32s(&data[0], 96, 8);
    pack->a_f1 = pow(2, -43) * bit32s(&data[0], 104, 16);
    pack->a_f0 = pow(2, -31) * bit32s(&data[0], 120, 22);
    pack->idoc = bit32s(&data[0], 142, 10);
    pack->c_rs = pow(2, -5) * bit32s(&data[0], 152, 16);
    pack->delta_n = pow(2, -43) * bit32s(&data[0], 168, 16) * pi;
    pack->M0 = pow(2, -31) * bit32s(&data[0], 184, 32) * pi;
    pack->c_uc = pow(2, -29) * bit32s(&data[0], 216, 16);
    pack->e = pow(2, -33) * bit32u(&data[0], 232, 32);
    pack->c_us = pow(2, -29) * bit32s(&data[0], 264, 16);
    pack->A_sqrt = pow(2, -19) * bit32u(&data[0], 280, 32); //9.8813129168249309e-324
    pack->t_oe = pow(2, 4) * bit32u(&data[0], 312, 16);
    pack->c_ic = pow(2, -29) * bit32s(&data[0], 328, 16);
    pack->OMEGA0 = pow(2, -31) * bit32s(&data[0], 344, 32) * pi;
    pack->c_is = pow(2, -29) * bit32s(&data[0], 376, 16);
    pack->i0 = pow(2, -31) * bit32s(&data[0], 392, 32) * pi; //0.31093034334480762
    pack->c_rc = pow(2, -5) * bit32s(&data[0], 424, 16);
    pack->omega = pow(2, -31) * bit32s(&data[0], 440, 32) * pi;    //0.18036843976005912
    pack->OMEGADOT = pow(2, -43) * bit32s(&data[0], 472, 24) * pi; //-2.6769839678308927e-09
    pack->t_GD = pow(2, -31) * bit32s(&data[0], 496, 8);
    pack->sv_health = bit32u(&data[0], 504, 6);
    pack->flag = bit32u(&data[0], 510, 1);
    pack->fit = bit32u(&data[0], 511, 1); // 0 - curve-fit interval is 4 hours
                                          //1 - curve-fit is greater than 4 hours  как кодировать? 4+0, 4+1 ?
}

xyz_coord algorithm(rtcm1019_pack pack19, rtcm1002_pack pack02, double tow, int st_num, sattelite *sats_arr)
{
    xyz_coord result;
    result.sat_id = pack02.satellite_id;
    //поправка по времени
    double t = tow - pack02.C1 / C;
    double dTs = t - pack19.t_oc;
    for (int i = 0; i < 2; i++)
    {
        dTs -= pack19.a_f0 + pack19.a_f1 * dTs + pack19.a_f2 * dTs * dTs;
    }
    dTs = pack19.a_f0 + pack19.a_f1 * dTs + pack19.a_f2 * dTs * dTs;

    for (int i = 0; i < st_num; i++)
    {
        if (sats_arr[i].prn == result.sat_id)
        {
            sats_arr[i].dT = dTs * 1000000;
        }
    }

    //individual satellite time
    double tbias = t - dTs;

    //time, elapsed since the reference epoch
    double tk = tbias - pack19.t_oe;

    double A = pow(pack19.A_sqrt, 2);
    double n0 = sqrt(GM / pow(A, 3)); //Computed mean motion
    double n = n0 + pack19.delta_n;   //Corrected mean motion

    double M = pack19.M0 + n * tk; //Mean anomaly 0.14682866 rad

    double E = M; //эксцентрическая аномалия
    double Ek = 0;
    for (int i = 0; i < 3; i++)
    {
        Ek = E;
        E = M + pack19.e * sin(Ek);
    }
    Ek = E;

    double v = atan2(sqrt(1.0 - pack19.e * pack19.e) * sin(E), cos(E) - pack19.e); //true anomaly
    double Fi = v + pack19.omega;                                                  // v + omega

    double du = pack19.c_uc * cos(2 * Fi) + pack19.c_us * sin(2 * Fi);
    double u = Fi + du; //argument of latitude

    double dr = pack19.c_rc * cos(2 * Fi) + pack19.c_rs * sin(2 * Fi);
    double r = A * (1 - pack19.e * cos(E)) + dr; //radius-vector

    double di = pack19.c_ic * cos(2 * Fi) + pack19.c_is * sin(2 * Fi);
    double i = pack19.i0 + pack19.idot*tk + di; //inclination

    double X_op = r * cos(u); //X position in the orbital plane
    double Y_op = r * sin(u); //Y position in the orbital plane

    double OMEGA = pack19.OMEGA0 + (pack19.OMEGADOT - omega_e) * tk - omega_e * pack19.t_oe; //longtitude of ascending node [rad]

    //geocentric sattelite coordinates ECEF
    double X_ecef = X_op * cos(OMEGA) - Y_op * sin(OMEGA) * cos(i);
    double Y_ecef = X_op * sin(OMEGA) + Y_op * cos(OMEGA) * cos(i);
    double Z_ecef = Y_op * sin(i);

    result.x = X_ecef;
    result.y = Y_ecef;
    result.z = Z_ecef;

    return result;
}
