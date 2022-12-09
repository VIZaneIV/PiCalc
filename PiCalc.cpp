#include <cmath>
#include <iostream>
#include <fstream>
#include "mpir/mpirxx.h"

using namespace std;

struct PQT
{
    mpz_class P, Q, T;
};

class Chudnovsky
{
  public:
    mpz_class A, B, C, D, E, C3_24;
    int DIGITS, PREC, N;
    double DIGITS_PER_TERM;
    clock_t t0, t1;
    PQT compPQT(int n1, int n2);     // Computer PQT (by BSA)
    mpf_class PI;

  public:
    Chudnovsky(int digitsNum);
    void compPi();
};

Chudnovsky::Chudnovsky(int digitsNum)
{
    DIGITS = digitsNum;
    A      = 13591409;
    B      = 545140134;
    C      = 640320;
    D      = 426880;
    E      = 10005;
    DIGITS_PER_TERM = 14.1816474627254776555;  // = log(53360^3) / log(10)
    C3_24  = C * C * C / 24;
    N      = DIGITS / DIGITS_PER_TERM;
    PREC   = DIGITS * log2(10);
}

/*
 * Compute PQT (by Binary Splitting Algorithm)
 */
PQT Chudnovsky::compPQT(int n1, int n2)
{
    int m;
    PQT res;

    if (n1 + 1 == n2) {
        res.P  = (2 * n2 - 1);
        res.P *= (6 * n2 - 1);
        res.P *= (6 * n2 - 5);
        res.Q  = C3_24 * n2 * n2 * n2;
        res.T  = (A + B * n2) * res.P;
        if ((n2 & 1) == 1) res.T = - res.T;
    } else {
        m = (n1 + n2) / 2;
        PQT res1 = compPQT(n1, m);
        PQT res2 = compPQT(m, n2);
        res.P = res1.P * res2.P;
        res.Q = res1.Q * res2.Q;
        res.T = res1.T * res2.Q + res1.P * res2.T;
    }

    return res;
}

void Chudnovsky::compPi()
{
    cout << "**** PI Computation ( " << DIGITS << " digits )" << endl;

    // Time (start)
    t0 = clock();

    // Compute Pi
    PQT PQT = compPQT(0, N);
    mpf_class pi(0, PREC);
    pi  = D * sqrt((mpf_class)E) * PQT.Q;
    pi /= (A * PQT.Q + PQT.T);
    PI = pi;

    // Time (end of computation)
    t1 = clock();
    cout << "TIME (COMPUTE): "
         << (double)(t1 - t0) / CLOCKS_PER_SEC
         << " seconds." << endl;
}

void computeLPSArray(const char* pat, int M, int* lps)
{
    int len = 0;
    lps[0] = 0;

    int i = 1;
    while (i < M) {
        if (pat[i] == pat[len]) {
            len++;
            lps[i] = len;
            i++;
        }
        else
        {
            if (len != 0) {
                len = lps[len - 1];
            }
            else
            {
                lps[i] = 0;
                i++;
            }
        }
    }
}

int KMPSearch(const char* pat, char* txt)
{
    int M = strlen(pat);
    int N = strlen(txt);

    int lps[M];

    computeLPSArray(pat, M, lps);

    int i = 0;
    int j = 0;
    while ((N - i) >= (M - j)) {
        if (pat[j] == txt[i]) {
            j++;
            i++;
        }

        if (j == M) {
            return i - j;
        }

        else if (i < N && pat[j] != txt[i]) {
            if (j != 0)
                j = lps[j - 1];
            else
                i = i + 1;
        }
    }
  return -1;
}

void writePi(Chudnovsky obj){
  // Writing to a file
    clock_t t0 = clock();
    ofstream ofs ("pi.txt");
    ofs.precision(obj.DIGITS + 1);
    ofs << obj.PI << endl;

    // Time (end of writing)
    clock_t t1 = clock();
    cout << "TIME (WRITE)  : "
         << (double)(t1 - t0) / CLOCKS_PER_SEC
         << " seconds." << endl;
}

char* readPi(string filename, int start, int end){
  ifstream file (filename);
  if(file.is_open())
  {
    file.seekg(start);
    char *s = new char[end - start + 1];
    file.read(s, end - start);
    return s;
  }
  return nullptr;
}

void f_pi(mpf_class pi, int start, int end){
  char* piSubStr = readPi("pi.txt", start, end);

  ofstream tablefile;
  tablefile.open("table.txt");
  for (int i = 0; i < end - start; i++){
    tablefile << i + start << "," << KMPSearch(to_string(i).c_str(), piSubStr) << endl;
  }
  tablefile.close();
}

int main()
{
  int n = 0;
  cout << "Type number of digits you want to compute:\n";
  cin >> n;

    try
    {
        Chudnovsky objMain(n);
        objMain.compPi();
        writePi(objMain);
        f_pi(objMain.PI, 10, 30);
    }
    catch (...) {
        cout << "ERROR!" << endl;
        return -1;
    }

    return 0;
}
