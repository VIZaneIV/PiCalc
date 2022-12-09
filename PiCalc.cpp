/***************************************************************
 * Computing pi by Binary Splitting Algorithm with GMP libarary.
 **************************************************************/
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include "mpir/mpirxx.h"

 // C++ program for implementation of KMP pattern searching
 // algorithm

using namespace std;

struct PQT
{
    mpz_class P, Q, T;
};

class Chudnovsky
{
    // Declaration
    mpz_class A, B, C, D, E, C3_24;  // GMP Integer
    int DIGITS, PREC, N;             // Integer
    double DIGITS_PER_TERM;          // Long
    clock_t t0, t1, t2;              // Time
    PQT compPQT(int n1, int n2);     // Computer PQT (by BSA)

public:
    Chudnovsky(long);                // Constructor
    void compPi();                   // Compute PI
};

/*
 * Constructor
 */
Chudnovsky::Chudnovsky(long _digits=100000000)
{
    // Constants
    DIGITS = _digits;
    A = 13591409;
    B = 545140134;
    C = 640320;
    D = 426880;
    E = 10005;
    DIGITS_PER_TERM = 14.1816474627254776555;  // = log(53360^3) / log(10)
    C3_24 = C * C * C / 24;
    N = DIGITS / DIGITS_PER_TERM;
    PREC = DIGITS * log2(10);
}

/*
 * Compute PQT (by Binary Splitting Algorithm)
 */
PQT Chudnovsky::compPQT(int n1, int n2)
{
    int m;
    PQT res;

    if (n1 + 1 == n2) {
        res.P = (2 * n2 - 1);
        res.P *= (6 * n2 - 1);
        res.P *= (6 * n2 - 5);
        res.Q = C3_24 * n2 * n2 * n2;
        res.T = (A + B * n2) * res.P;
        if ((n2 & 1) == 1) res.T = -res.T;
    }
    else {
        m = (n1 + n2) / 2;
        PQT res1 = compPQT(n1, m);
        PQT res2 = compPQT(m, n2);
        res.P = res1.P * res2.P;
        res.Q = res1.Q * res2.Q;
        res.T = res1.T * res2.Q + res1.P * res2.T;
    }

    return res;
}

/*
 * Compute PI
 */
void Chudnovsky::compPi()
{
    cout << "**** PI Computation ( " << DIGITS << " digits )" << endl;

    // Time (start)
    t0 = clock();

    // Compute Pi
    PQT PQT = compPQT(0, N);
    mpf_class pi(0, PREC);
    pi = D * sqrt((mpf_class)E) * PQT.Q;
    pi /= (A * PQT.Q + PQT.T);

    // Time (end of computation)
    t1 = clock();
    cout << "TIME (COMPUTE): "
        << (double)(t1 - t0) / CLOCKS_PER_SEC
        << " seconds." << endl;

    // Output
    ofstream ofs("pi.txt");
    ofs.precision(DIGITS + 1);
    ofs << pi << endl;

    // Time (end of writing)
    t2 = clock();
    cout << "TIME (WRITE)  : "
        << (double)(t2 - t1) / CLOCKS_PER_SEC
        << " seconds." << endl;
}

void computeLPSArray(char* pat, int M, int* lps);

// Prints occurrences of txt[] in pat[]
void KMPSearch(char* pat, char* txt)
{
    const int M = strlen(pat);
    const int N = strlen(txt);

    // create lps[] that will hold the longest prefix suffix
    // values for pattern
    //int lps[M]; 
    int* lps = new int[M];


    // Preprocess the pattern (calculate lps[] array)
    computeLPSArray(pat, M, lps);

    int i = 0; // index for txt[]
    int j = 0; // index for pat[]
    while ((N - i) >= (M - j)) {
        if (pat[j] == txt[i]) {
            j++;
            i++;
        }

        if (j == M) {
            printf("Found pattern at index %d ", i - j);
            j = lps[j - 1];
        }

        // mismatch after j matches
        else if (i < N && pat[j] != txt[i]) {
            // Do not match lps[0..lps[j-1]] characters,
            // they will match anyway
            if (j != 0)
                j = lps[j - 1];
            else
                i = i + 1;
        }
    }
    delete[] lps;
}

// Fills lps[] for given pattern pat[0..M-1]
void computeLPSArray(char* pat, int M, int* lps)
{
    // length of the previous longest prefix suffix
    int len = 0;

    lps[0] = 0; // lps[0] is always 0

    // the loop calculates lps[i] for i = 1 to M-1
    int i = 1;
    while (i < M) {
        if (pat[i] == pat[len]) {
            len++;
            lps[i] = len;
            i++;
        }
        else // (pat[i] != pat[len])
        {
            // This is tricky. Consider the example.
            // AAACAAAA and i = 7. The idea is similar
            // to search step.
            if (len != 0) {
                len = lps[len - 1];

                // Also, note that we do not increment
                // i here
            }
            else // if (len == 0)
            {
                lps[i] = 0;
                i++;
            }
        }
    }
}

void loadFiles() {
    fstream File;

    string pi;
    string line, substring;

    // Loading previously calculated Pi
    File.open("pi.txt", ios::in);
    if (File.is_open()) {
        getline(File, pi);
    }
    else {
        cout << "Failed to open the file with Pi" << endl;
        return;
    }

    int n = pi.length();
    char* piChar = new char[n + 1];
    for (int i = 0; i < sizeof(piChar); i++) {
        piChar[i] = pi[i];
    }

    File.close();

    //Loading and searching for substrings
    File.open("table.txt", ios::in);
    if (File.is_open()) {
        char* ptr;
        while (getline(File, line)) {
            cout << line << endl;
            std::size_t pos = line.find(',');
            if (pos != std::string::npos) {
                substring = line.substr(pos + 1);
                int m = substring.length();
                char* substringChar = new char[m + 1];

                for (int i = 0; i < sizeof(substring); i++) {
                    substringChar[i] = substring[i];
                }
                KMPSearch(piChar, substringChar);
                delete[] substringChar;
            }
        }
    }
    File.close();
    delete[] piChar;
}

int main()
{
    int power;
    cout << "Provide precision 10^x: ";
    cin >> power;
    try
    {
        // Instantiation
        Chudnovsky objMain(pow(10,power));

        // Compute PI
        objMain.compPi();
    }
    catch (...) {
        cout << "ERROR!" << endl;
        return -1;
    }

    //loadFiles();

    return 0;
}