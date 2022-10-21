#ifndef QUBITS_QUBIT_H
#define QUBITS_QUBIT_H

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

#define EPS 1e-9

class Qubit {
private:
    double a_0, a_1, p_0{}, p_1{};

    [[nodiscard]] bool Check_Correct() const;

    void Calculate_Chances();

    void Print_Qubit() const;

public:
    Qubit();

    Qubit(double a_0, double a_1);

    void Hadamard_Multiply();

    [[nodiscard]] double Qubit_Observe() const;

    static void Solve_Function(const int &a, const int &m);

    static void Solve_Function(const std::string &filename);

    static void Get_Matrix(const std::string &filename);

    static void Calculate_Function(int a, int m);

    static std::string Get_Arithmetic_Table(const int &a, const int &m);

    static std::vector<int> Transpose_Matrix(const std::vector<int> &vec);

    static std::vector<int> Multiply_Matrices(const std::vector<int> &vec, const std::vector<int> &vec1);
};

#endif //QUBITS_QUBIT_H
