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

    static void Get_Matrix(const std::string &filename);

    static std::vector<int> Transpose_Matrix(const std::vector<int> &vec);

    static std::vector<int> Multiply_Matrix(const std::vector<int> &vec);
};

#endif //QUBITS_QUBIT_H
