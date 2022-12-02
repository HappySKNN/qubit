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

    static std::string Get_Arithmetic_Table(const int &a, const int &m);

    static std::vector<int> Transpose_Matrix(const std::vector<int> &vec);

    static std::vector<int> Dense_Matrices_Multiply(const std::vector<int> &matrix_1, const std::vector<int> &matrix_2);

    static std::vector<std::vector<int>>
    Multiply_Matrices(std::vector<std::vector<int>> matrix_1, std::vector<std::vector<int>> matrix_2);

    static void Get_Qubit(const std::string &filename);

    static void Get_Qubit(const std::vector<int> &boolean_vector);

    static std::vector<std::vector<int>> Get_Vector_From_Matrix(const std::string &filename);

    static std::string Convert_Decimal_To_Binary(const int &number, const int &length);

    static std::vector<std::vector<int>> Tensor_Multiply(const std::vector<std::vector<int>> &matrix_1,
                                                         const std::vector<std::vector<int>> &matrix_2);

    static void Deutsch_Jozsa_Algorithm(const std::string &filename);
};

#endif //QUBITS_QUBIT_H
