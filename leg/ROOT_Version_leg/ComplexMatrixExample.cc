#include <Eigen/Dense>
#include <complex>
#include <iostream>

using Complex = std::complex<double>;
using CMatrix2 = Eigen::Matrix2cd;
using CMatrix4 = Eigen::Matrix4cd;
using CVector4 = Eigen::Vector4cd;

void ComplexMatrixExample()
{
    const Complex I(0.0, 1.0);

    // Example 2x2 spin-amplitude matrix:
    //
    //       initial proton
    //         up      down
    // Xi up   M00      M01
    // Xi down M10      M11

    CMatrix2 M;

    M << Complex(1.0,  0.2), Complex(0.3, -0.1),
         Complex(0.1,  0.4), Complex(0.8,  0.3);

    M *= Complex(0.,1.0); // Multiply by i
    // Pauli sigma_y
    CMatrix2 sigmaY;

    sigmaY << Complex(0.0, 0.0), -I,
              I,                  Complex(0.0, 0.0);

    // rho = M M^\dagger
    CMatrix2 rho = M * M.adjoint();

    double denominator = rho.trace().real();

    // Average over the initial proton spin
    double spinSum = 0.5 * denominator;

    // Xi polarization along y
    double Py =
        (rho * sigmaY).trace().real() / denominator;

    std::cout << "M =\n" << M << "\n\n";
    std::cout << "M dagger =\n" << M.adjoint() << "\n\n";
    std::cout << "rho = M M dagger =\n" << rho << "\n\n";
    std::cout << "spin sum = " << spinSum << "\n";
    std::cout << "Py = " << Py << "\n";
}
