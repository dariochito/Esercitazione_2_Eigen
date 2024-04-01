#include <iostream>
#include "Eigen/Eigen"
#include <iomanip>

using namespace std;
using namespace Eigen;


//The function SolutionPalu resolves the linear system with PALU decomposition and prints the solution and the relative error
Vector2d SolutionPalu (Matrix2d& A,Vector2d& b)
{
    Vector2d x_palu = A.fullPivLu().solve(b);
    cout << "The solution x with PALU decomposition is:" << endl;
    cout << "x = [";
    for(int i = 0; i < x_palu.size(); ++i)
    {
        cout << fixed << setprecision(16) << scientific << x_palu[i] << "," ;
    }
    cout << "]" << endl;

    Vector2d exact_solution = -Vector2d::Ones();

    double errRelPalu = (exact_solution - x_palu).norm() / exact_solution.norm();
    cout << fixed << setprecision(16) << scientific << "errRel_palu = " << errRelPalu << endl;

    return x_palu;

}

//The function SolutionPalu resolves the linear system with QR decomposition and prints the solution and the relative error
Vector2d SolutionQr(Matrix2d& A,Vector2d& b)
{
    Vector2d x_qr = A.householderQr().solve(b);
    cout << "The solution x with QR decomposition is:" << endl;
    cout << "x = [";
    for(int i = 0; i < x_qr.size(); ++i)
    {
        cout << fixed << setprecision(16) << scientific << x_qr[i] << "," ;
    }
    cout << "]" << endl;

    Vector2d exact_solution = -Vector2d::Ones();
    double errRelQr = (exact_solution - x_qr).norm() / exact_solution.norm();
    cout << fixed << setprecision(16) << scientific << "errRel_qr = " << errRelQr << endl;

    return x_qr;
}

int main()
{
    cout << "First linear system" << endl;
    Matrix2d A1;
    A1 << 5.547001962252291e-01, -3.770900990025203e-02, 8.320502943378437e-01, -9.992887623566787e-01;
    Vector2d b1;
    b1 << -5.169911863249772e-01, 1.672384680188350e-01;
    Vector2d x1_palu = SolutionPalu(A1,b1);
    Vector2d x1_qr = SolutionQr(A1,b1);

    cout << "\n";

    cout << "Second linear system" << endl;
    Matrix2d A2;
    A2 << 5.547001962252291e-01, -5.540607316466765e-01, 8.320502943378437e-01, -8.324762492991313e-01;
    Vector2d b2;
    b2 << -6.394645785530173e-04, 4.259549612877223e-04;
    Vector2d x2_palu = SolutionPalu(A2,b2);
    Vector2d x2_qr = SolutionQr(A2,b2);

    cout << "\n";

    cout << "Third linear system" << endl;
    Matrix2d A3;
    A3 << 5.547001962252291e-01, -5.547001955851905e-01, 8.320502943378437e-01,-8.320502947645361e-01;
    Vector2d b3;
    b3 << -6.400391328043042e-10, 4.266924591433963e-10;
    Vector2d x3_palu = SolutionPalu(A3,b3);
    Vector2d x3_qr = SolutionQr(A3,b3);


  return 0;
}
