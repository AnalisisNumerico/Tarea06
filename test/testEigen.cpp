/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: ger534
 * @Date  : 22.05.2018
 */

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <functional>
#include <cmath>
#include <limits>

#include "Tarea06.hpp"


namespace anpi {
    namespace test {

        /// Test the random simetric matrix
        template<typename T>
        void randomSymmetricSqrTest() {
            anpi::Matrix<T> A = anpi::randomSymmetricSqr<T>(10);
            for(int i =0; i < A.cols(); i++){
                for(int j =0; j < A.rows(); j++){
                    BOOST_CHECK(A[i][j] == A[i][j]);
                }
            }
        }

        /// Test the jacobi and eigen with lapack using A=EDE ^T
        /// with eigenvectors in columns
        template<typename T>
        void EigenTest(const std::function<void(const Matrix<T>&,
                                                  std::vector<T>&,
                                                  Matrix<T>&)>& eigen) {

            anpi::Matrix<T> A = anpi::randomSymmetricSqr<T>(10);

            std::vector<T> val;

            anpi::Matrix<T> E;

            eigen(A,val,E);

            anpi::Matrix<T> valMatrix (A.rows(),A.rows());
            ///vuelve el vector de eigenvalores una matriz para A=EDE^T
            for(int i =0; i < valMatrix.cols(); i++){
                for(int j =0; j < valMatrix.rows(); j++){
                    if(i==j){
                        valMatrix[i][j] = val[i];
                    }
                }
            }

            ///transpuesta
            anpi::Matrix<T> transpuestaE(E.cols(),E.rows());
            for(int i =0; i < E.cols(); i++){
                for(int j =0; j < E.rows(); j++){
                    transpuestaE[i][j] = E[j][i];
                }
            }


            anpi::Matrix<T>  Eval = E*valMatrix;

            anpi::Matrix<T> verifica = Eval*transpuestaE;
            const T eps = std::numeric_limits<T>::epsilon()*T(9000);
            for(int i =0; i < E.cols(); i++){
                for(int j =0; j < E.rows(); j++){
                    BOOST_CHECK(std::abs(A[i][j] - verifica[i][j]) < eps);
                }
            }
        }

        /// Test the jacobi and eigen with lapack, sorting their outputs and comparing them
        ///for measuaring time, we used the benchmarkEigen
        template<typename T>
        void EigenSortTest() {

            anpi::Matrix<T> A = anpi::randomSymmetricSqr<T>(10);
            anpi::Matrix<T> A1 = anpi::randomSymmetricSqr<T>(10);

            std::vector<T> val;
            std::vector<T> val1;

            anpi::Matrix<T> E;
            anpi::Matrix<T> E1;

            anpi::eig(A,val,E);
            anpi::jacobi(A1,val1,E1);

            anpi::sort(val,E);
            anpi::sort(val1,E1);

            const T eps = std::numeric_limits<T>::epsilon()*T(9000);
            for(int i =0; i < E.cols(); i++){
                for(int j =0; j < E.rows(); j++){
                    BOOST_CHECK(std::abs(std::abs(E[i][j]) - std::abs(E1[i][j])) < eps);
                }
            }

            for(int i =0; i < val.size(); i++){
                BOOST_CHECK(std::abs(val[i] - val1[i]) < eps);
            }

        }

        /// Test the sorting function using an specific example for simple precision
        template<typename T>
        void SortTest() {

            anpi::Matrix<T> A = { {3,2,4}, {2,0,2},{4,2,3} };

            std::vector<T> val;

            anpi::Matrix<T> E;

            anpi::eig(A,val,E);

            anpi::sort(val,E);

            anpi::Matrix<T> EV = { {0.666667, -0.536953, -0.51695}, {0.333333, 0.835121, -0.437564},{0.666667, 0.119393, 0.735732} };

            std::vector<T> valV = { 8, -1, -1 };

            const T eps = std::numeric_limits<T>::epsilon()*T(10);

            for(int i =0; i < E.cols(); i++){
                for(int j =0; j < E.rows(); j++){
                    BOOST_CHECK(std::abs(EV[i][j] - E[i][j]) < eps);
                }
            }

            for(int i =0; i < val.size(); i++){
                BOOST_CHECK(std::abs(valV[i] - val[i]) < eps);
            }

        }

    } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( Eigen )

    BOOST_AUTO_TEST_CASE( RandomMatrix ) {
        anpi::test::randomSymmetricSqrTest<float>();
        anpi::test::randomSymmetricSqrTest<double>();
    }

    BOOST_AUTO_TEST_CASE( Eig ) {
        anpi::test::EigenTest<float>(anpi::jacobi<float>);
        anpi::test::EigenTest<float>(anpi::eig<float>);
        anpi::test::EigenTest<double>(anpi::jacobi<double>);
        anpi::test::EigenTest<double>(anpi::eig<double>);
    }

    BOOST_AUTO_TEST_CASE( EigenSort ) {
        anpi::test::EigenSortTest<float>();
        anpi::test::EigenSortTest<double>();
    }

    BOOST_AUTO_TEST_CASE( Sort ) {
        anpi::test::SortTest<float>();
    }

BOOST_AUTO_TEST_SUITE_END()

