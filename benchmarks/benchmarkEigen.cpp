//
// Created by ger534 on 25/05/18.
//

//
// Created by ger534 on 22/04/18.
//

/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Pablo Alvarado
 * @date   29.12.2017
 */


#include <boost/test/unit_test.hpp>


#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <chrono>
#include <vector>
//#include "bits/MatrixArithmetic.hpp"

/**
 * Unit tests for the matrix class
 */
#include "benchmarkFramework.hpp"
#include "Matrix.hpp"
#include "Tarea06.hpp"

BOOST_AUTO_TEST_SUITE( Matrix )

    /// Benchmark for eigenvalues and eigenvectors
    template<typename T>
    class Eigen {
    protected:
        /// Maximum allowed size for the square matrices
        const size_t _maxSize;

        /// A large matrix holding
        anpi::Matrix<T> _data;

        /// State of the benchmarked evaluation
        anpi::Matrix<T> _a;

        anpi::Matrix<T> _E;

        std::vector<T> val;


    public:
        /// Construct
        Eigen(const size_t maxSize)
                : _maxSize(maxSize),_data(maxSize,maxSize,anpi::DoNotInitialize) {
        }

        /// Prepare the evaluation of given size
        void prepare(const size_t size) {
            assert (size<=this->_maxSize);
            this->_data=anpi::randomSymmetricSqr<T>(_maxSize);
            this->_a=anpi::randomSymmetricSqr<T>(size);
        }
    };

/// Provide the evaluation method for jacobi
    template<typename T>
    class Jacobi : public Eigen<T> {
    public:
        /// Constructor
        Jacobi(const size_t n) : Eigen<T>(n) { }

        // Evaluate add in-place
        inline void eval() {
            anpi::jacobi(this->_a,this->val,this->_E);
        }
    };

    /// Provide the evaluation method for eig using lapack
    template<typename T>
    class LaPack : public Eigen<T> {
    public:
        /// Constructor
        LaPack(const size_t n) : Eigen<T>(n) { }

        // Evaluate add on-copy
        inline void eval() {
            anpi::eig(this->_a,this->val,this->_E);
        }
    };


    BOOST_AUTO_TEST_CASE( JACOBI ) {

        std::vector<size_t> sizes = { 1,2,3,4,5,6,7,8,9,10 };
        const size_t n=sizes.back();
        const size_t repetitions=100;
        std::vector<anpi::benchmark::measurement> times;

        {
            Jacobi<float> baip(n);
            ANPI_BENCHMARK(sizes,repetitions,times,baip);
            ::anpi::benchmark::write("Jacobi presicion simple",times);
            ::anpi::benchmark::plotRange(times,"Jacobi presicion simple","b");
        }

        {
            Jacobi<double> baip(n);
            ANPI_BENCHMARK(sizes,repetitions,times,baip);
            ::anpi::benchmark::write("Jacobi presicion doble",times);
            ::anpi::benchmark::plotRange(times,"Jacobi presicion doble","r");
        }

        ::anpi::benchmark::show();
    }

    BOOST_AUTO_TEST_CASE( LAPACK ) {

        std::vector<size_t> sizes = {  1,2,3,4,5,6,7,8,9,10 };
        const size_t n=sizes.back();
        const size_t repetitions=100;
        std::vector<anpi::benchmark::measurement> times;

        {
            LaPack<float> baip(n);
            ANPI_BENCHMARK(sizes,repetitions,times,baip);
            ::anpi::benchmark::write("LaPack presicion simple",times);
            ::anpi::benchmark::plotRange(times,"LaPack presicion simple","b");
        }

        {
            LaPack<double> baip(n);
            ANPI_BENCHMARK(sizes,repetitions,times,baip);
            ::anpi::benchmark::write("LaPack presicion doble",times);
            ::anpi::benchmark::plotRange(times,"LaPack presicion doble","r");
        }

        ::anpi::benchmark::show();
    }
BOOST_AUTO_TEST_SUITE_END()