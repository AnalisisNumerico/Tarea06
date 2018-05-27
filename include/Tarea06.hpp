//
// Created by ger534 on 20/05/18.
//

#ifndef PROYECTO3ANALISIS_PACKADAPTER_HPP
#define PROYECTO3ANALISIS_PACKADAPTER_HPP

#include <lapacke.h>
#include "Matrix.hpp"
#include <iostream>
#include <limits>

namespace anpi {

template<typename T>
/**
 * función que construye una matriz cuadrada simetrica
 * @param N tamaño de la matriz NxN
 * @return matriz simetrica cuadrada
 */
anpi::Matrix<T> randomSymmetricSqr(const size_t N){
    anpi::Matrix<T> randomMatrix(N,N);
    std::srand(time(NULL));
    for (int i = 0; i < randomMatrix.cols(); i++) {
        for (int j = 0; j < randomMatrix.rows(); j++) {
            T randomNumber = std::rand() % 101;
            randomMatrix[i][j] = randomNumber;
            randomMatrix[j][i] = randomNumber;
        }
    }
    return randomMatrix;
}
int lp__syev(char jobz, char uplo,  int n  ,double* a, int lda  ,double* w) {
    return LAPACKE_dsyev(LAPACK_ROW_MAJOR,jobz,uplo,  n  ,a,  lda  , w);
}
int lp__syev(char jobz, char uplo,  int n  ,float* a, int lda  ,float* w) {
    return LAPACKE_ssyev(LAPACK_ROW_MAJOR,jobz,uplo,  n  ,a,  lda  , w);
}

template<typename T>
/**
 * función que utiliza la función syev de LAPACKE para resolver el eigensistema
 * @param A matriz simetrica cuadrada
 * @param val vector con los eigenvalores
 * @param E matriz con los eigenvectores
 */
void eig(const anpi::Matrix<T>&A, std::vector<T>&val,anpi::Matrix <T>&E){

    char jobz='V',uplo='U';
    int maxA = int(A.rows());
    int lda=maxA,n=maxA;
    E.allocate(maxA,maxA);
    T a[maxA*maxA] = {};

    int tam = 0;
    for (int i = 0; i < maxA; i++) {
        for (int j = 0; j < maxA; j++) {
            a[tam] = A[i][j];
            tam++;
        }
    }

    T w[A.rows()];

    lp__syev(jobz,uplo,  n  ,a,  lda  , w);

    tam = 0;
    for (int i = 0; i < maxA; i++) {
        for (int j = 0; j < maxA; j++) {
            E[i][j]=a[tam];
            tam++;
        }
    }

    int maxE = int(E.rows());
    for (int i = 0; i < maxE; i++) {
        val.push_back(w[i]);
    }
}

    template<typename T>
    void rotate(anpi::Matrix<T>& A,
                int i,
                int j,
                int k,
                int l,
                T tau,
                T s){

        T g = A[i][j];
        T h = A[k][l];

        A[i][j] = g - s * ( h + g * tau );
        A[k][l] = h + s * ( g - h * tau );

    }

template<typename T>
/**
 * función que utiliza la estrategia de Jacobi vista en clase para resolver el eigensistema
 * @param A matriz simetrica cuadrada
 * @param val vector con los eigenvalores
 * @param E matriz con los eigenvectores
 */
void jacobi(const anpi::Matrix<T>&A,std::vector<T>& val,anpi::Matrix <T>&E){

    anpi::Matrix<T> B;
    B = A;

    anpi::Matrix<T> b(1,A.rows());
    anpi::Matrix<T> z(1,A.rows());

    T tresh,theta,tau,t,sm,s,h,g,c;

    T eps = std::numeric_limits<T>::epsilon();

    E = anpi::Matrix<T>(A.rows(),A.rows());
    for(int i =0; i < int(E.rows()); ++i){
        E[i][i] = T(1.0);
    }

    for(int i =0; i < int(A.rows()); ++i){
        b[0][i]  = B[i][i];
        val.push_back(B[i][i]);
    }

    for(int i = 0; i <= 50; ++i){

        sm = T(0);

        for(int ip = 0; ip < int(A.rows()) -1; ++ip){
            for(int iq = ip+1; iq < int(A.rows()); ++iq){

                sm += (std::abs(B[ip][iq]));

            }
        }

        if(sm == T(0)){


            return;
        }

        if(i < 4){

            tresh = T(0.2) * sm / (T(A.rows())*T(A.rows()));

        }else{

            tresh = T(0);

        }

        for(int ip = 0; ip < int(A.rows())-1; ++ip){
            for(int iq = ip+1; iq < int(A.rows()); ++iq){

                g = T(100.0) * std::abs(B[ip][iq]);

                if(i > 4 && g <= eps * std::abs(val[ip])
                   && g <= eps * std::abs(val[iq])){

                    B[ip][iq] = T(0.0);

                }else if(std::abs(B[ip][iq]) > tresh){

                    h = val[iq] - val[ip];

                    if(g <= eps * std::abs(h)){

                        t = (B[ip][iq]) / h;

                    }else{

                        theta = T(0.5) * h / (B[ip][iq]);
                        t = T(1.0) / (std::abs(theta) + std::sqrt(T(1.0) + theta * theta));
                        if( theta < T(0)) t = -t;

                    }

                    c = T(1.0) / std::sqrt(T(1.0) + t*t);
                    s = t * c;
                    tau = s / (T(1.0) + c);
                    h = t * B[ip][iq];
                    z[0][ip] -= h;
                    z[0][iq] += h;
                    val[ip] -= h;
                    val[iq] += h;
                    B[ip][iq] = T(0);

                    for(int j = 0; j < ip; ++j){

                        rotate(B,j,ip,j,iq,tau,s);
                    }
                    for(int j = ip + 1; j < iq; ++j){
                        rotate(B,ip,j,j,iq,tau,s);
                    }
                    for(int j = iq + 1; j < int(A.rows()); ++j){
                        rotate(B,ip,j,iq,j,tau,s);
                    }
                    for(int j = 0; j < int(A.rows()); ++j){
                        rotate(E,j,ip,j,iq,tau,s);
                    }

                }
            }
        }

        for(int ip = 0; ip < int(A.rows()); ++ip){

            b[0][ip] += z[0][ip];
            val[ip] = b[0][ip];
            z[0][ip] = T(0);

        }
    }



    throw anpi::Exception("Limite de iteraciones maximas");

}
    /**
     * Función  que  ordena  los eigenvalores de forma descendiente
     * y que a su vez, de forma correspondiente, ordena
     * los eigenvectores
     * @param val vector con los eigenvalores
     * @param E matriz con los eigenvectores
     */
    template<typename T>
    void sort(std::vector<T>&val,anpi::Matrix <T> &E){
        for(int i = 0; i < val.size()-1; i++){
            if(val[i+1]-val[i] < 1e-07){
            }
            else{
                T temp = val[i+1];
                val[i+1] = val[i];
                val[i] = temp;
                std::vector<T> tempV (val.size());
                for(int j = 0; j < val.size(); j++){
                    tempV[j] = E[j][i+1];
                }
                for(int j = 0; j < val.size(); j++){
                    E[j][i+1] = E[j][i];
                }
                for(int j = 0; j < val.size(); j++){
                    E[j][i] = tempV[j];
                }
                i=-1;
            }
        }
    }
}


#endif //PROYECTO3ANALISIS_PACKADAPTER_HPP
