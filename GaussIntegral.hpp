/* Gauss Integral C++ version
 * Author: Ming Du 
 * Date July 28 2016
 *
 *  This code is a modyfied version of 
 *  "Gauss Integrals - A tool for describing protein backbone geometry  
 *  Authors: Peter R{\o}gen & Boris Fain"                
 *  
 *  in order to make it can be used for multithread program
 *
 *
 */
#ifndef GAUSSINTEGRAL_HPP
#define GAUSSINTEGRAL_HPP
#include "Utility.hpp"

class GaussIntegral
{
    template<class T>
    class SquareArray
    {
    public:
        SquareArray(int size)
        {
            this->size = size;
            base = new T[size*size];
        }
        ~SquareArray()
        {   
            delete[] base;
        }

        T& operator () (int row, int colum)
        {
            if(row>=size || colum>=size || row<0 || colum<0)
                throw std::out_of_range("SquareArray index out of range");
            return base[row*size+colum];
        }

        void resize(int newSize)
        {
            if(newSize <= size)
                return;
            delete [] base;
            base = new T[newSize* newSize];
        }
    private:
        T* base;
        int size;
    };


    public:
    GaussIntegral(int maxlen):unitvector(maxlen),omega(maxlen)
                              ,absomega(maxlen),partsum(maxlen)
                              ,abspartsum(maxlen)
    {
    
    }

    ~GaussIntegral()
    {
    }


    void setProtein(std::vector<POINT> * Pptr)
    {
        currentProtein = Pptr;
    }

    void test();


    //input: 
    void createUnitVector();

    private:
    SquareArray<POINT> unitvector;
    SquareArray<double> omega;
    SquareArray<double> absomega;
    SquareArray<double> partsum;
    SquareArray<double> abspartsum;

    std::vector<POINT>* currentProtein;
    
};

#endif
