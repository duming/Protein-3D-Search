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
#include "Cathdata.hpp"

class GaussIntegral
{
    template<class T>
    class SquareArray
    {
    public:
        SquareArray(int size)
        {
            this->maxSize = size;
            base = new T[maxSize*maxSize];
        }
        ~SquareArray()
        {   
            delete[] base;
        }

        inline T& operator () (int row, int colum)
        {
            //if(row>=size || colum>=size || row<0 || colum<0)
            //    throw std::out_of_range("SquareArray index out of range");
            return base[row*size+colum];
        }

        bool resize(int newSize)
        {
            size = newSize;
            if(newSize <= maxSize)
                return false;

            maxSize = newSize*1.5;
            delete [] base;
            base = new T[maxSize* maxSize];
            return true;
        }

        void printArry()
        {
            for(int i=0; i < size; i++)
            {
                for(int j = 0; j < size; j++)
                    std::cout<<base[i*size + j]<<'\t';
                std::cout<<std::endl;
                std::cout<<"#########################"<<std::endl;
            }
        }

    private:
        T* base;
        int size;
        int maxSize;
    };


    public:
    GaussIntegral(int maxlen):unitvector(maxlen),omega(maxlen)
                              ,absomega(maxlen),partsum(maxlen)
                              ,abspartsum(maxlen)
    {
        maxLen = maxlen; 
        currentProtein = NULL;
    }

    ~GaussIntegral()
    {
    }


    void setProtein(std::vector<POINT> * Pptr)
    {   
        int newSize = Pptr->size();
        unitvector.resize(newSize);
        omega.resize(newSize);
        absomega.resize(newSize);
        partsum.resize(newSize);
        if( abspartsum.resize(newSize))
            maxLen = newSize;
        currentProtein = Pptr;
        proteinLen = Pptr->size();
    }

    void test();


    //input: currentProtein
    //output: calculate all vector from point i to j 
    //for every i < j
    void createUnitVector();


    /* spangle: For  vec1, vec2, and, vec3 on the unit 2-sphere 
     *          spangle gives the exterior angle 
     *          spangle(vec1,vec2,vec3) between the geodesic segments 
     *          from vec1 to vec2 and from vec2 to vec3. 
     */
    double spangle(POINT const &op1, POINT const &op2, POINT const &op3);


    /* cerateomega creates half of the mesure of area on the unit 2-sphere
    *             from which a crossing is seen between 
    *             the line segment polyl[i] to  polyl[i+1] and 
    *             the line segment polyl[j] to polyl[j+1].
    *             Cf.  P. R{\o}gen & H.G. Bohr, A new family of global protein
    *             shape descriptors,  Mathematical Biosciences Volume 182,
    *             Issue 2, April 2003, Pages 167-181.
    */
    void createomega(void);

    
    /* createpartsums creats partsum[i][j]=writhe of the segment
     *                from polyl[i] to polyl[j+1]
     *
     *                using a dynamic programming method:
     *                let S(i,k) denote sum( W(j1,j2) ) for all i<=j1<j2<i+k+1
     *
     *                S(i,k) = 
     *                          0                                               if k<2
     *                          W(i,i+k)                                        if k==2
     *                          S(i,k-1) +[ S(i+1,k) -S(i+1,k-1)] +W(i,i+k)     otherwise
     *
     */
    void createpartsum(void);

    /* Calculation of "the mixed sums" equivalent the the 
     * writhe contribution from the inteval [a;b]
     * to the interval [c;d]. 
     * NB: [a;a] is the a'th line segment 
     * - the one from polyl[a] to polyl[a+1].
     */
    //when isabs = 1 calculate the absolute value
    double mixedsum(int a, int b, int c, int d, int isabs = 0);


    /*
     * following are 29 Gauss Integrals
     */
    //first order
    double int12(int isabs);
    //second order
    double int12_34(int isabs1, int isabs2);
    double int13_24(int isabs1, int isabs2);
    double int14_23(int isabs1, int isabs2);
    //third order
    double int12_34_56(void);
    double int12_35_46(void);
    double int12_36_45(void);
   
    double int13_24_56(void);
    double int13_25_46(void);
    double int13_26_45(void);
    
    double int14_23_56(void);
    double int14_25_36(void);
    double int14_26_35(void);
    
    double int15_23_46(void);
    double int15_24_36(void);
    double int15_26_34(void);
    
    double int16_23_45(void);
    double int16_24_35(void);
    double int16_25_34(void);


    //calculate all 29 gauss integrals in the order defined above
    //the second orders with parameters 0,0  0,1  1,0  1,1
    //return all value in the array that gptr pointed to
    void GaussAll(double * gptr);




    private:
    SquareArray<POINT> unitvector;
    SquareArray<double> omega;
    SquareArray<double> absomega;
    SquareArray<double> partsum;
    SquareArray<double> abspartsum;

    std::vector<POINT>* currentProtein;
    int proteinLen;
    int maxLen; 
    //the threshold that is considered as zero 
    static double Epsilon;
};

#endif
