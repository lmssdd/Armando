#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/remove.h>
#include <thrust/iterator/zip_iterator.h>

using namespace thrust;

// Make zip_iterator easy to use
    typedef thrust::device_vector< int >::iterator  IntDIter;
    typedef thrust::tuple< IntDIter, IntDIter >     IntDIterTuple2;
    typedef thrust::zip_iterator< IntDIterTuple2 >  ZipDIter;
    
// Make predicate easy to write
typedef thrust::tuple< int, int > IntTuple2;

// Predicate
struct isTuple2Negative
{
    __host__ __device__ bool operator() ( const IntTuple2& tup )
    {
        const int x = thrust::get<0>( tup );
        return ( x < 0 );
    }
};

int main(void) {

// Many vectors
    thrust::device_vector< int > vec0(10,0);
    thrust::device_vector< int > vec1(10,0);

    IntDIterTuple2 tb = thrust::make_tuple( vec0.begin(), vec1.begin() );
    IntDIterTuple2 te = thrust::make_tuple( vec0.end(), vec1.end() );
    ZipDIter zb = thrust::make_zip_iterator(tb);
    ZipDIter ze = thrust::make_zip_iterator(te);
    
    ZipDIter newEnd = ze;
    newEnd = thrust::remove(zb, ze, 1);
/*
// Remove elements in many vectors if element in vec0 is negative
    ZipDIter newEnd = thrust::remove_if(    thrust::make_zip_iterator( thrust::make_tuple( vec0.begin(), vec1.begin() ) ),
                                            thrust::make_zip_iterator( thrust::make_tuple( vec0.end(), vec1.end() ) ),
                                            isTuple2Negative() );
*/
// Erase the removed elements from the vectors
    IntDIterTuple2 endTuple = newEnd.get_iterator_tuple();
    vec0.erase( thrust::get<0>( endTuple ), vec0.end() );
    vec1.erase( thrust::get<1>( endTuple ), vec1.end() );

}
