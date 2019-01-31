#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>

#include <thrust/remove.h>
#include <thrust/partition.h>

#include <thrust/iterator/zip_iterator.h>

#include <iostream>

// Make zip_iterator easy to use
typedef thrust::device_vector< int >::iterator DVIint;
typedef thrust::tuple< DVIint, DVIint > t2int;
typedef thrust::zip_iterator<t2int> zipIter;

// Predicate
struct check_int : public thrust::unary_function<int>
{
    __host__ __device__
    bool operator()(int i)
    {
        return 0;
    }
};

int main(void)
{
    // initialize all ten integers of a device_vector to 1
    thrust::device_vector<int> D(10, 1);
    thrust::device_vector<int> I(10);

    t2int tb = thrust::make_tuple(D.begin(), I.begin());
    t2int te = thrust::make_tuple(D.end(), I.end());
    zipIter zb = thrust::make_zip_iterator(tb);
    zipIter ze = thrust::make_zip_iterator(te);
    
    // set the first seven elements of a vector to 9
    thrust::fill(D.begin(), D.begin() + 7, 9);
    thrust::sequence(I.begin(), I.end());
    
    // initialize a host_vector with the first five elements of D
    thrust::host_vector<int> H(D.begin(), D.begin() + 5);

    // set the elements of H to 0, 1, 2, 3, ...
    thrust::sequence(H.begin(), H.end());

    // copy all of H back to the beginning of D
    thrust::copy(H.begin(), H.end(), D.begin());

    // print D
    for(int i = 0; i < D.size(); i++)
        std::cout << "D[" << i << "] = " << D[i] << std::endl;
    
    //thrust::device_vector<int>::iterator new_end = thrust::remove(D.begin(), D.end(), 1);
    //DVIint new_end = thrust::remove(D.begin(), D.end(), 1);
    DVIint new_end = thrust::remove(D.begin(), D.end(), check_int());
    std::cout << "new_end = " << *new_end << std::endl;
    D.erase(new_end, D.end());
    
    // print D
    for(int i = 0; i < D.size(); i++)
        std::cout << "D[" << i << "] = " << D[i] << " I[" << i << "] = " << I[i] << std::endl;
    
    return 0;
}
