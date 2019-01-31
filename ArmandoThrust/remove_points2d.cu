#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/remove.h>
#include <thrust/random.h>

// This example generates random points in the 
// unit square [0,1)x[0,1) and then removes all 
// points where x^2 + y^2 > 1
//
// The x and y coordinates are stored in separate arrays
// and a zip_iterator is used to combine them together

typedef thrust::device_vector<float>::iterator  DVIfloat;
typedef thrust::tuple<DVIfloat, DVIfloat>     TDVIfloat;
typedef thrust::zip_iterator<TDVIfloat>  zip;

//template <typename T>
struct is_outside_circle
{
    template <typename Tuple>
    __host__ __device__
    bool operator()(const Tuple& tuple) const
    {
        // unpack the tuple into x and y coordinates
//        const T x = thrust::get<0>(tuple);
//        const T y = thrust::get<1>(tuple);
        float x = thrust::get<0>(tuple);
        float y = thrust::get<1>(tuple);

        if (x*x + y*y > 1)
            return true;
        else
            return false;
    }
};

int main(void)
{
    const size_t N = 20;

    // generate random points in the unit square on the host
    thrust::default_random_engine rng;
    thrust::uniform_real_distribution<float> u01(0.0f, 1.0f);
    thrust::host_vector<float> x(N);
    thrust::host_vector<float> y(N);
    thrust::device_vector<float> dx;
    thrust::device_vector<float> dy;
    
    for(size_t i = 0; i < N; i++)
    {
        x[i] = u01(rng);
        y[i] = u01(rng);
    }
    
    dx = x; dy = y;

    // print the initial points
    std::cout << std::fixed;
    std::cout << "Generated " << N << " points" << std::endl;
    for(size_t i = 0; i < N; i++)
        std::cout << "(" << x[i] << "," << y[i] << ")" << std::endl;
    std::cout << std::endl;
    
    TDVIfloat tb = thrust::make_tuple(dx.begin(), dy.begin());
    TDVIfloat te = thrust::make_tuple(dx.end(), dy.end());
    zip zb = thrust::make_zip_iterator(tb);
    zip ze = thrust::make_zip_iterator(te);

    // remove points where x^2 + y^2 > 1 and determine new array sizes
    /*
    size_t new_size = thrust::remove_if(thrust::make_zip_iterator(thrust::make_tuple(dx.begin(), dy.begin())),
                                        thrust::make_zip_iterator(thrust::make_tuple(dx.end(), dy.end())),
                                        is_outside_circle<float>())
                      - thrust::make_zip_iterator(thrust::make_tuple(dx.begin(), dy.begin()));
    */
    
//    zip zn = thrust::remove_if(zb, ze, is_outside_circle<float>());
    zip zn = thrust::remove_if(zb, ze, is_outside_circle());
    size_t new_size = zn - zb;
    
    // resize the vectors (note: this does not free any memory)
    dx.resize(new_size);
    dy.resize(new_size);
    
    thrust::copy(dx.begin(), dx.end(), x.begin());
    thrust::copy(dy.begin(), dy.end(), y.begin());
    x.resize(new_size);
    y.resize(new_size);
    
    x.push_back(1.0);

    // print the filtered points
    std::cout << "After stream compaction, " << new_size << " points remain" << std::endl;
    for(size_t i = 0; i < new_size; i++)
        std::cout << "(" << x[i] << "," << y[i] << ")" << std::endl;
    std::cout << "capacity " << x.capacity() << " " << dx.capacity() << std::endl;
    std::cout << "size " << x.size() << " " << dx.size() << std::endl;

    return 0;
}

