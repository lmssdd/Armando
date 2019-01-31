#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/remove.h>
#include <iostream>
#include <stdio.h>

// helper routine
template <typename String, typename Vector>
void print(const String& s, const Vector& v)
{
  std::cout << s << " [";
  for(size_t i = 0; i < v.size(); i++)
    std::cout << " " << v[i];
  std::cout << " ]\n";
}

// this functor returns true if the argument is odd, and false otherwise
template <typename T>
struct is_odd : public thrust::unary_function<T,bool>
{
    __host__ __device__
    bool operator()(T x)
    {
        return x % 2;
    }
};

struct is_ok
{
    __host__ __device__
    bool operator()(int x)
    {
        return x != 0;
    }
};

template <typename Iterator>
void print_range(const std::string& name, Iterator first, Iterator last)
{
    typedef typename std::iterator_traits<Iterator>::value_type T;

    std::cout << name << ": ";
    thrust::copy(first, last, std::ostream_iterator<T>(std::cout, " "));  
    std::cout << "\n";
}

__global__ void checkKeys(const int pn, const int* values, int* keys) {

    int ip;
	
    ip = threadIdx.x + blockDim.x * blockIdx.x;
    
    if (ip < pn) {
        if (values[ip] > 4) keys[ip] = 1;
        }
}

int main(void)
{
    // input size
    size_t N = 10;

    // define some types
    typedef thrust::device_vector<int> Vector;
    typedef Vector::iterator           Iterator;

    // allocate storage for array
    Vector values(N);

    // initialize array to [0, 1, 2, ... ]
    thrust::sequence(values.begin(), values.end());
    
    print_range("values", values.begin(), values.end());

    // another approach is to count the number of values that will 
    // be copied, and allocate an array of the right size
    size_t N_odd = thrust::count_if(values.begin(), values.end(), is_odd<int>());
    
    Vector small_output(N_odd);
    
    thrust::copy_if(values.begin(), values.end(), small_output.begin(), is_odd<int>());
    
    print_range("small_output", small_output.begin(), small_output.end());

	// test
    printf("%d\n", values.size());
    Vector keys(N);
    thrust::fill(keys.begin(), keys.end(), 0);
    
    int blocks = (values.size() + 256 - 1) / 256;

    checkKeys<<< blocks, 256 >>>
		(values.size(), 
		thrust::raw_pointer_cast(&values[0]), 
		thrust::raw_pointer_cast(&keys[0]));
	
    size_t nplus = thrust::count(keys.begin(), keys.end(), 1);
    printf("%d\n", nplus);
    
    print("values ", values);
    values.resize(values.size() +nplus);
    print("values ", values);

//    thrust::copy_if(values.begin(), values.begin() +10, values.begin() +10, is_ok());
    thrust::copy_if(values.begin(), values.begin() +10, keys.begin(), values.begin() +10, is_ok());
    print("values ", values);
    return 0;
}

