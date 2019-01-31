#include <iostream>
#include <vector>
#include <thrust/host_vector.h>
using namespace std;
using namespace thrust;

int main()
{
    vector<int> ivec;
    host_vector<int> hvec;
    for ( int i = 0; i < 100; ++i )
    {
        const int before = ivec.capacity();
        const int hbefore = hvec.capacity();
        ivec.push_back( i );
        hvec.push_back( i );
        const int after = ivec.capacity();
        const int hafter = hvec.capacity();
        
        cout << i << " " << ivec.size() << " " << hvec.size() << endl;
        if ( before != after )
            cout << i <<" Capacity: " << before << " " << after << endl;
        if ( hbefore != hafter )
            cout << i <<" hCapacity: " << hbefore << " " << hafter << endl;
    }
    return 0;
}
